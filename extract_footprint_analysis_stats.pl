#!/usr/bin/perl
#
# Extraction of statistics for footprint analysis. 
# This script uses the Perl rearrangement oo-framework to extract features
# of footprints and rearrangement clusters for downstream analysis. 
#
# The script assumes that the input rearrangements are properly filtered.
# No further filtering is performed on rearrangements, except that shard-
# bypassing rearrangements are removed as they contain redundant information.
#

########## LICENCE ##########
# Copyright (c) 2015 Genome Research Ltd.
#
# Author: Yilong Li <yl3@sanger.ac.uk>
#
# This file is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
########## LICENCE ##########

use strict;
use warnings FATAL => 'all';

use Getopt::Long;
use Pod::Usage qw(pod2usage);
use autodie qw(:all);
use Carp 'verbose';
use Scalar::Util qw(blessed);
use List::Util qw(sum min max);
$SIG { __DIE__ } = sub { Carp::confess(@_) };

use File::Basename;
use lib dirname($0);
use CopyNumberSegmentGenome;
use CopyNumberSegmentArray;
use CopyNumberSegment;
use RearrangementGroup;
use Footprint;

#
# Handle options
#
pod2usage(-exitval => 1, -output => \*STDERR) if @ARGV < 2;

# Set default options
my %opts = (
    # Basic options
    sample_name => "",
    acf    => "",
    ploidy => "",
    verbose => 0,
    help => 0,

    # Advanced options for extracting rearrangement patterns
    max_balanced_rg_dist => 500,
    max_balanced_rg_overlap => 100,  # Not exposed as an option yet
    max_foldback_distance => 50000,
    min_seg_size_for_small_del_or_td => 5000,  # Minimum segment size for determining CN
    min_cn_change => 0.3,  # Minimum CN change for a copy number or rg breakpoint
    min_cn_bkpt_seg_size => 50000,  # Segment bounded by copy number segmentation breakpoints only and smaller than this size will be merged with neighbouring segments. 
    keep_small_dels_and_tds => 1,
    max_shard_length => 10000,
    shard_bypassing_slop => 200,
);
GetOptions(
    "acf=f" => \&opt_handler,
    "ploidy=f" => \&opt_handler,
    "verbose:1" => \$opts{verbose},
    "help" => \$opts{help},
    "max_balanced_rg_dist=i" => \&opt_handler,
    "max_foldback_distance=i" => \&opt_handler,
    "min_seg_size_for_cn=i" => \&opt_handler,
    "min_cn_change=f" => \&opt_handler,
    "min_cn_bkpt_seg_size=i" => \&opt_handler,
    "shard_bypassing_slop=i" => \&opt_handler,
    "max_local_two_jump_distance=i" => \$opts{max_local_two_jump_distance},
    "sample_name=s" => \$opts{sample_name},

    # Output files
    "sv_classification=s" => \$opts{sv_classification},

    # Files for SV clusters info
    "sv_clusters_file=s" => \$opts{sv_clusters_file}
);

my $rgs_file = shift;
my $cn_file  = shift;
pod2usage(-msg => "File '$rgs_file' doesn't exist!", -exitval => 2) if !-e $rgs_file;
pod2usage(-msg => "File '$cn_file' doesn't exist!", -exitval => 2) if !-e $cn_file;
pod2usage(-msg => "Option -sample_name STRING is required!", -exitval => 2) if !$opts{sample_name};
pod2usage(-msg => "Option -acf FLOAT is required!", -exitval => 2) if !$opts{acf};
pod2usage(-msg => "Option -ploidy FLOAT is needed!", -exitval => 2) if !$opts{ploidy};
pod2usage(-msg => "Option -sv_clusters_file STRING is needed!", -exitval => 2) if !$opts{sv_clusters_file};
pod2usage(-msg => "Option -sv_classification FILE is needed!", -exitval => 2) if !$opts{sv_classification};
pod2usage() if $opts{help};


#
# Actual commands for the code
#

sub stderr_with_date {
    print STDERR "[" . `echo -n \`date\`` . "] " . $_[0];
}

# First read in all rearrangement and copy number data into a genome. 
stderr_with_date("Reading data in...\n") if $opts{verbose};
my $genome_of_cn_segs = CopyNumberSegmentGenome->new_from_BEDPE_and_CN_BED($rgs_file, $cn_file, %opts);

stderr_with_date("Sorting segments on each chromosome...\n") if $opts{verbose};
$genome_of_cn_segs->sort_segments;

stderr_with_date("Bypassing shards...\n") if $opts{verbose};
$genome_of_cn_segs->preprocess_rearrangement_and_copy_number_data(%opts, "no_shard_bypassing" => 1, "no_rg_filtering" => 1);

stderr_with_date("Reading in rearrangement clustering information...\n") if $opts{verbose};
$genome_of_cn_segs->read_events_from_clustering_file($opts{sv_clusters_file});

stderr_with_date("Cleaning rearrangement clusters...\n") if $opts{verbose};
$genome_of_cn_segs->normalise_sv_clustering(%opts);

stderr_with_date("Printing out sv classification to $opts{sv_classification}\n") if $opts{verbose};
open NAIVE_CLASSIFICATION, ">$opts{sv_classification}" or die $!;
select NAIVE_CLASSIFICATION;
$genome_of_cn_segs->print_footprint_cluster_classifications(%opts);
close NAIVE_CLASSIFICATION;

stderr_with_date("Done.\n") if $opts{verbose};
exit;


sub opt_handler {
    my ($opt_name, $opt_value) = @_;
    if ($opt_name eq "acf") {
        if ($opt_value < 0 || $opt_value > 100) {
            pod2usage("-acf value must be between 0 and 100!\n");
        }
        if ($opt_value > 1) { $opt_value /= 100; }
        $opts{acf} = $opt_value;
    }
    elsif ($opt_name eq "ploidy") {
        if ($opt_value <= 0) {
            pod2usage("-ploidy must be larger than 0!\n");
        }
        $opts{ploidy} = $opt_value;
    }
    elsif ($opt_name eq "max_balanced_rg_dist") {
        # if ($opt_value <= 0) {
        #     pod2usage("-max_balanced_rg_dist must be a non-negative integer!\n");
        # }
        $opts{max_balanced_rg_dist} = $opt_value;
    }
    elsif ($opt_name eq "max_foldback_distance") {
        if ($opt_value <= 0) {
            pod2usage("-max_foldback_distance must be a non-negative integer!\n");
        }
        $opts{max_foldback_distance} = $opt_value;
    }
    elsif ($opt_name eq "min_seg_size_for_cn") {
        if ($opt_value <= 0) {
            pod2usage("-min_seg_size_for_cn must be a non-negative integer!\n");
        }
        $opts{min_seg_size_for_cn} = $opt_value;
    }
    elsif ($opt_name eq "min_cn_change") {
        if ($opt_value <= 0) {
            pod2usage("-min_cn_change must be non-negative!\n");
        }
        $opts{min_cn_change} = $opt_value;
    }
    elsif ($opt_name eq "min_cn_bkpt_seg_size") {
        if ($opt_value <= 0) {
            pod2usage("-min_cn_bkpt_seg_size must be non-negative!\n");
        }
        $opts{min_cn_bkpt_seg_size} = $opt_value;
    }
    elsif ($opt_name eq "shard_bypassing_slop") {
        if ($opt_value <= 0) {
            pod2usage("-shard_bypassing_slop must be non-negative!\n");
        }
        $opts{shard_bypassing_slop} = $opt_value;
    }
    else {
        die;
    }
}


__END__


=head1 NAME


extract_footprint_analysis_stats.pl

This script extracts various information for statistical analysis of rearrangement
footprints and clusters. 

This script is not designed for general usage. 

=head1 SYNOPSIS

extract_footprint_analysis_stats.pl [OPTIONS] REARRANGEMENTS.BEDPE CN_SEGMENTS.BED

REARRANGEMENTS.BEDPE is a BEDPE file with rearrangement ID in column 7
and strands in columns 9 and 10. In addition, columns 15-18 must have
format /\d+ \(\d+\)/ indicating the positions of clipped reads at
rearrangement breakpoints plus number of reads supporting the
clipping positions. 

CN_SEGMENTS.BEDGRAPH is a BEDGRAPH file of copy number segments.
The fourth column corresponds to the absolute copy number of the segment.
The fifth and sixth columns correspond to copy number breakpoint type
(rearrangement vs. copy number segmentation breakpoint). The seventh
column is the number of windows in the current segment.

  Basic options:
    -help                       Print this message
    -acf FLOAT                  Aberrant cell fraction. Required
    -ploidy FLOAT               Tumour ploidy. Required
    -sample_name STRING         Sample name. Required
    -sv_clusters_file STRING    File of rearrangement clustering and footprint
                                results. 
    -verbose                    Print debugging messages

  Output file names - all required
    -arm_cn_stats_file STRING   Output file for arm-level CN data
    -size_1_footprints_file     File name for the size one footprints stats output.  
                        STRING 
    -size_2_footprints_file     File name for the size two footprints stats output.  
                        STRING 
    -size_3_footprints_file     File name for the size three footprints stats output.
                        STRING
    -size_4_footprints_file     File name for the size four footprints stats output.
                        STRING
    -size_5_footprints_file     File name for the size five footprints stats output.
                        StRING
    -misc_stats_file STRING     File name for misc footprint statistics. 

  Advanced options - Defaults were used for pancan
    -max_balanced_rg_dist INT       Maximum distance at which reciprocal
                                    rearrangements can still be considered balanced
                                    [1000]
    -max_balanced_rg_overlap INT    Maximum overlap between balanced breakpoints
                                    with some overlapping sequence. [100] 
    -max_foldback_distance INT      Maximum distance for fold-back type
                                    rearrangements to be considered as purely
                                    fold-back. [5000]
    -max_local_two_jump_distance    Maximum distance between two footprints of a
                               INT  two-jump event for it to still be considered
                                    local.
