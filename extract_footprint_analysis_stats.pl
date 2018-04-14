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
    # max_local_two_jump_distance => 1e6,  # Maximum distance between two footprints for a two-jump to still be called a two-jump

    # Options for library matching
    diploid_matching => 0,
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
    "arm_cn_stats_file=s" => \$opts{arm_cn_stats_file},
    "size_1_footprints_file=s" => \$opts{size_1_footprints_file},
    "size_2_footprints_file=s" => \$opts{size_2_footprints_file},
    "size_3_footprints_file=s" => \$opts{size_3_footprints_file},
    "size_4_footprints_file=s" => \$opts{size_4_footprints_file},
    "size_5_footprints_file=s" => \$opts{size_5_footprints_file},
    "size_6_footprints_file=s" => \$opts{size_6_footprints_file},
    "local_3_jumps_file=s" => \$opts{local_3_jumps_file},
    "local_4_jumps_file=s" => \$opts{local_4_jumps_file},
    "misc_stats_file=s" => \$opts{misc_stats_file},

    "diploid_matching" => \&opt_handler,

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
pod2usage(-msg => "Option -arm_cn_stats_file STRING is needed!", -exitval => 2) if !$opts{arm_cn_stats_file};
pod2usage(-msg => "Option -size_1_footprints_file STRING is needed!", -exitval => 2) if !$opts{size_1_footprints_file};
pod2usage(-msg => "Option -size_2_footprints_file STRING is needed!", -exitval => 2) if !$opts{size_2_footprints_file};
pod2usage(-msg => "Option -misc_stats_file FILE is needed!", -exitval => 2) if !$opts{misc_stats_file};
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

stderr_with_date("Extracting chromosome arm CN statistics...\n") if $opts{verbose};
open OUT, "> $opts{arm_cn_stats_file}" or die $!;
select OUT;
print_chr_arm_cn_stats(%opts);
close OUT;

stderr_with_date("Extracting statistics for size one footprints...\n") if $opts{verbose};
open OUT, "> $opts{size_1_footprints_file}" or die $!;
select OUT;
print_size_one_footprint_stats(%opts);
close OUT;

stderr_with_date("Extracting statistics for size two footprints...\n") if $opts{verbose};
open OUT, "> $opts{size_2_footprints_file}" or die $!;
select OUT;
# print_size_two_footprint_stats(%opts);
print_general_footprint_stats(%opts, size => 2);
close OUT;

stderr_with_date("Extracting statistics for size three footprints...\n") if $opts{verbose};
open OUT, "> $opts{size_3_footprints_file}" or die $!;
select OUT;
print_general_footprint_stats(%opts, size => 3);
close OUT;

stderr_with_date("Extracting statistics for size four footprints...\n") if $opts{verbose};
open OUT, "> $opts{size_4_footprints_file}" or die $!;
select OUT;
print_general_footprint_stats(%opts, size => 4);
close OUT;

stderr_with_date("Extracting statistics for size five footprints...\n") if $opts{verbose};
open OUT, "> $opts{size_5_footprints_file}" or die $!;
select OUT;
print_general_footprint_stats(%opts, size => 5);
close OUT;

stderr_with_date("Extracting statistics for size six footprints...\n") if $opts{verbose};
open OUT, "> $opts{size_6_footprints_file}" or die $!;
select OUT;
print_general_footprint_stats(%opts, size => 6);
close OUT;

stderr_with_date("Printing out misc footprint motif stats\n") if $opts{verbose};
open OUT, "> $opts{misc_stats_file}";
select OUT;
print_misc_stats($genome_of_cn_segs);
close OUT;

stderr_with_date("Printing out local 3-jumps\n") if $opts{verbose};
open OUT, "> $opts{local_3_jumps_file}";
select OUT;
print_local_n_jump_footprints($genome_of_cn_segs, 3, %opts);
close OUT;

stderr_with_date("Printing out local 4-jumps\n") if $opts{verbose};
open OUT, "> $opts{local_4_jumps_file}";
select OUT;
print_local_n_jump_footprints($genome_of_cn_segs, 4, %opts);
close OUT;

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
    elsif ($opt_name eq "diploid_matching") {
        print STDERR "Option -diploid_matching is not implemented yet. Exiting.\n";
        exit;
    }
    else {
        die;
    }
}

sub print_chr_arm_cn_stats {
    # This function traverses the global variable $genome_of_cn_segs and
    # prints the arm and telomere level copy numbers to currently open
    # file handle.
    print join("\t", "sample", "chr", "arm", "centromere_cn", "telomere_cn") . "\n";

    my %params = @_;

    for my $chr ($genome_of_cn_segs->chrs_array) {
        print join("\t", $params{sample_name}, $chr->name, "p", ($chr->p_centromere_cn or "NA"), ($chr->p_arm_telomere_cn or "NA")) . "\n";
        print join("\t", $params{sample_name}, $chr->name, "q", ($chr->q_centromere_cn or "NA"), ($chr->q_arm_telomere_cn or "NA")) . "\n";
    }
}

sub print_size_one_footprint_stats {
    # This function traverses the global variable $genome_of_cn_segs, extracts
    # all footprints with one breakpoints and outputs the relevant stats.
    # The stats are printed to the currently open file handle. 
    print join(
        "\t",
        "sample",
        "cluster_id",
        "cluster_size",
        "footprint_id",
        "chr",
        "pos",
        "distance_to_mate",
        "orientation",
        "rg_id",
        "rg_type",
        "rg_side_cn",
        "non_rg_side_cn",
        "rg_side_arm_cn",
        "non_rg_side_arm_cn",
        "arm_has_stable_cn",
        "p_arm_cn",
        "q_arm_cn",
    ) . "\n";

    my %params = @_;
    my $bkpt;
    for my $cluster (values %{$genome_of_cn_segs->{rg_clusters}}) {
        for my $footprint ($cluster->footprints_array) {
            if ($footprint->size != 1) {
                next;
            }

            ($bkpt) = $footprint->rg_ends_array;
            print join(
                "\t",
                $params{sample_name},
                $cluster->id,
                $cluster->size,
                $footprint->id,
                $bkpt->chr_name,
                $bkpt->pos,
                $bkpt->distance_to_mate,
                $bkpt->dir,
                $bkpt->id,
                $bkpt->rg->rg_type_s(%params),
                ($bkpt->rg_side_cn(%params) or "NA"),
                ($bkpt->non_rg_side_cn(%params) or "NA"),
                ($bkpt->rg_side_arm_cn or "NA"),
                ($bkpt->non_rg_side_arm_cn or "NA"),
                $bkpt->chr->has_stable_arm_cn($bkpt->is_on_p_arm ? 'p' : 'q'),
                ($footprint->chrom->p_arm_telomere_cn or "NA"),
                ($footprint->chrom->q_arm_telomere_cn or "NA"),
            ) . "\n";
        }
    }
}

sub print_size_two_footprint_stats {
    # This function traverses the global variable $genome_of_cn_segs, extracts
    # all footprints with two breakpoints and outputs the relevant stats.
    # The stats are printed to the currently open file handle. 
    print join(
        "\t",
        "sample",
        "cluster_id",
        "cluster_size",
        "footprint_id",
        "chr",
        "start",
        "end",
        "rg_ids",
        "low_end_relative_cn",
        "high_end_relative_cn",
        "type",  # Primary type
        "detailed_type",  # Additional info of the footprint type
        "distance",
        "orientations",
        "upstream_dist",
        "downstream_dist",
    ) . "\n";

    my %params = @_;
    my ($bkpt_1, $bkpt_2, $us_end, $ds_end);
    for my $cluster (values %{$genome_of_cn_segs->{rg_clusters}}) {
        for my $footprint ($cluster->footprints_array) {
            if ($footprint->size != 2) {
                next;
            }

            # Flip the two rearrangement breakpoints to right order
            ($bkpt_1, $bkpt_2) = $footprint->rg_ends_array;
            if (
                $bkpt_1->segment->is_bal_rg_overlap(%params) &&
                $bkpt_1->balanced_bkpt_partner_rg_end(%params) == $bkpt_2 &&
                $bkpt_1->is_rev
            ) {
                ($bkpt_1, $bkpt_2) = ($bkpt_2, $bkpt_1);
            }
            elsif ($bkpt_1->pos > $bkpt_2->pos) {
                ($bkpt_1, $bkpt_2) = ($bkpt_2, $bkpt_1);
            }

            $us_end = $bkpt_1->closest_us_rg_end(within => 3e9, exclude_rgs => []);
            $us_end = defined($us_end) ? $bkpt_1->pos - $us_end->pos + 1 : "Inf";
            $ds_end = $bkpt_2->closest_ds_rg_end(within => -1, exclude_rgs => []);
            $ds_end = defined($ds_end) ? $ds_end->pos - $bkpt_2->pos + 1 : "Inf";

            print join(
                "\t",
                $params{sample_name},
                $cluster->id,
                $cluster->size,
                $footprint->id,
                $footprint->chr_name,
                $footprint->first_pos,
                $footprint->last_pos,
                join(",", map { $_->id } ($bkpt_1, $bkpt_2)),  # rg IDs
                ($bkpt_1->cn_relative_to_arm(%params) or "NA"),
                ($bkpt_2->cn_relative_to_arm(%params) or "NA"),
                $footprint->type(%params),
                $footprint->detailed_type(%params),
                1 + abs($bkpt_1->pos - $bkpt_2->pos),
                join("", map { $_->dir } ($bkpt_1, $bkpt_2)),
                $us_end,
                $ds_end,
            ) . "\n";
        }
    }
}

sub print_footprint_details {
    # Erm, print details of a footprint

    my $footprint = shift;
    my $cluster = $footprint->cluster;
    my %params = @_;

    my $rg_pattern;
    my($us_bkpt, $ds_bkpt, $us_dist, $ds_dist, $bkpt_1, $bkpt_2, @rg_ends,
            %rg_seen, @rg_cn_changes, $us_segment, $ds_segment);

    @rg_ends = $footprint->sorted_rg_ends_array(%params);
    if ($rg_ends[0]->is_fwd and !($rg_ends[0]->is_bal_rg_overlap(%params))) {
        $us_segment = $rg_ends[0]->segment;
    }
    else {
        $us_segment = $rg_ends[0]->segment->prev_seg;
    }
    if ($rg_ends[-1]->is_rev and !($rg_ends[-1]->is_bal_rg_overlap(%params))) {
        $ds_segment = $rg_ends[-1]->segment;
    }
    else {
        $ds_segment = $rg_ends[-1]->segment->next_seg;
    }

    $us_bkpt = $bkpt_1 = $rg_ends[0];
    do {
        $us_bkpt = $us_bkpt->closest_us_rg_end(within => 3e9, exclude_rgs => []);
    } while (defined($us_bkpt) and $us_bkpt->footprint == $bkpt_1->footprint);
    $us_dist = defined($us_bkpt) ? $bkpt_1->pos - $us_bkpt->pos + 1 : "Inf";

    $ds_bkpt = $bkpt_2 = $rg_ends[-1];
    do {
        $ds_bkpt = $ds_bkpt->closest_ds_rg_end(within => 3e9, exclude_rgs => []);
    } while (defined($ds_bkpt) and $ds_bkpt->footprint == $bkpt_2->footprint);
    $ds_dist = defined($ds_bkpt) ? $ds_bkpt->pos - $bkpt_2->pos + 1 : "Inf";

    # Print the copy number changes of each rearrangement
    %rg_seen = ();
    @rg_cn_changes = ();
    for my $bkpt (@rg_ends) {
        next if $rg_seen{$bkpt->rg->id};
        $rg_seen{$bkpt->rg->id} = 1;
        if ($bkpt->rg->weighted_avg_cn_change_across_rg(%params)) {
            push @rg_cn_changes, $bkpt->rg->weighted_avg_cn_change_across_rg(%params);
        }
        else {
            push @rg_cn_changes, "NA";
        }
    }

    # Footprint detailed type
    my $detailed_type = $footprint->detailed_type(%params);
    if (defined($detailed_type) and $detailed_type eq "shard_cycle:2") {
        if ($cluster->cluster_classification(%params) eq "TD_after_unbal_transloc") {
            $detailed_type = "TD_after_unbal_transloc";
        }
    }

    $rg_pattern = $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1);
    if (defined($rg_pattern)) {
        print join(
            "\t",
            $params{sample_name},
            $cluster->id,
            $cluster->size,
            $footprint->id,
            # scalar($footprint->rgs_array),
            $cluster->n_rgs,
            $footprint->chr_name,
            $footprint->first_pos,
            $footprint->last_pos,
            join(",", keys($footprint->rgs)),  # rg IDs
            $us_dist,
            $ds_dist,
            $detailed_type,
            $rg_pattern->{rg_string},
            $rg_pattern->{is_inverted},
            join(",", @{$rg_pattern->{distances}}),
            join(",", map( { defined($_) ? sprintf("%.3f", $_) : "NA" } @{$rg_pattern->{relative_cns}})),
            join(",", map( { $_ eq "NA" ? "NA" : sprintf("%.3f", $_) } @{$rg_pattern->{cn_changes}})),
            join(",", map( { $_ eq "NA" ? "NA" : sprintf("%.3f", $_) } @{$rg_pattern->{cn_change_variances}})),
            join(",", @rg_cn_changes),
            ($footprint->chrom->p_arm_telomere_cn or "NA"),
            ($footprint->chrom->q_arm_telomere_cn or "NA"),
            ($us_segment->cn or "NA") . "," . ($ds_segment->cn or "NA"),
            $us_segment->length . "," . $ds_segment->length,
        ) . "\n";
    }
    else {
        print join(
            "\t",
            $params{sample_name},
            $cluster->id,
            $cluster->size,
            $footprint->id,
            scalar($footprint->rgs_array),
            $footprint->chr_name,
            $footprint->first_pos,
            $footprint->last_pos,
            join(",", keys($footprint->rgs)),  # rg IDs
            $us_dist,
            $ds_dist,
            $footprint->detailed_type(%params),
            ("NA")x6,
            join(",", @rg_cn_changes),
            ($footprint->chrom->p_arm_telomere_cn or "NA"),
            ($footprint->chrom->q_arm_telomere_cn or "NA"),
            ($us_segment->cn or "NA") . "," . ($ds_segment->cn or "NA"),
            $us_segment->length . "," . $ds_segment->length,
        ) . "\n";
    }
}

sub print_footprint_details_header {
    print join(
        "\t",
        "sample",
        "cluster_id",
        "cluster_size",
        "footprint_id",
        "n_rgs",
        "chr",
        "start",
        "end",
        "rg_ids",
        "distances_to_closest_us_rg",
        "distances_to_closest_ds_rg",
        "footprint_type",
        "rg_pattern",
        "pattern_is_inverted",
        "distances",
        "relative_cns",
        "cn_changes",
        "cn_change_variances",
        "rg_cn_changes",
        "p_arm_cn",
        "q_arm_cn",
        "neighbour_segment_cns",
        "neighbour_segment_lengths",
    ) . "\n";
}

sub print_general_footprint_stats {
    # This function traverses the global variable $genome_of_cn_segs, extracts
    # all footprints with three breakpoints and outputs the relevant stats.
    # The stats are printed to the currently open file handle. 

    print_footprint_details_header();
    my %params = @_;
    my $rg_pattern;
    my($us_bkpt, $ds_bkpt, $us_dist, $ds_dist, $bkpt_1, $bkpt_2, @rg_ends,
            %rg_seen, @rg_cn_changes, $us_segment, $ds_segment);

    for my $cluster (values %{$genome_of_cn_segs->{rg_clusters}}) {
        for my $footprint ($cluster->footprints_array) {
            if ($footprint->size != $params{size}) {
                next;
            }

            print_footprint_details($footprint, %params);
        }
    }
}


sub print_misc_stats {
    # This function prints misc stats into a single file, with rows of each
    # data type prepended with a tag. 
    my $genome = shift;

    print_bkpt_single_shard_stats(%opts, tag => "TRANSLOC_SINGLE_SHARD_STATS");
    find_inversion_nested_in_tds(%opts, tag => "INVERSION_WITHIN_TD_STATS");
    cn_shift_of_local_two_jump_chrs(%opts, tag => "LOCAL_TWO_JUMP_CHR_CN_SHIFT");
    cn_shift_of_size_two_footprint_chrs(%opts, tag => "TWO_BKPT_FOOTPRINT_CHR_CN_SHIFT");
    size_3_footprint_partners($genome, %opts, tag => "SIZE_3_FOOTPRINT_PARTNERS");
    four_bkpt_footprint_partners($genome, %opts, tag => "FOUR_BKPT_FOOTPRINT_PARTNERS");
    print_bal_olap_bkpt_footprints($genome, %opts, tag => "BAL_BKPT_OLAP");
    shard_transloc_footprint_partners($genome, %opts, tag => "SHARD_TRANSLOC_PARTNERS");
    print_chromoplexy_events($genome, %opts, tag => "CHROMOPLEXY");
}

sub print_bkpt_single_shard_stats {
    my %params = @_;

    # Print header
    print join(
        "\t",
        $params{tag},
        "sample_name",
        "footprint_pattern",
        "cluster_id",
        "footprint_id",
        "chr",
        "shard_start",
        "shard_end",
        "dir",
        "shard_size",
        "distance_to_unbal_bkpt",
        "cluster_size",
        "footprint_size",
        "low_end_mate",
        "high_end_mate",
        "low_end_cluster_size",
        "high_end_cluster_size",
        "cn_changes",
        "relative_cns",
    ) . "\n";

    my ($footprint, $cluster, $rg_pattern, @sorted_bkpts, $start, $end, $dist_to_bkpt, $low_end_partner, $high_end_partner, $dist);
    my @footprints_of_interest = ();
    for my $cluster (values %{$genome_of_cn_segs->{rg_clusters}}) {
        for my $footprint ($cluster->footprints_array) {
            next if $footprint->size == 0;  # Peeled off footprints
            $rg_pattern = $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1);
            if (
                defined($rg_pattern) and (
                    $rg_pattern->{rg_string} eq "1-/1+" or     # Single shard
                    $rg_pattern->{rg_string} eq "0+/2-/2+" or  # Inserted shard
                    $rg_pattern->{rg_string} eq "0+,2+/2-" or  # Inversion shard
                    $rg_pattern->{rg_string} eq "0+,2-/2+" or  # Deletion shard
                    $rg_pattern->{rg_string} eq "1-,2-/2+" or  # Fold-back translocation
                    $rg_pattern->{rg_string} eq "1-/1+/2+"
                )
            ) {
                push @footprints_of_interest, $footprint;
            }
        }
    }

    my @out_fields;
    for $footprint (@footprints_of_interest) {
        $rg_pattern = $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1);

        # Get the start and end of the inserted shard
        @sorted_bkpts = $footprint->sorted_rg_ends_array(%params);
        if ($rg_pattern->{rg_string} eq "1-/1+") {
            ($start, $end) = map { $_->pos } @sorted_bkpts;
            $dist = "NA";
            $low_end_partner = $sorted_bkpts[0]->mate;
            $high_end_partner = $sorted_bkpts[1]->mate;
        }
        elsif (!$rg_pattern->{is_inverted}) {
            if ($rg_pattern->{rg_string} eq "1-/1+/2+") {
                ($start, $end) = ($sorted_bkpts[0]->pos, $sorted_bkpts[1]->pos);
                $dist = max(1, $sorted_bkpts[2]->pos - $sorted_bkpts[1]->pos);
                $low_end_partner = $sorted_bkpts[0]->mate;
                $high_end_partner = $sorted_bkpts[1]->mate;
            }
            else {
                # This applies to all of:
                # - 0+/2-/2+
                # - 0+,2+/2-
                # - 0+,2-/2+
                # - 1-,2-/2+
                ($start, $end) = ($sorted_bkpts[1]->pos, $sorted_bkpts[2]->pos);
                $dist = max(1, $sorted_bkpts[1]->pos - $sorted_bkpts[0]->pos);
                $low_end_partner = $sorted_bkpts[1]->mate;
                $high_end_partner = $sorted_bkpts[2]->mate;
            }
        }
        else {
            if ($rg_pattern->{rg_string} eq "1-/1+/2+") {
                ($start, $end) = ($sorted_bkpts[1]->pos, $sorted_bkpts[2]->pos);
                $dist = max(1, $sorted_bkpts[1]->pos - $sorted_bkpts[0]->pos);
                $low_end_partner = $sorted_bkpts[1]->mate;
                $high_end_partner = $sorted_bkpts[2]->mate;
            }
            else {
                # This applies to all of:
                # - 0+/2-/2+
                # - 0+,2+/2-
                # - 0+,2-/2+
                # - 1-,2-/2+
                ($start, $end) = ($sorted_bkpts[0]->pos, $sorted_bkpts[1]->pos);
                $dist = max(1, $sorted_bkpts[2]->pos - $sorted_bkpts[1]->pos);
                $low_end_partner = $sorted_bkpts[0]->mate;
                $high_end_partner = $sorted_bkpts[1]->mate;
            }
        }

        @out_fields = ($params{tag});
        push @out_fields, $params{sample_name};
        push @out_fields, $rg_pattern->{rg_string};
        push @out_fields, $footprint->cluster->id;
        push @out_fields, $footprint->id;
        push @out_fields, $footprint->chr_name;
        push @out_fields, $start;
        push @out_fields, $end;
        push @out_fields, ($rg_pattern->{is_inverted} ? "-" : "+");
        push @out_fields, max($end - $start, 1);
        push @out_fields, $dist;
        push @out_fields, $footprint->cluster->size;
        push @out_fields, $footprint->size;
        push @out_fields, ($low_end_partner->chr_name . ":" . $low_end_partner->pos . ":" . $low_end_partner->dir);
        push @out_fields, ($high_end_partner->chr_name . ":" . $high_end_partner->pos . ":" . $high_end_partner->dir);
        push @out_fields, $low_end_partner->footprint->cluster->size;
        push @out_fields, $high_end_partner->footprint->cluster->size;
        push @out_fields, join(",", map { $_ or "NA" } @{$rg_pattern->{cn_changes}});
        push @out_fields, join(",", map { $_ or "NA" } @{$rg_pattern->{relative_cns}});

        print join("\t", @out_fields) . "\n";
    }
    
}

sub find_inversion_nested_in_tds {
    # Find the following cases:
    # 1. SV clusters with exactly three SVs in a single chromosome, if their
    #    simple rearrangement pattern is 1-,3+/1+,2+/2-,3-.
    # 2. Direct inversion SVs nested within a TD.

    my $rg_pattern;
    my %params = @_;

    my $MAX_TD_SIZE = 1e7;

    # Inversion-within-TD SV clusters
    for my $cluster (values %{$genome_of_cn_segs->{rg_clusters}}) {
        if ($cluster->n_rgs != 3) {
            next;
        }
        if (!$cluster->is_intra_chromosomal) {
            next;
        }
        my @rg_ends_pos = map {$_->pos} $cluster->rg_ends_array;
        if (max(@rg_ends_pos) - min(@rg_ends_pos) > $MAX_TD_SIZE) {
            next;
        }
        print STDERR $cluster->id . " " . $cluster->n_rgs . " " . join(":", map { $_->id } $cluster->rgs_array) . "\n";
        $rg_pattern = Footprint->simple_rg_pattern_from_rgs(
            [$cluster->rgs_array],
            %params
        );
        if (
                defined($rg_pattern) &&
                $rg_pattern->{rg_string} eq "1-,5+/1+,3+/3-,5-"
        ) {
            print join(
                "\t",
                $params{tag},
                $params{sample_name},
                $cluster->id,
                join(",", map {$_->id} $cluster->rgs_array),
            ) . "\n";
        }
    }

    # Direct inversion within a TD rearrangement
    my @td_rgs = ();
    my $td_rg;
    sub td_size {
        my $td = shift;
        return abs($td->high_end->pos - $td->low_end->pos);
    }
    for my $cluster (values %{$genome_of_cn_segs->{rg_clusters}}) {
        if (
                $cluster->cluster_classification(%params) eq "td" and
                td_size($cluster->rgs_array) <= $MAX_TD_SIZE
        ) {
            push @td_rgs, $cluster->rgs_array;
        }
    }
    for my $cluster (values %{$genome_of_cn_segs->{rg_clusters}}) {
        if ($cluster->cluster_classification(%params) ne "direct_inversion") {
            next;
        }
        for $td_rg (@td_rgs) {
            if (
                    $td_rg->low_end->chr_name eq (($cluster->rg_ends_array)[0])->chr_name  and
                    scalar(grep {$td_rg->low_end->pos > $_->pos} $cluster->rg_ends_array) == 0  and
                    scalar(grep {$td_rg->high_end->pos < $_->pos} $cluster->rg_ends_array) == 0
            ) {
                print join(
                    "\t",
                    $params{tag},
                    $params{sample_name},
                    $cluster->id . "," . $td_rg->cluster->id,
                    join(",", map {$_->id} $cluster->rgs_array) . "," . join(",", map {$_->id} $td_rg->cluster->rgs_array),
                ) . "\n";
            }
        }
    }
}

sub cn_shift_of_local_two_jump_chrs {
    # For each local two-jump type, find all chromosomes that have this
    # pattern, then calculate their average copy number shift relative
    # to genome-wide chromosomal copy number.
    my %params = @_;

    my $MAX_PATTERN_CHR_COUNT = 5;

    my $genome_cn = $genome_of_cn_segs->average_chr_cn;

    for my $pattern ("direct_inversion", "inversion_gain_loss",
            "inverted_duplication", "dup_trp_dup", "fb_then_fb", "inter") {
        my $diff_sum = 0;
        my $diff_count = 0;
        my %chrs_seen = ();
        for my $cluster (grep {$_->cluster_classification(%params) eq $pattern}
                values(%{$genome_of_cn_segs->{rg_clusters}})) {
            next if exists($chrs_seen{($cluster->footprints_array)[0]->chrom});
            $chrs_seen{($cluster->footprints_array)[0]->chrom} = 1;
            $diff_sum += ($cluster->footprints_array)[0]->chrom->mean_cn - $genome_cn;
            $diff_count += 1;
        }
        if ($diff_count > 0 and $diff_count <= $MAX_PATTERN_CHR_COUNT) {
            print join(
                "\t",
                $params{tag},
                $params{sample_name},
                $pattern,
                $diff_sum / $diff_count,
            ) . "\n";
        }
    }
}

sub cn_shift_of_size_two_footprint_chrs {
    # For each size two footprint involving two rearrangements, find all
    # chromosomes that have this pattern, then calculate their average copy
    # number shift relative to genome-wide chromosomal copy number. 
    my %params = @_;

    my $MAX_PATTERN_CHR_COUNT = 5;

    my $genome_cn = $genome_of_cn_segs->average_chr_cn;
    my($cluster, $footprint, $rg_pattern);
    my %chrs_seen;

    for my $rg_string ("0+/1+", "0+/2-", "1-/1+") {
        my $diff_sum = 0;
        my $diff_count = 0;
        %chrs_seen = ();
        for $cluster (values(%{$genome_of_cn_segs->{rg_clusters}})) {
            for $footprint ($cluster->footprints_array) {
                $rg_pattern = $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1);
                next if !defined($rg_pattern);
                next if $rg_pattern->{rg_string} ne $rg_string;
                next if exists($chrs_seen{$footprint->chrom});
                $chrs_seen{$footprint->chrom} = 1;
                $diff_sum += $footprint->chrom->mean_cn - $genome_cn;
                $diff_count += 1;
            }
        }
        if ($diff_count > 0 and $diff_count <= $MAX_PATTERN_CHR_COUNT) {
            print join(
                "\t",
                $params{tag},
                $params{sample_name},
                $rg_string,
                $diff_sum / $diff_count,
            ) . "\n";
        }
    }
}

sub size_3_footprint_partners {
    # The footprint of the outreaching rearrangement and other information
    # related to size 3 footprints. 
    my $genome = shift;
    my %params = @_;
    my ($cluster, $footprint);

    sub print_info_of_size_3_footprint_partners {
        my $footprint = shift;
        my %params = @_;
        
        # Get the outreaching breakpoint
        my $outreaching_bkpt;
        my $rg_pattern = $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1);
        my $rg_string = $rg_pattern->{rg_string};
        if ($rg_string eq "A+^C+/C-") {
            $outreaching_bkpt = ($footprint->sorted_rg_ends_array(%params))[1];
        }
        elsif ($rg_string eq "B-^C-/C+") {
            if (not $rg_pattern->{is_inverted}) {
                $outreaching_bkpt = ($footprint->sorted_rg_ends_array(%params))[2];
            }
            else {
                $outreaching_bkpt = ($footprint->sorted_rg_ends_array(%params))[0];
            }
        }
        elsif ($rg_string eq "A+/C-/C+") {
            if (not $rg_pattern->{is_inverted}) {
                $outreaching_bkpt = ($footprint->sorted_rg_ends_array(%params))[0];
            }
            else {
                $outreaching_bkpt = ($footprint->sorted_rg_ends_array(%params))[2];
            }
        }
        elsif ($rg_string eq "B-/B+/C+") {
            if (not $rg_pattern->{is_inverted}) {
                $outreaching_bkpt = ($footprint->sorted_rg_ends_array(%params))[2];
            }
            else {
                $outreaching_bkpt = ($footprint->sorted_rg_ends_array(%params))[0];
            }
        }
        else {
            return;
        }

        my $partner_fp = $outreaching_bkpt->mate->footprint;
        my $partner_rg_pattern = $partner_fp->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1);

        my @sorted_rg_ends = $footprint->sorted_rg_ends_array(%params);
        if ($rg_pattern->{is_inverted}) {
            @sorted_rg_ends = reverse(@sorted_rg_ends);
        }

        print join(
            "\t",
            $params{tag},
            $params{sample_name},
            $footprint->id,
            $footprint->cluster->id,
            $footprint->cluster->size,
            $footprint->cluster->n_rgs,
            $footprint->chr_name,
            $footprint->first_pos,
            $footprint->last_pos,
            join(",", map {$_->rg->id} @sorted_rg_ends),
            $rg_string,
            $rg_pattern->{is_inverted},
            $partner_fp->id,
            $partner_fp->size,
            ($partner_fp->chr_name eq $footprint->chr_name ? "=" : $partner_fp->chr_name),
            $partner_fp->first_pos,
            $partner_fp->last_pos,
            $partner_fp->detailed_type(%params),
            (defined($partner_rg_pattern) ? $partner_rg_pattern->{rg_string} : "NA"),
            (defined($partner_rg_pattern) ? $partner_rg_pattern->{is_inverted} : "NA"),
        )  . "\n";
    }

    # Print header
    print join(
        "\t",
        $params{tag},
        "sample_name",
        "footprint_id",
        "cluster_id",
        "n_footprints_in_cluster",
        "n_svs_in_cluster",
        "footprint_chr",
        "footprint_start",
        "footprint_end",
        "footprint_rgs",
        "rg_pattern",
        "pattern_is_inverted",
        "partner_footprint_id",
        "partner_footprint_n_bkpts",
        "partner_footprint_chr",
        "partner_footprint_start",
        "partner_footprint_end",
        "partner_footprint_type",
        "partner_footprint_pattern",
        "partner_footprint_is_inverted",
    ), "\n";

    for $cluster (values(%{$genome_of_cn_segs->{rg_clusters}})) {
        for $footprint ($cluster->footprints_array) {
            next unless defined($footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1));

            print_info_of_size_3_footprint_partners($footprint, %params);
        }
    }
}

sub four_bkpt_footprint_partners {
    # The footprint of the outreaching rearrangement and other information
    # related to four footprints. 

    my $genome = shift;
    my %params = @_;
    my ($cluster, $footprint);

    sub print_info_of_four_bkpt_footprints {
        my $footprint = shift;
        my %params = @_;
        my @rg_ends = $footprint->sorted_rg_ends_array(%params);
        my $footprint_pattern = $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1);
        @rg_ends = reverse(@rg_ends) if $footprint_pattern->{is_inverted};
        my $partner_info;
        my($first_end, $second_end);

        # Only interested in the following footprint patterns
        if ($footprint_pattern->{rg_string} eq "A+^C+/C-/E-") {
            ($first_end, $second_end) = @rg_ends[1, 3];
        }
        elsif ($footprint_pattern->{rg_string} eq "A+^D+/C-/D-") {
            ($first_end, $second_end) = @rg_ends[1, 2];
        }
        elsif ($footprint_pattern->{rg_string} eq "B-^D-/B+/D+") {
            ($first_end, $second_end) = @rg_ends[1, 3];
        }
        elsif ($footprint_pattern->{rg_string} eq "B-^C-/C+/D+") {
            ($first_end, $second_end) = @rg_ends[2, 3];
        }
        elsif ($footprint_pattern->{rg_string} eq "A+^E-/C-/C+") {
            ($first_end, $second_end) = @rg_ends[1, 2];
        }
        elsif ($footprint_pattern->{rg_string} eq "A+/C-/C+/E-") {
            ($first_end, $second_end) = @rg_ends[1, 2];
        }
        else {
            return;
        }

        $partner_info = join(
            "\t",
            $params{sample_name} . ":" . $first_end->mate->footprint->id,
            $first_end->mate->footprint->detailed_type(%params),
            $params{sample_name} . ":" . $second_end->mate->footprint->id,
            $second_end->mate->footprint->detailed_type(%params),
        );

        print join(
            "\t",
            $params{tag},
            $params{sample_name},
            $footprint->id,
            $footprint->cluster->id,
            $footprint->cluster->size,
            $footprint->cluster->n_rgs,
            $footprint->chr_name,
            $footprint->first_pos,
            $footprint->last_pos,
            $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1)->{rg_string},
            $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1)->{is_inverted},
            $partner_info,
        )  . "\n";
    }

    # Print header
    print join(
        "\t",
        $params{tag},
        "sample_name",
        "footprint_id",
        "cluster_id",
        "n_footprints_in_cluster",
        "n_svs_in_cluster",
        "footprint_chr",
        "footprint_start",
        "footprint_end",
        "rg_pattern",
        "pattern_is_inverted",
        "partner_1_footprint_id",
        "partner_1_footprint_pattern",
        "partner_2_footprint_id",
        "partner_2_footprint_pattern",
    ) . "\n";

    for $cluster (values(%{$genome_of_cn_segs->{rg_clusters}})) {
        for $footprint ($cluster->footprints_array) {
            next unless defined($footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1));
            print_info_of_four_bkpt_footprints($footprint, %params);
        }
    }
}


sub shard_transloc_footprint_partners {
    # Partners of A+/C-/C+ footprints

    my $genome = shift;
    my %params = @_;
    my ($cluster, $footprint);

    sub print_info_of_shard_transloc_footprints {
        my $footprint = shift;
        my %params = @_;
        my @rg_ends = $footprint->sorted_rg_ends_array(%params);
        my $footprint_pattern = $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1);
        @rg_ends = reverse(@rg_ends) if $footprint_pattern->{is_inverted};
        my $partner_info;

        $partner_info = join(
            "\t",
            $params{sample_name} . ":" . $rg_ends[0]->mate->footprint->id,
            $rg_ends[0]->mate->footprint->detailed_type(%params),
            $params{sample_name} . ":" . $rg_ends[1]->mate->footprint->id,
            $rg_ends[1]->mate->footprint->detailed_type(%params),
            $params{sample_name} . ":" . $rg_ends[2]->mate->footprint->id,
            $rg_ends[2]->mate->footprint->detailed_type(%params),
        );

        print join(
            "\t",
            $params{tag},
            $params{sample_name},
            $footprint->id,
            $footprint->cluster->id,
            $footprint->cluster->size,
            $footprint->cluster->n_rgs,
            $footprint->chr_name,
            $footprint->first_pos,
            $footprint->last_pos,
            $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1)->{rg_string},
            $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1)->{is_inverted},
            $partner_info,
        )  . "\n";
    }

    # Print header
    print join(
        "\t",
        $params{tag},
        "sample_name",
        "footprint_id",
        "cluster_id",
        "n_footprints_in_cluster",
        "n_svs_in_cluster",
        "footprint_chr",
        "footprint_start",
        "footprint_end",
        "rg_pattern",
        "pattern_is_inverted",
        "partner_1_footprint_id",
        "partner_1_footprint_pattern",
        "partner_2_footprint_id",
        "partner_2_footprint_pattern",
        "partner_3_footprint_id",
        "partner_3_footprint_pattern",
    ) . "\n";

    for $cluster (values(%{$genome_of_cn_segs->{rg_clusters}})) {
        for $footprint ($cluster->footprints_array) {
            next unless defined($footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1));
            my $rg_pattern = $footprint->rg_pattern(%params, no_bal_bkpts => 1, symbol_seg_names => 1);
            next unless $rg_pattern->{rg_string} eq "A+/C-/C+";
            print_info_of_shard_transloc_footprints($footprint, %params);
        }
    }
} # sub shard_transloc_footprint_partners


sub print_bal_olap_bkpt_footprints {
    my $genome = shift;
    my %params = @_;

    for my $cluster (values(%{$genome->{rg_clusters}})) {
        for my $footprint ($cluster->footprints_array) {
            if ($footprint->type(%params) eq "balanced_olap") {
                print(join(
                    "\t",
                    $params{tag},
                    $params{sample_name},
                    $footprint->chr_name,
                    $footprint->start,
                    $footprint->end
                ), "\n");
            }
        }
    }
}


sub print_chromoplexy_events {
    # Print out all chromoplexy events
    my $genome = shift;
    my %params = @_;

    # Print header
    print(join(
        "\t",
        $params{tag},
        "sample",
        "cluster_id",
        "chromoplexy_event_id",
        "balanced_breakpoint_count",
        "templated_insertion_count",
        "chain_or_cycle",
        "footprint_index",
        "footprint_type",
        "footprint_chrom",
        "footprint_start",
        "footprint_end",
        "sv_1",
        "sv_2",
    ) , "\n");

    # Algorithm: go through all footprints to find balanced breakpoints. Then
    # walk through them to find chromoplexy chains and cycles.
    my %seen_footprints = ();
    my $chromoplexy_event_id = 0;
    for my $cluster (values($genome->{rg_clusters})) {
    for my $footprint ($cluster->footprints_array) {
        if (not $footprint->is_balanced_type(%params)) {
            next;
        }
        if (exists($seen_footprints{$footprint->id})) {
            next;
        }

        print STDERR "FOUND A BALANCED BREAKPOINT: " . $footprint->to_s(%params) . "\n";

        # Now we have a balanced breakpoint.
        my $starting_footprint = $footprint;
        $seen_footprints{$starting_footprint->id} = 1;
        my @all_footprints = ($starting_footprint);

        # Record whether we should jump from low or high end next
        my $next_jump_from = 'low';

        # Go backwards, then forwards
        my @bkpts = ();
        my $new_footprint;
        my $is_cycle = 0;
        while (
                $all_footprints[0]->is_balanced_type(%params) or
                $all_footprints[0]->is_shard_type(%params)
        ) {
            @bkpts = $all_footprints[0]->sorted_rg_ends_array(%params);
            if ($next_jump_from eq 'low') {
                $new_footprint = $bkpts[0]->mate->footprint;
                if (($new_footprint->sorted_rg_ends_array(%params))[0] eq $bkpts[0]->mate) {
                    $next_jump_from = "high";
                }
                else {
                    $next_jump_from = "low";
                }
            }
            else {
                $new_footprint = $bkpts[1]->mate->footprint;
                if (($new_footprint->sorted_rg_ends_array(%params))[0] eq $bkpts[1]->mate) {
                    $next_jump_from = "high";
                }
                else {
                    $next_jump_from = "low";
                }
            }

            $seen_footprints{$new_footprint->id} = 1;
            if ($new_footprint eq $starting_footprint) {
                # We found a cycle
                $is_cycle = 1;
                last;
            }
            else {
                unshift @all_footprints, $new_footprint;
            }
        }

        # Only go forward if the current thing is not a cycle
        $next_jump_from = 'high';
        if (!$is_cycle) {
            while (
                    $all_footprints[-1]->is_balanced_type(%params) or
                    $all_footprints[-1]->is_shard_type(%params)
            ) {
                @bkpts = $all_footprints[-1]->sorted_rg_ends_array(%params);
                if ($next_jump_from eq 'low') {
                    $new_footprint = $bkpts[0]->mate->footprint;
                    if (($new_footprint->sorted_rg_ends_array(%params))[0] eq $bkpts[0]->mate) {
                        $next_jump_from = "high";
                    }
                    else {
                        $next_jump_from = "low";
                    }
                }
                else {
                    $new_footprint = $bkpts[1]->mate->footprint;
                    if (($new_footprint->sorted_rg_ends_array(%params))[0] eq $bkpts[1]->mate) {
                        $next_jump_from = "high";
                    }
                    else {
                        $next_jump_from = "low";
                    }
                }
                push @all_footprints, $new_footprint;
                $seen_footprints{$new_footprint->id} = 1;
            }
        }

        # Add first and last footprint to seen footprints
        $seen_footprints{$all_footprints[-1]->id} = 1;
        $seen_footprints{$all_footprints[$#all_footprints]->id} = 1;

        # Print the current chromoplexy cycle
        my $footprint_idx = -1;
        while (++$footprint_idx <= $#all_footprints) {
            my $footprint_type;
            if ($all_footprints[$footprint_idx]->is_balanced_type(%params)) {
                $footprint_type = "balanced";
            }
            elsif ($all_footprints[$footprint_idx]->is_shard_type(%params)) {
                $footprint_type = "shard";
            }
            elsif ($all_footprints[$footprint_idx]->type(%params) eq "single") {
                $footprint_type = "single";
            }
            else {
                $footprint_type = "other";
            }
            print(join(
                "\t",
                $params{tag},
                $params{sample_name},
                $starting_footprint->cluster->id,
                $chromoplexy_event_id,
                scalar(grep {$_->is_balanced_type(%params)} @all_footprints),
                scalar(grep {$_->is_shard_type(%params)} @all_footprints),
                ($is_cycle ? "cycle" : "chain"),
                $footprint_idx,
                $footprint_type,
                $all_footprints[$footprint_idx]->chr_name,
                $all_footprints[$footprint_idx]->first_pos,
                $all_footprints[$footprint_idx]->last_pos,
                (
                    ($footprint_type =~ /^(balanced|shard)$/) ?
                    ($all_footprints[$footprint_idx]->sorted_rg_ends_array(%params))[0]->rg->id :
                    "NA"
                ),
                (
                    ($footprint_type =~ /^(balanced|shard)$/) ?
                    ($all_footprints[$footprint_idx]->sorted_rg_ends_array(%params))[1]->rg->id :
                    "NA"
                ),
            ), "\n");
        }
        $chromoplexy_event_id++;
    }
    }
}


sub print_local_n_jump_footprints {
    # Here local n-jumps are defined as those with at least one inversion-type
    # SV, where either all n SVs in a single footprint, or where all n SVs
    # are within 5Mb.

    my $genome = shift;
    my $target_n_rgs = shift;
    my %params = @_;

    print_footprint_details_header();
    for my $cluster (values($genome->{rg_clusters})) {
        if (
                $cluster->size == 1  and
                $cluster->n_rgs == $target_n_rgs
        ) {
            print_footprint_details($cluster->footprints_array, %params);
        }
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
