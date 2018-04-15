#!/usr/bin/perl
#
# This script removes shards that are bypassed by other breakpoints. 
#
# Usage: perl bypass_shards_from_bedpe.pl rgs.bedpe rg_cns.txt
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

# Set default options
my %opts = (
    # Basic options
    verbose => 0,
    help => 0,

    # Advanced options for extracting rearrangement patterns
    max_balanced_rg_dist => 500,  # Distance for breakpoint to be considered balanced
    max_balanced_rg_overlap => 100,  # Distance for breakpoint to be considered overlap-balanced
    max_shard_length => 2000,  # Only need to worry about shards that can be bypassed
    shard_bypassing_slop => 200,  # Distance between 
    min_cn_bkpt_seg_size => 50000,
    min_cn_change => 0.3,
);

my $rgs_file = shift;
my $cn_file  = shift;


sub stderr_with_date {
    print STDERR "[" . `echo -n \`date\`` . "] " . $_[0];
}

# First read in all rearrangement and copy number data into a genome. 
stderr_with_date("Reading data in...\n") if $opts{verbose};
my $genome_of_cn_segs = CopyNumberSegmentGenome->new_from_BEDPE_and_CN_BED($rgs_file, $cn_file, %opts);

stderr_with_date("Sorting segments on each chromosome...\n") if $opts{verbose};
$genome_of_cn_segs->sort_segments;

stderr_with_date("Bypassing shards...\n") if $opts{verbose};
my %removed_rgs = $genome_of_cn_segs->preprocess_rearrangement_and_copy_number_data(%opts, "no_shard_bypassing" => 0, "no_rg_filtering" => 1);

print STDERR <<EOF;
Removed rearrangements:
Number of iterations: $removed_rgs{iter}
Copy number segmentation breakpoints with CN change < $opts{min_cn_change} : $removed_rgs{cn_removed}
Copy number segments smaller than $opts{min_cn_bkpt_seg_size}: $removed_rgs{small_cn_bkpt_seg}
shard_bypassing_svs => $removed_rgs{shard_bypassing_svs}
no_cn_change_svs => $removed_rgs{no_cn_change_svs}

EOF

stderr_with_date("Outputting rearrangements...\n") if $opts{verbose};
my @sorted_rgs = sort { $a->low_end->chr_name cmp $b->low_end->chr_name || $a->high_end->pos <=> $b->high_end->pos} values(%{$genome_of_cn_segs->{rg_of_id}});
for my $rg (@sorted_rgs) {
    print join("\t", @{$rg->{orig_data}}) . "\n";
}
