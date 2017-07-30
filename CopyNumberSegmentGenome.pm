use strict;
use warnings FATAL => 'all';

use lib '/nfs/users/nfs_y/yl3/programs/scripts/pancan_related_scripts/';
use CopyNumberSegmentArray;
use CopyNumberSegment;
use Rearrangement;
use RearrangementConnectedComponent;
use RearrangementGroup;
use Footprint;
use FootprintCluster;
use List::Util qw(sum min max);

package CopyNumberSegmentGenome;

sub new {
    my $class = shift;
    my %params = @_;
    
    return bless {
        chrs => [],
        rg_of_id => ($params{rg_of_id} or undef)
    }, $class;
}

sub new_from_BEDPE_and_CN_BED {
    my $class = shift;
    my $rgs_file = shift;
    my $cn_file = shift;
    my %params = @_;
    my %rg_of_id;
    my @F;
    my(@reads_min_max_pos, @reads_clipped);

    open RGS, $rgs_file;
    while (<RGS>) {
        chomp;
        @F = split /\t/;
        @reads_min_max_pos = @reads_clipped = ();
        for (@F[12..15]) {
        # for (@F[14..17]) {  ## JUST A QUICK FIX!!!
            /^(\d+) \((\d+)\)$/ or die;
            push @reads_min_max_pos, $1;
            push @reads_clipped, $2;
        }
        $rg_of_id{$F[6]} = Rearrangement->new(
            low_end_dir => $F[8],
            high_end_dir => $F[9],
            id => $F[6],
            orig_data => [@F],
            reads_pos => \@reads_min_max_pos,
            reads_clip => \@reads_clipped,
        );
    }
    close RGS;

    my ($chr, $start_pos, $end_pos, $cn, $low_end_bkpt, $high_end_bkpt, $n_win);
    my (%segments_of_chr, $low_end, $high_end);
    if ($cn_file =~ /\.gz$/) {
        open CN, "gunzip -c $cn_file |";
    }
    else {
        open CN, $cn_file;
    }
    while (<CN>) {
        chomp;
        ($chr, $start_pos, $end_pos, $cn, $low_end_bkpt, $high_end_bkpt, $n_win) = split /\t/;
        # $start_pos++;  # Convert 0-based BEDPE format start back to 1-based. # DOESN'T NEED TO BE DONE

        if ($low_end_bkpt =~ /:/ && !exists($rg_of_id{substr($low_end_bkpt, 0, -2)})) {
            # warn "Rearrangement $low_end_bkpt in CN file not found in RG file. Changing into cn_bkpt";
            print STDERR "$low_end_bkpt ";
            $low_end_bkpt = "cn_bkpt";
        }
        if ($high_end_bkpt =~ /:/ && !exists($rg_of_id{substr($high_end_bkpt, 0, -2)})) {
            # warn "Rearrangement $high_end_bkpt in CN file not found in RG file. Changing into cn_bkpt";
            print STDERR "$high_end_bkpt ";
            $high_end_bkpt = "cn_bkpt";
        }

        if (!(exists $segments_of_chr{$chr})) {
            $segments_of_chr{$chr} = CopyNumberSegmentArray->new(name => $chr);
        }

        if ($n_win eq "NA") {
            print STDERR "Warning: encountered segment $chr:$start_pos-$end_pos (size: " . ($end_pos - $start_pos + 1) . ") with n_win 'NA'. n_win set to 1.\n" if $params{verbose} >= 2;
            $n_win = 1;
        }
        $segments_of_chr{$chr}->add_copy_number_segment(
            CopyNumberSegment->new(
                chr           => $segments_of_chr{$chr},
                start         => $start_pos,
                end           => $end_pos,
                cn            => ($cn eq "NA" ? undef : $cn),
                n_win         => $n_win,
                low_end_bkpt  => $low_end_bkpt,
                high_end_bkpt => $high_end_bkpt,
                rg_of_id      => \%rg_of_id,
            )
        );
    }
    close CN;

    # Sanity check: after reading rearrangement data in, each rearrangement end
    # in %rg_of_id should be associated with a copy number segment end.
    for (keys %rg_of_id) {
        if (
            !defined($rg_of_id{$_}->low_end->{segment_end}) &&
            !defined($rg_of_id{$_}->high_end->{segment_end})
        ) {
            print STDERR "Rearrangement '$_' was not found in RG_CNS file and will be removed from \%rg_of_id.\n" if $params{verbose} >= 2;
            delete $rg_of_id{$_};
        }
        elsif (!defined($rg_of_id{$_}->low_end->{segment_end})) {
            warn "Rearrangement $_ low end was not assigned to any copy number segment!";
        }
        elsif (!defined($rg_of_id{$_}->high_end->{segment_end})) {
            warn "Rearrangement $_ high end was not assigned to any copy number segment!";
        }
    }

    return bless {
        chrs => \%segments_of_chr,
        rg_of_id => \%rg_of_id,
        connected_components => undef,
    }, $class;
}

sub read_events_from_clustering_file {
    # Annotate rearrangements with cluster and footprint information.
    # $self->{rg_clusters} is added and computed. 
    #
    # This function reads in a clustering result file, which current is
    # generated using clustering_index.R. This input file is in BEDPE format
    # with the following augmented columns.
    # 11: Rearrangement cluster ID
    # 12: Number of rearrangements in the rearrangement cluster (not used in
    #     this function)
    # 13: Footprint ID of the low end breakpoint
    # 14: Footprint ID of the high end breakpoint
    # 15: Start and end coordinates of the low end footprint
    # 16: Start and end coordinates of the high end footprint

    my $self = shift;
    my $clusters_file = shift;
    open IN, $clusters_file or die $!;
    my @F;
    my($rg_id, $cluster_id, $footprint_l_id, $footprint_h_id, $footprint_l_coord, $footprint_h_coord);
    my($footprint_start, $footprint_end);
    my($cluster_count, $footprint_count) = (0, 0);
    $self->{rg_clusters} = {};
    while (<IN>) {
        chomp;
        @F = split /\t/;
        ($rg_id, $cluster_id, $footprint_l_id, $footprint_h_id, $footprint_l_coord, $footprint_h_coord) = @F[6, 10, 12..15];
        # ($cluster_id, $footprint_l_id, $footprint_h_id, $footprint_l_coord, $footprint_h_coord) = @F[12, 14..17];  # DELETE LATER

        # First the low end. If the footprint cluster or footprint does not exist yet, then need to be created. 
        if (!exists($self->{rg_clusters}->{$cluster_id})) {
            $self->add_rg_cluster(cluster_id => $cluster_id, genome => $self);
            $cluster_count++;
        }

        # Low end footprint
        if (!exists($self->{rg_clusters}->{$cluster_id}->footprints->{$footprint_l_id})) {
            ($footprint_start, $footprint_end) = split "-", $footprint_l_coord;
            $self->{rg_clusters}->{$cluster_id}->add_footprint(
                Footprint->new(
                    chrom => $self->{chrs}->{$F[0]},
                    id => $footprint_l_id,
                    start => $footprint_start,
                    end => $footprint_end,
                    cluster => $self->{rg_clusters}->{$cluster_id},
                )
            );
            $footprint_count++;
        }
        $self->{rg_clusters}->{$cluster_id}->footprints->{$footprint_l_id}->add_rg_ends($self->rg_of_id->{$rg_id}->low_end);
        # $self->rg_of_id->{$rg_id}->low_end->{footprint} = $self->{rg_clusters}->{$cluster_id}->footprints->{$footprint_l_id};

        # Then same thing for high end. 
        if (!exists($self->{rg_clusters}->{$cluster_id}->footprints->{$footprint_h_id})) {
            ($footprint_start, $footprint_end) = split "-", $footprint_h_coord;
            $self->{rg_clusters}->{$cluster_id}->add_footprint(
                Footprint->new(
                    chrom => $self->{chrs}->{$F[3]},
                    id => $footprint_h_id,
                    start => $footprint_start,
                    end => $footprint_end,
                    cluster => $self->{rg_clusters}->{$cluster_id},
                )
            );
            $footprint_count++;
        }
        $self->{rg_clusters}->{$cluster_id}->footprints->{$footprint_h_id}->add_rg_ends($self->rg_of_id->{$rg_id}->high_end);
        # $self->rg_of_id->{$rg_id}->high_end->{footprint} = $self->{rg_clusters}->{$cluster_id}->footprints->{$footprint_h_id};
    }
    close IN;

    print STDERR sprintf("Read in %s clusters and %s footprints.\n", $cluster_count, $footprint_count);
}


#
# Getter methods
#
sub rg_of_id {
    my $self = shift;
    if (!defined($self->{rg_of_id})) {
        die "Attempted to call rg_of_id() to get undefined $self\->{rg_of_id}!";
    }
    return $self->{rg_of_id};
}

sub chrs {
    my $self = shift;
    # if (!defined($self->{chrs})) {
    #     die "Attempted to call $self\->chrs() to get undefined $self\->{chrs}!";
    # }
    return $self->{chrs};
}

sub chrs_array {
    my $self = shift;
    return(values %{$self->{chrs}});
}

sub average_chr_cn {
    my $self = shift;
    my $count = 0;
    my $sum = 0;
    for my $chr ($self->chrs_array) {
        $sum += $chr->mean_cn;
        $count++;
    }
    return $sum / $count;
}

sub connected_components {
    my $self = shift;
    if (!defined($self->{connected_components})) {
        die "Attempted to call connected_components() to get undefined $self\->{connected_components}!";
    }
    return $self->{connected_components};
}

sub connected_components_array {
    my $self = shift;
    return @{$self->connected_components};
}

sub chr_by_name {
    my $self = shift;
    my $chr_name = shift;
    if (!exists($self->chrs->{$chr_name})) {
        warn "No match for chromosome name '$chr_name' in $self\->chr_by_name()";
        return undef;
    }
    else {
        return $self->chrs->{$chr_name};
    }
}
#
# End of getter methods
#


#
# Worker subroutines
#
sub sort_segments {
    my $self = shift;
    for my $chr ($self->chrs_array) {
        $chr->sort_segments;
    }
}

sub delete_rg {
    # Deletes $rg_id from $self->{rg_of_id}, but before that makes all the
    # necessary updates in the affected copy number segment data etc. 
    my $self = shift;
    my $rg_id = shift;
    my $cur_rg = $self->rg_of_id->{$rg_id};
    $cur_rg->low_end->segment_end->{bkpt} = undef;
    $cur_rg->low_end->segment_end->{boundary} = "cn_bkpt";
    $cur_rg->low_end->segment_end_neighbour->{boundary} = "cn_bkpt";
    $cur_rg->high_end->segment_end->{bkpt} = undef;
    $cur_rg->high_end->segment_end->{boundary} = "cn_bkpt";
    $cur_rg->high_end->segment_end_neighbour->{boundary} = "cn_bkpt";
    delete($self->rg_of_id->{$rg_id});
}

sub add_rg_cluster {
    # Add a new rg cluster 
    my $self = shift;
    my %params = @_;
    my $cluster_id = $params{cluster_id};

    if (exists($self->{rg_clusters}->{$cluster_id})) {
        die "Attempted to add rearrangement cluster '$cluster_id' but it is also part of the current CopyNumberSegmentGenome";
    }

    $self->{rg_clusters}->{$cluster_id} = FootprintCluster->new(id => $cluster_id, genome => $self);
}

sub filter_rgs_by_no_cn_change {
    my $self = shift;
    my %params = @_;
    if (!exists($params{min_cn_change})) {
        die "Attempted to call $self\->filter_rgs_by_no_cn_change() without a parameter 'min_cn_change'";
    }
    if (!exists($params{keep_small_dels_and_tds})) {
        die "Attempted to call $self\->filter_rgs_by_no_cn_change() without a parameter 'keep_small_dels_and_tds'";
    }
    if (!exists($params{max_shard_length})) {
        die "Attempted to call $self\->filter_rgs_by_no_cn_change() without a parameter 'max_shard_length'";
    }

    my($cur_rg, $has_cn_change, $both_ends_balanced, @rgs_to_delete);
    for my $rg_id (keys %{$self->rg_of_id}) {
        $cur_rg = $self->rg_of_id->{$rg_id};

        if ($params{keep_small_dels_and_tds} && ($cur_rg->is_small_td(%params) || $cur_rg->is_small_del(%params))) {
            next;
        }

        # If the segments are too short for copy number to be estimated properly
        if (
            !$cur_rg->is_foldback(%params) &&
            (
                (   
                    $cur_rg->low_end->segment->length < 1000 ||
                    $cur_rg->low_end->segment_end_neighbour->segment->length < 1000
                ) &&
                (   
                    $cur_rg->high_end->segment->length < 1000 ||
                    $cur_rg->high_end->segment_end_neighbour->segment->length < 1000
                )
            )
        ) {
            $has_cn_change = 1;
        }
        elsif (defined($cur_rg->weighted_avg_cn_change_across_rg(%params))) {
            $has_cn_change = $cur_rg->weighted_avg_cn_change_across_rg(%params) >= $params{min_cn_change};
        }
        else {
            ## The fact that we got here means that the current rearrangement
            ## is a fold-back but had the exact same breakpoint with some other
            ## rearrangments. 
            $has_cn_change = 1;
        }

        if (!$cur_rg->is_foldback(%params)) {
            $both_ends_balanced = defined($cur_rg->low_end->balanced_bkpt_partner_rg_end(%params)) &&
                                  defined($cur_rg->high_end->balanced_bkpt_partner_rg_end(%params));

            if (!$has_cn_change && !$both_ends_balanced) {
                push @rgs_to_delete, $rg_id;
            }
        }
        elsif (
            !$has_cn_change &&
            !($cur_rg->low_end->dir eq "+" && defined($cur_rg->high_end->balanced_bkpt_partner_rg_end(%params))) &&
            !($cur_rg->low_end->dir eq "-" && defined($cur_rg->low_end->balanced_bkpt_partner_rg_end(%params)))
        ) {
            push @rgs_to_delete, $rg_id;
        }
    }

    for my $rg_id (@rgs_to_delete) {
        $self->delete_rg($rg_id);
    }

    return scalar(@rgs_to_delete);
}

sub remove_shard_sequence_bypassing_rgs {
    my $self = shift;
    my %params = @_;
    my @rgs_to_delete = ();
    my $rg_id;
    for $rg_id (keys %{$self->rg_of_id}) {
        if ($self->rg_of_id->{$rg_id}->is_shard_bypassing(%params)) {
            push @rgs_to_delete, $rg_id;
        }
    }

    for $rg_id (@rgs_to_delete) {
        select STDERR;
        print "Bypassing shard removed: ";
        $self->{rg_of_id}->{$rg_id}->print;
        select STDOUT;
        $self->delete_rg($rg_id);
    }

    return scalar(@rgs_to_delete);
}

sub remove_small_cn_bkpt_segments {
    my $self = shift;
    my %params = @_;
    for my $chr ($self->chrs_array) {
        $chr->remove_small_cn_bkpt_segments(%params);
    }
}

sub remove_cn_bkpts_without_cn_change {
    my $self = shift;
    my %params = @_;
    for my $chr ($self->chrs_array) {
        $chr->remove_cn_bkpts_without_cn_change(%params);
    }
}

sub normalise_single_bp_segments {
    my $self = shift;
    my %params = @_;
    for my $chr ($self->chrs_array) {
        $chr->normalise_single_bp_segments(%params);
    }
}

sub preprocess_rearrangement_and_copy_number_data {
    my $self = shift;
    my %params = @_;
    $self->sort_segments;

    my $changed = 1;
    my $iterations = 0;
    my($shard_bypassing_svs_removed, $cn_removed,
        $small_cn_bkpt_seg_removed, $no_cn_change_sv_removed) = (0, 0, 0, 0);
    while ($changed) {
        $changed = 0;
        $iterations++;
        if ($params{verbose} >= 2) {
            print STDERR "  Iterations: $iterations...\n";
        }
        $_ = $self->remove_cn_bkpts_without_cn_change(%params);
        $cn_removed += $_;
        $changed = 1 if $_;

        $_ = $self->remove_small_cn_bkpt_segments(%params);
        $small_cn_bkpt_seg_removed += $_;
        $changed = 1 if $_;

        if (!$params{no_shard_bypassing}) {
            $_ = $self->remove_shard_sequence_bypassing_rgs(%params);
            $shard_bypassing_svs_removed += $_;
            $changed = 1 if $_;
        }

        if (!$params{no_rg_filtering}) {
            $_ = $self->filter_rgs_by_no_cn_change(%params);
            $no_cn_change_sv_removed += $_;
            $changed = 1 if $_;
        }

        $_ = $self->normalise_single_bp_segments(%params);
        $changed = 1 if $_;
    }

    return (
        iter => $iterations,
        cn_removed => $cn_removed,
        small_cn_bkpt_seg => $small_cn_bkpt_seg_removed,
        shard_bypassing_svs => $shard_bypassing_svs_removed,
        no_cn_change_svs => $no_cn_change_sv_removed,
    );
}


#
# Creating new RG clusters from RGs - fixing the clustering
#

sub create_sv_cluster_from_rg {
    # Given a single rg, create an SV cluster and the footprints
    # out of it.
    my $self = shift;
    my $rg = shift;
    my $cluster_id = shift;
    my %params = @_;

    # Remove traces of old footprints from the rg
    my($bkpt_l, $bkpt_h) = sort(
        {$a->chr_name cmp $b->chr_name  or  $a->pos <=> $b->pos}
        ($rg->low_end, $rg->high_end)
    );
    if (defined($bkpt_l->{footprint})) {
        $bkpt_l->footprint->remove_rg_ends(%params, rg_ends => [$bkpt_l]);
    }
    if (defined($bkpt_h->{footprint})) {
        $bkpt_h->footprint->remove_rg_ends(%params, rg_ends => [$bkpt_h]);
    }

    # Create the new SV cluster, then new footprints
    $self->add_rg_cluster(%params, cluster_id => $cluster_id);

    my $footprint_l_id = "$cluster_id.0";
    $self->{rg_clusters}->{$cluster_id}->add_footprint(
        Footprint->new(
            chrom => $bkpt_l->chr,
            id => $footprint_l_id,
            start => $bkpt_l->pos,
            end => $bkpt_l->pos,
            cluster => $self->{rg_clusters}->{$cluster_id},
        )
    );
    $bkpt_l->{footprint} = undef;
    $self->{rg_clusters}->{$cluster_id}->footprints->{$footprint_l_id}->add_rg_ends($bkpt_l);
    $bkpt_l->{footprint} = $self->{rg_clusters}->{$cluster_id}->footprints->{$footprint_l_id};

    my $footprint_h_id = "$cluster_id.1";
    $self->{rg_clusters}->{$cluster_id}->add_footprint(
        Footprint->new(
            chrom => $bkpt_h->chr,
            id => $footprint_h_id,
            start => $bkpt_h->pos,
            end => $bkpt_h->pos,
            cluster => $self->{rg_clusters}->{$cluster_id},
        )
    );
    $bkpt_h->{footprint} = undef;
    $self->{rg_clusters}->{$cluster_id}->footprints->{$footprint_h_id}->add_rg_ends($bkpt_h);
    $bkpt_h->{footprint} = $self->{rg_clusters}->{$cluster_id}->footprints->{$footprint_h_id};
}

sub create_sv_cluster_from_footprints {
    # Given an array of footprints, create a new SV cluster.
    # Also removes the footprints from their existing cluster. 
    my $self = shift;
    my $cluster_id = shift;
    my @footprints = @{shift()};

    $self->add_rg_cluster(cluster_id => $cluster_id);
    my $new_cluster = $self->{rg_clusters}->{$cluster_id};
    for (@footprints) {
        if (exists($_->{cluster})) {
            $_->cluster->remove_footprint($_);
        }
        $_->{cluster} = $new_cluster;
    }
    $new_cluster->add_footprint(@footprints);
}

sub remove_flanking_dels_and_tds_from_footprint {
    # Peel off flanking TDs and deletions from a footprint. 
    my $self = shift;
    my $footprint = shift;
    my %params = @_;

    my $cluster = $footprint->cluster;
    my @bkpts;
    my($is_peelable_type, $is_flanking_event);
    my $peeled_id = 0;
    my($cluster_id, $footprint_l_id, $footprint_h_id);
    my $neighbour_bkpt;

    # If the footprint is literally a deletion or a TD, then just remove it.
    if (
            $footprint->size == 2 and
            $footprint->n_rgs == 1 and
            (
                ($footprint->rgs_array)[0]->is_del_type or
                ($footprint->rgs_array)[0]->is_td_type
            )
     ) {
        if ($cluster->size == 1) {
            return;
        }

        print STDERR sprintf("In SV cluster '%s', footprint '%s' is a simple deletion/TD. Removed.\n", $cluster->id, $footprint->id);
        print STDERR ($footprint->to_s(%params) . "\n");
        print STDERR "Created cluster " . $footprint->id . "\n";
        $self->create_sv_cluster_from_rg($footprint->rgs_array, $footprint->id, %params);

        return;
    }

    # Peel from the 5'-end
    my $keep_peeling = 1;
    while ($keep_peeling) {
        # Continue until first breakpoint is not a deletion or
        # TD that can be peeled. 
        $keep_peeling = 0;

        @bkpts = $footprint->sorted_rg_ends_array(%params);
        if (not @bkpts) {
            next;
        }
        $is_peelable_type = ($bkpts[0]->rg->is_del_type(%params) or $bkpts[0]->rg->is_td_type(%params));
        if (not $is_peelable_type) {
            next;
        }

        # Now a TD/deletion is to be peeled off. Need to create a new
        # cluster and footprints for the event, while updating the
        # remaining.
        if (@bkpts >= 3  and  $bkpts[0]->rg->id eq $bkpts[1]->rg->id) {
            # Both ends in the same footprint. Create a new cluster, and
            # add the two footprints. 
            $cluster_id = $footprint->id . "." . $peeled_id;
            # $footprint->remove_rg_ends(%params, rg_ends => [@bkpts[0, 1]]);
            $self->create_sv_cluster_from_rg(($bkpts[0])->rg, $cluster_id, %params);
            $peeled_id++;
            $keep_peeling = 1;

            print STDERR sprintf("From SV cluster '%s' footprint '%s', removed breakpoints (TD or del):\n", $cluster->id, $footprint->id);
            print STDERR ($bkpts[0]->to_s . "\n");
            print STDERR ($bkpts[1]->to_s . "\n");
            print STDERR "Created cluster $cluster_id\n";
        }
        elsif (
                @bkpts >= 2 and
                defined($bkpts[0]->closest_us_rg_end(%params, within => -1, exclude_rgs => [])) and
                $bkpts[0]->rg->id eq ($bkpts[0]->closest_us_rg_end(%params, within => -1, exclude_rgs => []))->rg->id and
                ($bkpts[0]->closest_us_rg_end(%params, within => -1, exclude_rgs => []))->footprint->size == 1
        ) {
            # One end in footprint, the other end just outside forming a solitary footprint
            $cluster_id = $footprint->id . "." . $peeled_id;
            # $footprint->remove_rg_ends(%params, rg_ends => [$bkpts[0]]);
            # $neighbour_bkpt = $bkpts[0]->closest_us_rg_end(%params, within => -1, exclude_rgs => []);
            # $neighbour_bkpt->footprint->remove_rg_ends(%params, rg_ends => [$neighbour_bkpt]);
            $self->create_sv_cluster_from_rg(($bkpts[0])->rg, $cluster_id, %params);
            $peeled_id++;
            $keep_peeling = 1;

            print STDERR sprintf("From SV cluster '%s' footprint '%s', removed breakpoint:\n", $cluster->id, $footprint->id);
            print STDERR ($bkpts[0]->to_s . "\n");
            print STDERR sprintf("Removed solitary partner of the above from footprint '%s':\n", $bkpts[0]->closest_us_rg_end(%params, within => -1, exclude_rgs => [])->footprint->id);
            print STDERR ($bkpts[0]->closest_us_rg_end(%params, within => -1, exclude_rgs => [])->to_s . "\n");
            print STDERR "Created cluster $cluster_id\n";
        }
    }

    # Peel from the 3'-end
    $keep_peeling = 1;
    while ($keep_peeling) {
        # Continue until first breakpoint is not a deletion or
        # TD that can be peeled. 
        $keep_peeling = 0;

        @bkpts = $footprint->sorted_rg_ends_array(%params);
        if (not @bkpts) {
            next;
        }
        $is_peelable_type = ($bkpts[-1]->rg->is_del_type(%params) or $bkpts[-1]->rg->is_td_type(%params));
        if (not $is_peelable_type) {
            next;
        }

        # Now a TD/deletion is to be peeled off. Need to create a new
        # cluster and footprints for the event, while updating the
        # remaining.
        if (@bkpts >= 3  and  $bkpts[-2]->rg->id eq $bkpts[-1]->rg->id) {
            # add the two footprints. 
            $cluster_id = $footprint->id . "." . $peeled_id;
            # $footprint->remove_rg_ends(%params, rg_ends => [@bkpts[-2, -1]]);
            $self->create_sv_cluster_from_rg(($bkpts[-1])->rg, $cluster_id, %params);
            $peeled_id++;
            $keep_peeling = 1;

            print STDERR sprintf("From SV cluster '%s' footprint '%s', removed breakpoints (TD or del):\n", $cluster->id, $footprint->id);
            print STDERR ($bkpts[-2]->to_s . "\n");
            print STDERR ($bkpts[-1]->to_s . "\n");
            print STDERR "Created cluster $cluster_id\n";
        }
        elsif (
                @bkpts >= 2 and
                defined($bkpts[-1]->closest_ds_rg_end(%params, within => -1, exclude_rgs => [])) and
                $bkpts[-1]->rg->id eq ($bkpts[-1]->closest_ds_rg_end(%params, within => -1, exclude_rgs => []))->rg->id and
                ($bkpts[-1]->closest_ds_rg_end(%params, within => -1, exclude_rgs => []))->footprint->size == 1
        ) {
            # One end in footprint, the other end just outside forming a solitary footprint
            $cluster_id = $footprint->id . "." . $peeled_id;
            # $footprint->remove_rg_ends(%params, rg_ends => [$bkpts[-1]]);
            # $neighbour_bkpt = $bkpts[-1]->closest_ds_rg_end(%params, within => -1, exclude_rgs => []);
            # $neighbour_bkpt->footprint->remove_rg_ends(%params, rg_ends => [$neighbour_bkpt]);
            $self->create_sv_cluster_from_rg(($bkpts[-1])->rg, $cluster_id, %params);
            $peeled_id++;
            $keep_peeling = 1;

            print STDERR sprintf("From SV cluster '%s' footprint '%s', removed breakpoint:\n", $cluster->id, $footprint->id);
            print STDERR ($bkpts[-1]->to_s . "\n");
            print STDERR sprintf("Removed solitary partner of the above from footprint '%s':\n", $bkpts[-1]->closest_ds_rg_end(%params, within => -1, exclude_rgs => [])->footprint->id);
            print STDERR ($bkpts[-1]->closest_ds_rg_end(%params, within => -1, exclude_rgs => [])->to_s . "\n");
            print STDERR "Created cluster $cluster_id\n";
        }
    }

    # If the remaining $footprint contains only a single TD or deletion type
    # rearrangement, then go ahead and separate it out.
    if (
            $footprint->n_rg_ends == 2 and
            $footprint->n_rgs == 1 and
            (
                ($footprint->rgs_array)[0]->is_del_type or
                ($footprint->rgs_array)[0]->is_td_type
            )
    ) {
        print STDERR sprintf(
            "From cluster %s, footprint %s is a single TD or del (%s). Move to its own cluster %s.\n",
            $footprint->cluster->id,
            $footprint->id,
            ($footprint->rgs_array)[0]->id,
            $footprint->id . "." . $peeled_id
        );

        $self->create_sv_cluster_from_rg(
            $footprint->rgs_array,
            $footprint->id . "." . $peeled_id,
            %params
        );
    }
}

sub remove_flanking_dels_and_tds {
    # Go through each footprint and for each of them, separate off flanking
    # TDs and deletions. 
    my $self = shift;
    my %params = @_;
    my @rg_clusters = values($self->{rg_clusters});
    my @footprints;
    for my $cluster (@rg_clusters) {
        @footprints = $cluster->footprints_array;
        for my $footprint (@footprints) {
            $self->remove_flanking_dels_and_tds_from_footprint($footprint, %params);
        }
    }
}

sub remove_all_del_clusters {
    # SV clusters that have only deletion-type events are separated into
    # individual deletion SV clusters. 
    my $self = shift;
    my %params = @_;
    my $rg;
    my $del_idx = 0;
    my $cluster_id;
    for my $cluster ($self->{rg_clusters}) {
        if (!$cluster->only_has_deletions) {
            next;
        }
        for $rg ($cluster->rgs_array) {
            $cluster_id = ($cluster->id . "." . $del_idx);
            $self->create_sv_cluster_from_rg($rg, $cluster_id, %params);
            print STDERR sprintf("From deletion-only SV cluster '%s', discarded following deletion:\n", $cluster->id);
            print STDERR $rg->to_s . "\n";
            print STDERR "Created cluster $cluster_id\n";
        }
    }
}

sub separate_out_shard_and_balanced_bkpt_cycles_from_cluster {
    # Given an SV cluster, find footprint/shard cycles and
    # separate them out of the cluster. 
    my $self = shift;
    my $cluster = shift;
    my %params = @_;

    # Some helper functions
    sub print_transferred {
        # Print out the list of footprints transferred from a cluster
        # into a new cluster.
        my $old_cluster_id = shift;
        my $new_cluster_id = shift;
        my $max_balanced_rg_overlap = shift;
        my @footprints = @_;
        print STDERR sprintf(
            "Following footprints were transferred from cluster %s to a new cluster %s:\n",
            $old_cluster_id,
            $new_cluster_id
        );
        for (@footprints) {
            print STDERR $_->to_s(max_balanced_rg_overlap => $max_balanced_rg_overlap) . "\n";
        }
    }  # End of print_transferred()

    sub get_balanced_or_shard_footprint_cycle {
        # Returns two values. 
        # 1. Whether the input footprint is part of a cycle made of
        #    templated insertions and balanced-type breakpoints.
        # 2. The list of traversed footprints. 
        my $initial_footprint = shift;
        my %params = @_;

        if (not ($initial_footprint->is_balanced_type(%params) or $initial_footprint->is_shard_type(%params))) {
            return(0, [$initial_footprint]);
        }

        my @footprint_list = ($initial_footprint);
        my $rg_end = (($initial_footprint->rg_ends_array)[0])->mate;
        my $footprint = $rg_end->footprint;

        # Loop. Either we loop back to the initial footprint, in which case
        # we have a cycle, or we stop looping. 
        while ($footprint->is_balanced_type(%params) or $footprint->is_shard_type(%params)) {
            if ($footprint == $initial_footprint) {
                # We went through the full cycle - create new cluster. 
                return(
                    1,
                    \@footprint_list,
                );
            }
            else {
                # Keep cycling
                push @footprint_list, $footprint;
                if ($footprint->is_balanced_type(%params)) {
                    $rg_end = $rg_end->balanced_bkpt_partner_rg_end_by_footprint(%params);
                }
                else {
                    $rg_end = $rg_end->shard_partner_rg_end_by_footprint(%params);
                }
                $rg_end = $rg_end->mate;
                $footprint = $rg_end->footprint;
            }
        }

        return(
            0,
            \@footprint_list,
        );
    }  # End of sub get_balanced_or_shard_footprint_cycle()


    my $cluster_id;
    my $new_cluster_idx = 0;
    my %cur_seen_footprints;
    my($is_cycle, $footprint_list_p);

    FOOTPRINT_CYCLE_SEARCH:
    while (1) {
        for my $initial_footprint ($cluster->footprints_array) {
            next if exists($cur_seen_footprints{$initial_footprint});
            ($is_cycle, $footprint_list_p) = get_balanced_or_shard_footprint_cycle($initial_footprint, %params);

            # Mark new footprints that we saw, so we can avoid having to traverse
            # them again. 
            for (@{$footprint_list_p}) {
                $cur_seen_footprints{$_} = 1;
            }

            if ($is_cycle and scalar(@{$footprint_list_p}) < $cluster->size) {
                $cluster_id = $cluster->id . "." . $new_cluster_idx;
                $new_cluster_idx++;
                print_transferred($cluster->id, $cluster_id, $params{max_balanced_rg_overlap}, @{$footprint_list_p});
                $self->create_sv_cluster_from_footprints($cluster_id, \@{$footprint_list_p});
                next FOOTPRINT_CYCLE_SEARCH;
            }
        }

        last;
    }
}  # End of sub separate_out_shard_and_balanced_bkpt_cycles_from_cluster()

sub separate_out_shard_clusters {
    # Separate out shard and balanced breakpoint cycles into
    # separate SV clusters.
    my $self = shift; 
    my %params = @_;
    for my $cluster (values($self->{rg_clusters})) {
        $self->separate_out_shard_and_balanced_bkpt_cycles_from_cluster($cluster, %params);
    }
}

sub normalise_single_rg_clusters {
    # If an SV cluster only has a single rearrangement, make sure it is in two
    # separate footprints. 

    my $self = shift;
    my %params = @_;

    # Ensure all single rg clusters have exactly two footprints
    for my $cluster (values(%{$self->{rg_clusters}})) {
        next if $cluster->n_rgs > 1;
        next if $cluster->size == 2;
        
        # Easiest way: recreate the cluster
        my $cluster_id = $cluster->id;
        my ($rg) = $cluster->rgs_array;
        $self->{rg_clusters}->{$cluster_id}->destroy;
        $self->create_sv_cluster_from_rg($rg, $cluster_id, %params);
    }
}

sub normalise_simple_dels_and_tds_clusters {
    # Separate TDs and deletions from SV clusters of type
    # 'simple_dels_and_tds' into individual clusters. 

    my $self = shift;
    my %params = @_;

    # Deal with clusters with only deletions and TDs
    my $rg;
    for my $cluster (values(%{$self->{rg_clusters}})) {
        if ($cluster->cluster_classification(%params) eq "simple_dels_and_tds") {
            print STDERR "Breaking down cluster " . $cluster->id . " of type "
                       . "simple_dels_and_tds into individual clusters for "
                       . "each rearrangement (n = " . $cluster->n_rgs . ").\n";
            my $idx = 0;
            while (exists($self->{rg_clusters}->{$cluster->id . ".$idx"})) {
                $idx++;
            }
            for $rg ($cluster->rgs_array) {
                $self->create_sv_cluster_from_rg(
                    $rg,
                    $cluster->id . ".$idx",
                    %params,
                );
                $idx++;
            }
            $cluster->destroy;
        }
    }
}

# sub find_inv_local_two_jump_clusters {
sub find_local_clusters {
    # OLD: If current cluster involves exactly 2 inversion-type SVs within 5Mb,
    # If the current cluster is contained within 5Mb,
    # then just bundle everything into a single footprint. 
    my $self = shift;
    my %params = @_;

    my $MAX_SIZE = 5e6;
    
    my %all_chrs = ();
    my $rg_end;
    for my $cluster (values(%{$self->{rg_clusters}})) {
        # # Two SVs exactly
        # At least two SVs
        if ($cluster->n_rgs == 1) {
            next;
        }

        # Require more than one footprint
        if (scalar($cluster->footprints_array) == 1) {
            next;
        }

        # All SVs on the same chromosome
        if (!$cluster->is_intra_chromosomal) {
            next;
        }

        # All breakpoints within 5Mb
        my @rg_ends_array = $cluster->rg_ends_array;
        if (List::Util::max(map {$_->pos} @rg_ends_array) - List::Util::min(map {$_->pos} @rg_ends_array) > $MAX_SIZE) {
            next;
        }

        # Now just merge the footprints
        my @footprints = $cluster->footprints_array;
        for (@footprints[1..$#footprints]) {
            for $rg_end ($_->rg_ends_array) {
                $rg_end->{footprint} = undef;
                $footprints[0]->add_rg_ends($rg_end);
            }
            $cluster->remove_footprint($_);
        }
        print STDERR sprintf(
            # "Cluster %s is a local inverted two-jump <= %s bp. Merging footprints %s into %s (%s).\n",
            "Cluster %s is a local event contained within %s bp. Merging footprints %s into %s (%s).\n",
            $cluster->id,
            $MAX_SIZE,
            join(", ", map {$_->id} @footprints[1..$#footprints]),
            $footprints[0]->id,
            $footprints[0]->to_s(%params)
        );
    }
}

sub rescue_shards_and_balanced_bkpt_footprints {
    # Sometimes shards and balanced breakpoints are separated into individual
    # footprints. Rescue them. 
    #
    # Exact algorithm:
    # If two adjacent footprints with exactly one breakpoint each (1) have their
    # breakpoints oriented in +- or -+ orientation, (2) are within $MAX_DISTANCE
    # and the next closest breakpoints of this cluster are further than
    # $DIST_MULTIPLIER times their distance away, then merge the footprints.
    my $self = shift;
    my %params = @_;

    my $MAX_DISTANCE = 5e6;
    my $DIST_MULTIPLIER = 3;  # Next closest breakpoint must be this times distance away
    
    my $cluster;
    my %footprints_of_chr;
    for $cluster (values %{$self->{rg_clusters}}) {
        if ($cluster->n_rgs == 1) {
            next;
        }

        # Collect all footprints by chromosome
        %footprints_of_chr = ();
        for my $footprint ($cluster->footprints_array) {
            push @{$footprints_of_chr{$footprint->chr_name}}, $footprint;
        }

        # Go through footprints of each chromosome and merge footprints
        my @footprints;
        my $cur_distance;
        my $rg_end;
        for my $chr_name (keys %footprints_of_chr) {
            @footprints = @{$footprints_of_chr{$chr_name}};
            @footprints = sort {$a->first_pos <=> $b->first_pos} @footprints;
            # print STDERR "Cluster " . $cluster->id . ", ", scalar(@footprints) . " footprints\n";
            # print STDERR $_->to_s(%params) . "\n" for @footprints;
            
            $_ = 0;
            while ($_ < $#footprints) {
                $cur_distance = $footprints[$_+1]->first_pos - $footprints[$_]->last_pos;
                my $both_footprints_are_singletons = (($footprints[$_]->size == 1)  and  ($footprints[$_+1]->size == 1));
                my $distance_within_cutoff = $cur_distance <= $MAX_DISTANCE;
                my $distance_to_upstream_ok
                        = ($_ == 0 or ($footprints[$_]->first_pos - $footprints[$_-1]->last_pos > $cur_distance * $DIST_MULTIPLIER));
                my $distance_to_downstream_ok
                        = ($_ == $#footprints-1 or ($footprints[$_+2]->first_pos - $footprints[$_+1]->last_pos > $cur_distance * $DIST_MULTIPLIER));
                my $opposing_orientations
                        = (($footprints[$_]->rg_ends_array)[0]->is_fwd xor ($footprints[$_+1]->rg_ends_array)[0]->is_fwd);
                my $different_rg_ids
                        = ($footprints[$_]->rg_ends_array)[0]->rg->id ne ($footprints[$_+1]->rg_ends_array)[0]->rg->id;
                if (
                        $both_footprints_are_singletons  and
                        $distance_within_cutoff  and
                        $distance_to_upstream_ok  and
                        $distance_to_downstream_ok  and
                        $opposing_orientations  and
                        $different_rg_ids
                ) {
                    # print STDERR ($footprints[$_]->size == 1) . "\n";
                    # print STDERR ($footprints[$_+1]->size == 1) . "\n";
                    # print STDERR $cur_distance . " " . ($cur_distance <= $MAX_DISTANCE) . "\n";
                    # print STDERR ($_ == 0 or $footprints[$_]->first_pos - $footprints[$_-1]->last_pos > $cur_distance * $DIST_MULTIPLIER) . "\n";
                    # print STDERR ($_ == $#footprints-1 or $footprints[$_+2]->first_pos - $footprints[$_+1]->last_pos > $cur_distance * $DIST_MULTIPLIER) . "\n";
                    print STDERR sprintf(
                        "Merging singleton footprints containing %s and %s together\n",
                        ($footprints[$_]->rg_ends_array)[0]->to_s,
                        ($footprints[$_+1]->rg_ends_array)[0]->to_s
                    );

                    $rg_end = ($footprints[$_+1]->rg_ends_array)[0];
                    $footprints[$_+1]->remove_rg_ends(%params, rg_ends => [$rg_end]);
                    $footprints[$_]->add_rg_ends($rg_end);
                    # $cluster->remove_footprint($footprints[$_+1]);  # This is done as part of Footprint->remove_rg_ends()
                    splice(@footprints, $_+1);
                }

                $_++;
            }
        }
    }
}

sub rescue_inversion_shards {
    # Rescue inversion-shard-type events. 
    # For each balanced breakpoint, if one of the SV is an inversion-type
    # event where the mate is in a singleton footprint right adjacent
    # to the balanced breakpoint within $MAX_DISTANCE, then merge the
    # singleton footprint into the balanced breakpoint footprint.
    my $self = shift;
    my %params = @_;

    my $MAX_DISTANCE = 5e6;

    # Collect all footprints by chromosome
    my $cluster;
    my %footprints_of_chr;
    for $cluster (values %{$self->{rg_clusters}}) {
        %footprints_of_chr = ();
        for my $footprint ($cluster->footprints_array) {
            push @{$footprints_of_chr{$footprint->chr_name}}, $footprint;
        }

        # Go through footprints of each chromosome and merge footprints
        my @footprints;
        my @rg_ends;
        for my $chr_name (keys %footprints_of_chr) {
            @footprints = @{$footprints_of_chr{$chr_name}};
            @footprints = sort {$a->first_pos <=> $b->first_pos} @footprints;
            
            $_ = 0;
            while ($_ <= $#footprints) {
                @rg_ends = $footprints[$_]->sorted_rg_ends_array(%params);
                if (not $footprints[$_]->is_balanced_type(%params)) {
                    $_++;
                    next;
                }

                if (
                        # Merging with footprint upstream?
                        $_ > 0  and
                        $footprints[$_-1]->size == 1  and
                        $footprints[$_]->first_pos - $footprints[$_-1]->last_pos < $MAX_DISTANCE  and
                        ($footprints[$_-1]->rg_ends_array)[0]->rg->is_mm_type  and
                        ($footprints[$_-1]->rg_ends_array)[0]->rg->id eq $rg_ends[1]->rg->id

                ) {
                    print STDERR sprintf(
                        "Merging solo footprint %s with footprint %s to rescue inversion shard\n",
                        ($footprints[$_-1]->rg_ends_array)[0]->to_s,
                        $footprints[$_]->to_s(%params)
                    );

                    ($footprints[$_-1]->rg_ends_array)[0]->{footprint} = undef;
                    $footprints[$_]->add_rg_ends($footprints[$_-1]->rg_ends_array);
                    ($footprints[$_-1]->rg_ends_array)[0]->{footprint} = $footprints[$_];
                    $cluster->remove_footprint($footprints[$_-1]);
                }

                if (
                        # Merging with footprint downstream?
                        $_ < $#footprints  and
                        $footprints[$_+1]->size == 1  and
                        $footprints[$_+1]->first_pos - $footprints[$_]->last_pos < $MAX_DISTANCE  and
                        ($footprints[$_+1]->rg_ends_array)[0]->rg->is_pp_type  and
                        ($footprints[$_+1]->rg_ends_array)[0]->rg->id eq $rg_ends[0]->rg->id

                ) {
                    print STDERR sprintf(
                        "Merging solo footprint %s with footprint %s to rescue inversion shard\n",
                        ($footprints[$_+1]->rg_ends_array)[0]->to_s,
                        $footprints[$_]->to_s(%params)
                    );

                    ($footprints[$_+1]->rg_ends_array)[0]->{footprint} = undef;
                    $footprints[$_]->add_rg_ends($footprints[$_+1]->rg_ends_array);
                    ($footprints[$_+1]->rg_ends_array)[0]->{footprint} = $footprints[$_];
                    $cluster->remove_footprint($footprints[$_+1]);
                }

                $_++;
            }
        }
    }
}

sub normalise_sv_clustering {
    my $self = shift;
    my %params = @_;
    $self->remove_flanking_dels_and_tds(%params);
    $self->separate_out_shard_clusters(%params);
    $self->normalise_single_rg_clusters(%params);
    # $self->find_inv_local_two_jump_clusters(%params);  # Merge clusters with two inversion-type event in a <5Mb span
    $self->find_local_clusters(%params);  # Merge clusters contained in a <5Mb span
    $self->rescue_shards_and_balanced_bkpt_footprints(%params);  # Rescue shard or balanced breakpoint footprints
    $self->rescue_inversion_shards(%params);
    $self->normalise_simple_dels_and_tds_clusters(%params);
}

#
# End of functions for improving SV clusters
#


# Function for marking balanced rearrangement pairs.
sub mark_reciprocal_rearrangements {
    my $self = shift;
    my %params = @_;
    my @rearrangements = values %{$self->rg_of_id};
    my $i;
    my @reciprocal_rg_ends;
    for $i (0..($#rearrangements-1)) {
        ## Start the reciprocal rearrangement search from the low end.
        if ($rearrangements[$i]->low_end->is_fwd) {
            @reciprocal_rg_ends = $rearrangements[$i]->low_end->get_us_effector_rg_ends(%params);
        }
        else {
            @reciprocal_rg_ends = $rearrangements[$i]->low_end->get_ds_effector_rg_ends(%params);
        }
        for (@reciprocal_rg_ends) {
            if ($rearrangements[$i]->id ge $_->id) {
                next;
            }
            elsif ($rearrangements[$i]->is_reciprocal_with($_->rg)) {
                $rearrangements[$i]->add_reciprocal_rgs($_->rg);
                $_->rg->add_reciprocal_rgs($rearrangements[$i]);
            }
        }
    }
}

sub compute_connected_components {
    ## A breadth-first search algorithm for connected components of
    ## rearrangements. Currently rearrangements become connected only
    ## through being balanced. 
    my $self = shift;
    my %params = @_;

    my %unseen_rearrangements = %{$self->rg_of_id};
    my $current_component;
    my @connected_components = ();  ## A list of all connected components. 
    my @find_neighbours_of = ();    ## A helper list for keeping tab of unanalysed rearrangements
                                    ## in the current component. FIFO. 
    my ($cur_rg_id, $cur_rg, $partner_rg_end);

    for my $rg_id (keys %unseen_rearrangements) {
        # Has this been already 'seen' in previous iterations?
        if (!exists($unseen_rearrangements{$rg_id})) {
            next;
        }
        $cur_rg = $unseen_rearrangements{$rg_id};

        # Start breadth-first search with $rg_id. First initiate a new 'current component'. 
        $current_component = RearrangementConnectedComponent->new(
            rgs => { $rg_id => $cur_rg }
        );
        die if defined($cur_rg->component);
        $cur_rg->set_component($current_component);
        @find_neighbours_of = ($rg_id);
        delete $unseen_rearrangements{$rg_id};

        # Then iterate through the breadth-first search. 
        while (@find_neighbours_of) {
            $cur_rg_id = shift @find_neighbours_of;
            $cur_rg    = $current_component->rgs->{$cur_rg_id};

            # Find potential neighbours for low end...
            $partner_rg_end = $cur_rg->low_end->balanced_bkpt_partner_rg_end(%params);
            if (
                defined($partner_rg_end) &&
                exists($unseen_rearrangements{$partner_rg_end->id})
            ) {
                push @find_neighbours_of, $partner_rg_end->id;
                delete $unseen_rearrangements{$partner_rg_end->id};
                $current_component->add_rgs($partner_rg_end->rg);
                die if defined($partner_rg_end->rg->component);
                $partner_rg_end->rg->set_component($current_component);
            }

            # ... and high end.
            $partner_rg_end = $cur_rg->high_end->balanced_bkpt_partner_rg_end(%params);
            if (
                defined($partner_rg_end) &&
                exists($unseen_rearrangements{$partner_rg_end->id})
            ) {
                push @find_neighbours_of, $partner_rg_end->id;
                delete $unseen_rearrangements{$partner_rg_end->id};
                $current_component->add_rgs($partner_rg_end->rg);
                die if defined($partner_rg_end->rg->component);
                $partner_rg_end->rg->set_component($current_component);
            }
        }

        # Store the current component and prepare for the next available rearrangement. 
        push @connected_components, $current_component;
    }

    $self->{connected_components} = \@connected_components;
}

sub print_rearrangements_in_bedpe {
    my $self = shift;
    my %params = @_;
    my $rg;
    for (keys %{$self->rg_of_id}) {
        $rg = $self->rg_of_id->{$_};
        print join(
            "\t",
            @{$rg->orig_data}[0..5],
            $rg->id,
            $rg->orig_data->[7],
            $rg->low_end->dir,
            $rg->high_end->dir,
            @{$rg->orig_data}[10..$#{$rg->orig_data}],
            ($rg->low_end->cn_across_bkpt(%params) || "NA"),
            ($rg->high_end->cn_across_bkpt(%params) || "NA"),
            $rg->low_end->segment->cn_s,
            $rg->low_end->segment_end_neighbour->cn_s,
            $rg->high_end->segment->cn_s,
            $rg->high_end->segment_end_neighbour->cn_s,
        ) . "\n";
    }
}

sub print_rg_cns_bedpe {
    my $self = shift;
    my %params = @_;
    for my $chr (sort {$a->name cmp $b->name} values($self->chrs)) {
        $chr->print_rg_cns_bedpe;
    }
}

sub print_naive_classifications {
    # Old implementation that prints classifications based on $self->{complex_events}
    #

    my $self = shift;
    my %params = @_;
    my ($c, $classification, $event_id);
    $event_id = 0;
    
    my %components_in_complex_regions = ();
    my($event, $component);
    my($rg1, $rg2);
    for $event (@{$self->{complex_events}}) {
        if ($event->components_array == 1 && (($event->components_array)[0])->rgs_array == 2) {
            ($rg1, $rg2) = (($event->components_array)[0])->rgs_array;
            if (
                $rg1->is_part_of_shard_cycle(%params) &&
                $rg1->low_end->shard_partner(%params)->rg->id eq $rg2->id
            ) { 
                print join(
                    "\t",
                    $rg1->low_end->chr_name,
                    $rg1->low_end->pos - 1,
                    $rg1->low_end->pos,
                    $rg1->high_end->chr_name,
                    $rg1->high_end->pos - 1,
                    $rg1->high_end->pos,
                    (($params{sample_name} ? "$params{sample_name}:" : "") . "$event_id:" . $rg1->id),
                    1,
                    $rg1->low_end->dir,
                    $rg1->high_end->dir,
                    "two-jump",
                    (
                        "(" . $rg1->low_end->segment->cn_s  . "," . $rg1->low_end->segment_end_neighbour->cn_s . ")/" .
                        "(" . $rg1->high_end->segment->cn_s . "," . $rg1->high_end->segment_end_neighbour->cn_s . ")"
                    ),
                ) . "\n";
                
                print join(
                    "\t",
                    $rg2->low_end->chr_name,
                    $rg2->low_end->pos - 1,
                    $rg2->low_end->pos,
                    $rg2->high_end->chr_name,
                    $rg2->high_end->pos - 1,
                    $rg2->high_end->pos,
                    (($params{sample_name} ? "$params{sample_name}:" : "") . "$event_id:" . $rg2->id),
                    1,
                    $rg2->low_end->dir,
                    $rg2->high_end->dir,
                    "two-jump",
                    (
                        "(" . $rg2->low_end->segment->cn_s  . "," . $rg2->low_end->segment_end_neighbour->cn_s . ")/" .
                        "(" . $rg2->high_end->segment->cn_s . "," . $rg2->high_end->segment_end_neighbour->cn_s . ")"
                    ),
                ) . "\n";
                
                $event_id++;
                next;
            }  
        }

        for $component ($event->components_array) {
            for ($component->rgs_array) {
                print join(
                    "\t",
                    $_->low_end->chr_name,
                    $_->low_end->pos - 1,
                    $_->low_end->pos,
                    $_->high_end->chr_name,
                    $_->high_end->pos - 1,
                    $_->high_end->pos,
                    (($params{sample_name} ? "$params{sample_name}:" : "") . "$event_id:" . $_->id),
                    1,
                    $_->low_end->dir,
                    $_->high_end->dir,
                    "chromothripsis",
                    (
                        "(" . $_->low_end->segment->cn_s  . "," . $_->low_end->segment_end_neighbour->cn_s . ")/" .
                        "(" . $_->high_end->segment->cn_s . "," . $_->high_end->segment_end_neighbour->cn_s . ")"
                    ),
                ) . "\n";
            }

            $components_in_complex_regions{$component} = 1;
        }

        $event_id++;
    }

    for my $c ($self->connected_components_array) {
        if (exists($components_in_complex_regions{$c})) {
            next;
        }

        $classification = $c->naive_classification(%params);
        for ($c->rgs_array) {
            print join(
                "\t",
                $_->low_end->chr_name,
                $_->low_end->pos - 1,
                $_->low_end->pos,
                $_->high_end->chr_name,
                $_->high_end->pos - 1,
                $_->high_end->pos,
                (($params{sample_name} ? "$params{sample_name}:" : "") . "$event_id:" . $_->id),
                1,
                $_->low_end->dir,
                $_->high_end->dir,
                $classification,
                (
                    "(" . $_->low_end->segment->cn_s  . "," . $_->low_end->segment_end_neighbour->cn_s . ")/" .
                    "(" . $_->high_end->segment->cn_s . "," . $_->high_end->segment_end_neighbour->cn_s . ")"
                ),
            ) . "\n";
        }

        $event_id++;
    }
}

sub print_classification_data {
    # Some old implementation as well...
    my $self = shift;
    my %params = @_;
    my ($c, $classification, $event_id);
    my($rg1, $rg2, @cur_rgs);

    sub print_rgs_of_event {
        my $sample_name = shift;
        my $eid = shift;
        my $etype = shift;
        my %params = %{shift()};
        for my $rg (@_) {
            print join(
                "\t",
                $rg->low_end->chr_name,
                $rg->low_end->pos - 1,
                $rg->low_end->pos,
                $rg->high_end->chr_name,
                $rg->high_end->pos - 1,
                $rg->high_end->pos,
                (($sample_name ? "$sample_name:" : "") . "$eid:" . $rg->id),
                1,
                $rg->low_end->dir,
                $rg->high_end->dir,
                $etype,
                (  # Copy numbers at each SV
                    "(" . $rg->low_end->segment->cn_s  . "," . $rg->low_end->segment_end_neighbour->cn_s . ")/" .
                    "(" . $rg->high_end->segment->cn_s . "," . $rg->high_end->segment_end_neighbour->cn_s . ")"
                ),
                $rg->rg_type_s(%params),
                $rg->rg_dist
            ) . "\n";
        }
    }

    my $sample_name = ($params{sample_name} || "");
    my ($cur_component, $event_type);
    for $event_id (keys %{$self->{rgs_of_event}}) {
        @cur_rgs = values %{$self->{rgs_of_event}->{$event_id}};
        $cur_component = RearrangementConnectedComponent->new(%params, rgs => $self->{rgs_of_event}->{$event_id});
        $event_type = $cur_component->naive_classification(%params);
        print_rgs_of_event($sample_name, $event_id, $event_type, \%params, @cur_rgs);
    }
}

sub print_unique_breakpoints {
    my $self = shift;
    my %params = @_;

    if (!$params{lr_classification}) {
        die "File for -lr_classification is required for $0\->print_unique_breakpoints!";
    }

    open IN, $params{lr_classification} or die $!;
    my @F;
    while (<IN>) {
        chomp;
        @F = split /\t/;
        $self->rg_of_id->{$F[0]}->{is_lr} = ($F[1] eq "lr");
    }
    close IN;

    my @out_bkpts;
    my $i;
    my @involved_rgs;
    for my $chr (sort { $a->name cmp $b->name } $self->chrs_array) {
        @out_bkpts = $chr->get_unique_breakpoints(%params);
        $i = 0;
        while ($i <= $#out_bkpts) {
            chomp(@F = split /\t/, $out_bkpts[$i]);
            @involved_rgs = split /,/, $F[2];
            @involved_rgs = grep { !($self->rg_of_id->{substr($_, -1)}->{is_lr}) } @involved_rgs;
            if (@involved_rgs) {
                $out_bkpts[$i] = "$F[0]\t$F[1]\t" . join(",", @involved_rgs) . "\n";
                $i++;
            }
            else {
                splice(@out_bkpts, $i);
            }
        }

        print(@out_bkpts);
    }
}


sub find_rg_pattern_motifs {
    my $self = shift;
    my %params = @_;
    if (!exists($params{within})) {
        die "Attempted to call $self\->find_rg_pattern_motifs() without a parameter 'within'";
    }
    if (!exists($params{away_from})) {
        die "Attempted to call $self\->find_rg_pattern_motifs() without a parameter 'away_from'";
    }
    if (!exists($params{size})) {
        die "Attempted to call $self\->find_rg_pattern_motifs() without a parameter 'size'";
    }

    # Strategy:
    # Store found motifs in a hash.
    # Then go through each rearrangement in turn.
    # Find neighbours of each rearrangment, and grow the size of each group
    # until the desired size is found. 
    my %selected_groups = ();
    my %selected_rgs;
    my $chr;
    my $seg;
    my @rg_ends;
    my $group;
    my($rg_end, $rg);
    my $have_non_involved_rgs_nearby;
    my @neighbour_candidates;

    sub rg_neighbours {
        my $rg = shift;
        my %params = shift;
        my %neighbours;
        for ($rg->low_end->closest_rg_end(%params), $rg->high_end->closest_rg_end(%params)) {
            $neighbours{$_->rg->id} = $_->rg;
        }
        return values(%neighbours);
    }

    sub iterate_further_rgs {
        my %included_rgs = @{shift()};
        my $next_rg_to_be_included = shift;
        my $depth = shift;
        my %params = @_;
        my $rg;

        $included_rgs{$next_rg_to_be_included->id} = $next_rg_to_be_included;
        if ($depth == 1) {
            # Check if any rearrangement end of the current included rearrangements are
            # too close to a non-included rearrangement. 
            for $rg (values %included_rgs) {
                if (
                    defined($rg->low_end->closest_rg_end(%params, exclude_rgs => [keys %included_rgs])) ||
                    defined($rg->high_end->closest_rg_end(%params, exclude_rgs => [keys %included_rgs]))
                ) {
                    return();
                }
            }

            return(
                join(",", sort(keys %included_rgs)),
                RearrangementGroup->new(
                    components => [map { $_->component } values(%included_rgs)]
                )
            );
        }
        else {
            # Find new potential neighbours and iterate further. 
            my @iterated_groups = ();
            my %neighbours = ();

            for $rg (values %included_rgs) {
                for (
                    $rg->low_end->closest_rg_end(%params, exclude_rgs => [keys %included_rgs]),
                    $rg->high_end->closest_rg_end(%params, exclude_rgs => [keys %included_rgs])
                ) {
                    if (defined($_)) {
                        $neighbours{$_->rg->id} = $_->rg;
                    }
                }
            }

            for $rg (values %neighbours) {
                push @iterated_groups, iterate_further_rgs([%included_rgs], $rg, $depth-1, %params);
            }

            return @iterated_groups;
        }
    }

    for (values $self->rg_of_id) {
        %selected_groups = (%selected_groups, iterate_further_rgs([], $_, $params{size}, %params));
    }

#     for $chr ($self->chrs_array) {
#         @rg_ends = ();
#         for $seg ($chr->segments_array) {
#             if (defined($seg->low_end->bkpt)) {
#                 push @rg_ends, $seg->low_end->bkpt;
#             }
#             if (defined($seg->high_end->bkpt)) {
#                 push @rg_ends, $seg->high_end->bkpt;
#             }
#         }
#         
#         # Strategy for each encountered rearrangement end:
#         # 1. Collect an array of $params{size} successive rearrangements
#         # 2. Check that none of the involved rearrangements ends are within
#         #    $params{away_from} base pairs of a non-involved rearrangement end. 
#         # 3. Make sure rearrangements are within a window of $params{within} size
#         my $next_idx;
#         GROUP_CANDIDATE:
#         for (0..($#rg_ends - $params{size})) {
#             %selected_rgs = ();
#             $next_idx = $_-1;
#             while (scalar(keys %selected_rgs) < $params{size} && $next_idx+1 < $#rg_ends) {
#                 $next_idx++;
#                 if ($rg_ends[$next_idx]->segment_end->pos - $rg_ends[$_]->segment_end->pos > $params{within}) {
#                     next GROUP_CANDIDATE;
#                 }
#                 $selected_rgs{$rg_ends[$next_idx]->rg->id} = $rg_ends[$next_idx]->rg;
#             }
#             if (scalar(keys %selected_rgs) < $params{size}) {
#                 next;
#             }
# 
#             # Check if the current group is already handled/found
#             for $rg_end (@rg_ends[$_..($_+$params{size}-1)]) {
#                 $selected_rgs{$rg_end->rg->id} = $rg_end->rg;
#             }
#             if (exists($selected_groups{ join(',', sort(keys %selected_rgs)) })) {
#                 next;
#             }
# 
#             # Make sure none of the involved rearrangements are within near
#             # non-involved rearrangement ends.
#             $have_non_involved_rgs_nearby = 0;
#             for $rg (values %selected_rgs) {
#                 if (
#                     defined($rg->low_end->closest_rg_end(within => $params{away_from}, exclude_rgs => [keys %selected_rgs])) ||
#                     defined($rg->high_end->closest_rg_end(within => $params{away_from}, exclude_rgs => [keys %selected_rgs]))
#                 ) {
#                     $have_non_involved_rgs_nearby = 1;
#                     last;
#                 }
#             }
#             if ($have_non_involved_rgs_nearby) {
#                 next;
#             }
# 
#             # Finally store the group
#             $selected_groups{join(',', sort(keys %selected_rgs))} = RearrangementGroup->new(
#                 components => [map { $_->component } values(%selected_rgs)]
#             );
# 
#             # Sanity check - should never happen if $params{within} < $params{away_from}
#             if (scalar($selected_groups{join(',', sort(keys %selected_rgs))}->rgs_array) > $params{size}) {
#                 die;
#             }
#         }
#     }
# 

    map { $_->get_normalised_rg_patterns(%params) } values(%selected_groups);
    return values(%selected_groups);
}

sub print_footprint_cluster_classifications {
    # Classification function based on the footprints/SV clustering-based method.
    # For each SV cluster, print the cluster and footprint classifications. 

    my $self = shift;
    my %params = @_;
    my($cluster_type, $footprint_type, @rgs, $rg);

    sub rg_footprint {
        # Helper function to return the lexicographically smaller footprint
        # ID of a rearrangement.
        my $rg = $_[0];
        my $id_low = $rg->low_end->footprint->id;
        my $id_high = $rg->high_end->footprint->id;
        return($id_low lt $id_high ? $id_low : $id_high);
    }

    sub add_sample_name {
        my $string = shift;
        my %params = @_;
        if (!exists($params{sample_name})) {
            return $string;
        }
        else {
            return $params{sample_name} . ":" . $string;
        }
    }

    for my $cluster (sort { $a->id cmp $b->id } values($self->{rg_clusters})) {
        # Sanity check
        if ($cluster->size == 0) {
            die "Empty cluster " . $cluster->id . " not deleted!";
        }
        if (grep {$_->size == 0} $cluster->footprints_array) {
            die "Cluster " . $cluster->id . " has an empty footprint!";
        }

        $cluster_type = $cluster->cluster_classification(%params);
        @rgs = $cluster->rgs_array;
        @rgs = sort { rg_footprint($a) cmp rg_footprint($b) } @rgs;  # Order by footprint ID
        for $rg (@rgs) {
            my($low_end_footprint_rg_string, $high_end_footprint_rg_string);

            my $low_end_rg_pattern = $rg->low_end->footprint->rg_pattern(%params, symbol_seg_names => 1);
            if ($rg->low_end->footprint->size <= 6  and  defined($low_end_rg_pattern)) {
                $low_end_footprint_rg_string = 
                    $low_end_rg_pattern->{rg_string} . " " .
                    join(
                        "/",
                        map(
                            {$_ eq "NA" ? "NA" : sprintf "%.1f", $_}
                            @{$low_end_rg_pattern->{cn_changes}}
                        )
                    );
            }
            else {
                $low_end_footprint_rg_string = "NA";
            }

            my $high_end_rg_pattern = $rg->high_end->footprint->rg_pattern(%params, symbol_seg_names => 1);
            if ($rg->high_end->footprint->size <= 6  and  defined($high_end_rg_pattern)) {
                $high_end_footprint_rg_string = 
                    $high_end_rg_pattern->{rg_string} . " " .
                    join(
                        "/",
                        map(
                            {$_ eq "NA" ? "NA" : sprintf "%.1f", $_}
                            @{$high_end_rg_pattern->{cn_changes}}
                        )
                    )
            }
            else {
                $high_end_footprint_rg_string = "NA";
            }

            print join(
                "\t",
                @{$rg->orig_data}[0..5],
                add_sample_name($rg->orig_data->[6], %params),
                @{$rg->orig_data}[7..11],
                add_sample_name($cluster->id, %params),
                $cluster->n_rgs,
                $cluster->size,
                add_sample_name($rg->low_end->footprint->id, %params),
                add_sample_name($rg->high_end->footprint->id, %params),
                $cluster_type,
                $rg->low_end->footprint->detailed_type(%params),
                $rg->high_end->footprint->detailed_type(%params),
                $low_end_footprint_rg_string,
                $high_end_footprint_rg_string,
            ) . "\n";
        }
    }
}
#
# End of worker subroutines
#


1;
