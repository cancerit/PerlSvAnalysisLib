use warnings FATAL => 'all';
use strict;

use File::Basename;
use lib dirname($0);
use RearrangementEnd;
use Scalar::Util qw(blessed);
use Carp;

package Rearrangement;

sub new {
    my $class = shift;
    my %params = @_;
    my $rg_low_end = RearrangementEnd->new(
        end => "low",
        dir => $params{low_end_dir},
        min_reads_pos => $params{reads_pos}->[0],
        min_reads_clip => $params{reads_clip}->[0],
        max_reads_pos => $params{reads_pos}->[1],
        max_reads_clip => $params{reads_clip}->[1],
    );
    my $rg_high_end = RearrangementEnd->new(
        end => "high",
        dir => $params{high_end_dir},
        min_reads_pos => $params{reads_pos}->[2],
        min_reads_clip => $params{reads_clip}->[2],
        max_reads_pos => $params{reads_pos}->[3],
        max_reads_clip => $params{reads_clip}->[3],
    );  

    my $rg = bless {
        low_end  => $rg_low_end,
        high_end => $rg_high_end,
        id => $params{id},
        reciprocal_rgs => {},
        orig_data => $params{orig_data},
    }, $class;
    $rg_low_end->{rg}  = $rg;
    $rg_high_end->{rg} = $rg;

    return $rg;
}

sub copy {
    my $self = shift;
    my $rg = bless {
        low_end => $self->low_end->copy,
        high_end => $self->high_end->copy,
        id => $self->id,
        reciprocal_rgs => $self->reciprocal_rgs,
        orig_data => $self->orig_data,
    };
    $rg->low_end->{rg} = $rg;
    $rg->high_end->{rg} = $rg;
    return $rg;
}

#
# getter methods
#
sub low_end {
    my $self = shift;
    if (!defined($self->{low_end})) {
        die "Attempted to call low_end() to get undefined Rearrangement->{low_end}!";
    }
    return $self->{low_end};
}

sub high_end {
    my $self = shift;
    if (!defined($self->{high_end})) {
        die "Attempted to call high_end() to get undefined Rearrangement->{high_end}!";
    }
    return $self->{high_end};
}

sub id {
    my $self = shift;
    if (!defined($self->{id})) {
        die "Attempted to call id() to get undefined Rearrangement->{id}!";
    }
    return $self->{id};
}

sub reciprocal_rgs {
    my $self = shift;
    return $self->{reciprocal_rgs};
}

sub reciprocal_rgs_array {
    my $self = shift;
    return (values %{$self->reciprocal_rgs});
}

sub reciprocal_rg_by_idx {
    my $self = shift;
    my $idx = shift;
    if (!defined($idx)) {
        return ($self->reciprocal_rgs_array)[0];
    }
    if ($idx < -1 or $idx > scalar($self->reciprocal_rgs_array) - 1) {
        die;
    }
    return ($self->reciprocal_rgs_array)[$idx];
}

sub cluster {
    my $self = shift;

    # Sanity check
    if ($self->low_end->footprint->cluster != $self->high_end->footprint->cluster) {
        die;
    }

    return $self->low_end->footprint->cluster;
}

sub orig_data {
    my $self = shift;
    return $self->{orig_data};
}
#
# End of getter methods
#


#
# Helper methods
#
sub to_s {
    my $self = shift;
    return $self->id . "\t" . $self->low_end->to_s . "\t" . $self->high_end->to_s . "\t" . (defined($self->component) ? $self->component : "Component: undef");
}

sub print {
    my $self = shift;
    print $self->to_s . "\n";
}

sub add_reciprocal_rgs {
    my $self = shift;

    for (@_) {
        if (!defined($_) || Scalar::Util::blessed($_) ne "Rearrangement") {
            die "A argument of type 'Rearrangement' is needed for $self\->add_reciprocal_rg()";
        }

        if (exists($self->reciprocal_rgs->{$_->id})) {
            Carp::carp "Attempted to add a pre-existing rearrangement in $self\->add_reciprocal_rg()";
        }
        else {
            $self->{reciprocal_rgs}->{$_->id} = $_;
        }
    }
}

sub set_component {
    my $self = shift;
    my $component = shift;
    if (!defined($component)) {
        die "Argument for component is needed for $self\->set_component()";
    }
    if (!defined(Scalar::Util::blessed $component) or Scalar::Util::blessed($component) ne "RearrangementConnectedComponent") {
        die "Argument must be of type 'RearrangementConnectedComponent' in $self\->set_component()";
    }
    $self->{component} = $component;
}

sub component {
    my $self = shift;
    if (!exists($self->{component})) {
        return undef;
    }
    return $self->{component};
}

sub is_foldback {
    my $self = shift;
    my %params = @_;
    if (!exists($params{max_foldback_distance})) {
        die "Parameter 'max_foldback_distance' is needed for $self\->is_foldback()";
    }
    
    if ($self->low_end->segment->chr != $self->high_end->segment->chr) {
        return 0;
    }

    if ($self->low_end->dir eq "+" && $self->high_end->dir eq "+") {
        if ($self->low_end->segment->is_bal_rg_overlap(%params)) {
            # Is a balanced rearrangement -> can't be fold-back
            return 0;
        }

        my $is_normal_fb =
            !defined($self->low_end->segment->next_seg->low_end->bkpt) &&
            defined($self->low_end->segment->next_seg->high_end->bkpt) &&
            $self->low_end->segment->next_seg->high_end->bkpt->eq($self->high_end) &&
            $self->low_end->segment->next_seg->length <= $params{max_foldback_distance};
        my $is_fb_with_overlap_bal =
            !defined($self->low_end->segment->next_seg->low_end->bkpt) &&
            !defined($self->low_end->segment->next_seg->high_end->bkpt) &&
            defined($self->low_end->segment->next_seg->next_seg) &&
            $self->low_end->segment->next_seg->next_seg->is_bal_rg_overlap(%params) &&
            defined($self->low_end->segment->next_seg->next_seg->high_end->bkpt) &&
            $self->low_end->segment->next_seg->next_seg->high_end->bkpt->eq($self->high_end) &&
            $self->low_end->segment->next_seg->length + $self->low_end->segment->next_seg->next_seg->length <= $params{max_foldback_distance};
        return($is_normal_fb || $is_fb_with_overlap_bal);
    }
    elsif ($self->low_end->dir eq "-" && $self->high_end->dir eq "-") {
        if ($self->high_end->segment->is_bal_rg_overlap(%params)) {
            # Is a balanced rearrangement -> can't be fold-back
            return 0;
        }

        my $is_normal_fb =
            !defined($self->high_end->segment->prev_seg->high_end->bkpt) &&
            defined($self->high_end->segment->prev_seg->low_end->bkpt) &&
            $self->high_end->segment->prev_seg->low_end->bkpt->eq($self->low_end) &&
            $self->high_end->segment->prev_seg->length <= $params{max_foldback_distance};
        my $is_fb_with_overlap_bal =
            !defined($self->high_end->segment->prev_seg->high_end->bkpt) &&
            !defined($self->high_end->segment->prev_seg->low_end->bkpt) &&
            defined($self->high_end->segment->prev_seg->prev_seg) &&
            $self->high_end->segment->prev_seg->prev_seg->is_bal_rg_overlap(%params) &&
            defined($self->high_end->segment->prev_seg->prev_seg->low_end->bkpt) &&
            $self->high_end->segment->prev_seg->prev_seg->low_end->bkpt->eq($self->low_end) &&
            $self->high_end->segment->prev_seg->length + $self->high_end->segment->prev_seg->prev_seg->length <= $params{max_foldback_distance};
        return($is_normal_fb || $is_fb_with_overlap_bal);
    }
    else {
        return 0;
    }
}

sub is_intra_chr {
    my $self = shift;
    return $self->low_end->chr_name eq $self->high_end->chr_name;
}

sub is_inter_chr {
    my $self = shift;
    return !$self->is_intra_chr;
}

sub is_del_type {
    my $self = shift;
    return(
        $self->is_intra_chr &&
        $self->low_end->is_fwd &&
        $self->high_end->is_rev
    );
}

sub is_td_type {
    my $self = shift;
    return(
        $self->is_intra_chr &&
        $self->low_end->is_rev &&
        $self->high_end->is_fwd
    );
}

sub is_pp_type {
    my $self = shift;
    return(
        $self->is_intra_chr &&
        $self->low_end->is_fwd &&
        $self->high_end->is_fwd
    );
}

sub is_mm_type {
    my $self = shift;
    return(
        $self->is_intra_chr &&
        $self->low_end->is_rev &&
        $self->high_end->is_rev
    );
}

sub is_inv_type {
    my $self = shift;
    return($self->is_pp_type || $self->is_mm_type);
}

sub rg_type_s {
    my $self = shift;
    my %params = @_;
    if ($self->is_inter_chr) {
        return "inter";
    }
    else {
        if    ($self->is_del_type) { return "del"; }
        elsif ($self->is_td_type)  { return "td"; }
        elsif ($self->is_foldback(%params)) { return "fb"; }
        elsif ($self->is_inv_type) { return "inv"; }
        else { die; }
    }
}

sub rg_dist {
    my $self = shift;
    if ($self->is_inter_chr) {
        return "inter";
    }
    else {
        return abs($self->high_end->pos - $self->low_end->pos);
    }
}

sub is_small_del {
    my $self = shift;
    my %params = @_;
    if (!exists($params{min_seg_size_for_small_del_or_td})) {
        die "Parameter 'min_seg_size_for_small_del_or_td' is needed for $self\->is_small_del()";
    }

    if (
        $self->low_end->dir eq "+" &&
        defined($self->low_end->segment->next_seg) &&
        defined($self->low_end->segment->next_seg->next_seg) &&
        defined($self->low_end->segment->next_seg->next_seg->low_end->bkpt) &&
        $self->low_end->mate->eq($self->low_end->segment->next_seg->next_seg->low_end->bkpt) &&
        $self->low_end->segment->next_seg->length < $params{min_seg_size_for_small_del_or_td}
    ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub is_small_td {
    my $self = shift;
    my %params = @_;
    if (!exists($params{min_seg_size_for_small_del_or_td})) {
        die "Parameter 'min_seg_size_for_small_del_or_td' is needed for $self\-is_small_td()";
    }

    if (
        $self->low_end->segment->length < $params{min_seg_size_for_small_del_or_td} &&
        $self->low_end->dir eq "-" &&
        defined($self->low_end->segment->high_end->bkpt) &&
        $self->low_end->mate->eq($self->low_end->segment->high_end->bkpt)
    ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub is_shard_bypassing {
    my $self = shift;
    my %params = @_;

    ## Algorithm:
    ## For each rearrangement end neighbour for the index rearrangement low
    ## end, find a shard path that leads to the other end of the
    ## rearrangement. If found return 1. 
    my @target_rg_ends = $self->high_end->neighbours_array(%params);
    my @source_rg_ends = $self->low_end->neighbours_array(%params);
    my($src_rg_end, $cur_rg_end);
    my @bypassed_rgs;
    for my $src_rg_end (@source_rg_ends) {
        if (!$src_rg_end->mate->is_on_shard(%params)) {
            next;
        }
        $cur_rg_end = $src_rg_end;
        @bypassed_rgs = ($src_rg_end->rg);
        while ($cur_rg_end->mate->is_on_shard(%params)) {
            $cur_rg_end = $cur_rg_end->mate->shard_partner(%params);
            push @bypassed_rgs, $cur_rg_end->rg;

            ## Ignore shard cycles
            if ($cur_rg_end == $src_rg_end) {
                last;
            }

            ## Detect complete bridges
            if (grep { $_ == $cur_rg_end->mate } @target_rg_ends) {
                print(
                    STDERR
                    sprintf(
                        "REARRANGEMENT (%s) bypasses (%s)\n",
                        $self->id,
                        join(",", map {$_->id} @bypassed_rgs)
                    )
                );
                return 1;
            }
        }
    }

    return 0;
}

sub weighted_avg_cn_change_across_rg {
    my $self = shift;
    my %params = @_;

    # Fold-back rearrangements and "normal" rearrangements
    # are handled separately. 
    if ($self->is_foldback(%params)) {
        my($cn_before, $cn_after);
        if ($self->low_end->dir eq "+") {
            $cn_before = $self->low_end->segment->cn;
            $cn_after = $self->high_end->segment_end_neighbour->segment->cn;
        }
        else {
            $cn_before = $self->high_end->segment->cn;
            $cn_after = $self->low_end->segment_end_neighbour->segment->cn;
        }

        if (defined($cn_before) && defined($cn_after)) {
            return $cn_before - $cn_after;
        }
        else {
            return undef;
        }
    }

    # If we got this far, then the current rearrangement is not fold-back,
    # and we actually have to compute the weighted average CN change. 
    if (defined($self->low_end->cn_across_bkpt(%params)) && defined($self->high_end->cn_across_bkpt(%params))) {
        my $low_end_cn_change      = $self->low_end->cn_across_bkpt(%params);
        my $low_end_cn_change_var  = $self->low_end->relative_cn_diff_var(%params);
        my $high_end_cn_change     = $self->high_end->cn_across_bkpt(%params);
        my $high_end_cn_change_var = $self->high_end->relative_cn_diff_var(%params);
        my $weighted_cn_change = $low_end_cn_change/$low_end_cn_change_var + $high_end_cn_change/$high_end_cn_change_var;
        my $total_weight = 1/$low_end_cn_change_var + 1/$high_end_cn_change_var;
        return($weighted_cn_change / $total_weight);
    }
    elsif (defined($self->low_end->cn_across_bkpt(%params))) {
        return $self->low_end->cn_across_bkpt(%params);
    }
    elsif (defined($self->high_end->cn_across_bkpt(%params))) {
        return $self->high_end->cn_across_bkpt(%params);
    }
    else {
        return undef;
    }
}

sub weighted_avg_cn_change_var_across_rg {
    # Copy number change across the SV breakpoint - weighted by the lengths
    # of the associated segments at each end of the SV. 

    my $self = shift;
    my %params = @_;

    if ($self->is_foldback(%params)) {
        return $self->low_end->relative_cn_diff_var(%params);
    }
    else {
        my $low_end_var = $self->low_end->relative_cn_diff_var(%params);
        my $high_end_var = $self->high_end->relative_cn_diff_var(%params);
        if (!defined($low_end_var) && !defined($high_end_var)) {
            return undef;
        }
        elsif (!defined($low_end_var)) {
            return $high_end_var;
        }
        elsif (!defined($high_end_var)) {
            return $low_end_var;
        }
        else {
            return( 1 / ( 1/$low_end_var + 1/$high_end_var ) );
        }
    }
}

sub is_reciprocal_with {
    my $self = shift;
    my $target = shift;
    if (!defined($target) || Scalar::Util::blessed($target) ne "Rearrangement") {
        die "A argument of type 'Rearrangement' is needed for $self\->is_reciprocal_with()";
    }

    # Low end vs. low end
    if (
        $self->low_end->is_reciprocal_with($target->low_end) &&
        $self->high_end->is_reciprocal_with($target->high_end)
    ) {
        return 1;
    }
    if (
        $self->high_end->is_reciprocal_with($target->low_end) &&
        $self->low_end->is_reciprocal_with($target->high_end)
    ) {
        return 1;
    }

    return 0;
}

sub is_part_of_shard_cycle {
    my $self = shift;
    my %params = @_;
    my($mate, $shard_count) = $self->low_end->traverse_across_shards(%params);

    ## Returns number of shards if is part of shard cycle. Otherwise
    ## returns 0. 
    if (!defined($mate)) {
        return $shard_count;
    }
    else {
        return 0;
    }
}

sub is_long_range {
    my $self = shift;
    if ($self->low_end->chr ne $self->high_end->chr) {
        return 1;
    }
    elsif ($self->low_end->dir eq $self->high_end->dir) {
        return 1;
    }
    elsif ($self->high_end->pos - $self->low_end->pos > 10e6) {
        return 1;
    }
    else {
        return 0;
    }
}
#
# End of helper methods
#


1;
