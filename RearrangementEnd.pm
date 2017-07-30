# RearrangementEnd.pm

use warnings FATAL => 'all';
use strict;
use Scalar::Util qw(refaddr blessed);
use List::Util qw(min max);
use Footprint;

package RearrangementEnd;

sub new {
    my $class = shift;
    my %params = @_; 
    return bless {
        end            => $params{end},  # 'low' or 'high'
        dir            => $params{dir},
        segment_end    => ($params{segment_end} || undef),
        rg             => ($params{rg} || undef),
        min_reads_pos  => ($params{min_reads_pos} || undef),
        min_reads_clip => ($params{min_reads_clip} || undef),
        max_reads_pos  => ($params{max_reads_pos} || undef),
        max_reads_clip => ($params{max_reads_clip} || undef),
        footprint      => ($params{footprint} || undef),  # Footprint ID of the footprint in question
    }, $class;
}

sub copy {
    my $self = shift;
    return RearrangementEnd->new(%{$self});
}

#
# Getter methods
#
sub end {
    my $self = shift;
    if (!defined($self->{end})) {
        die "Attempted to call $self\->end() to get undefined $self\->{end}!";
    }
    if ($self->{end} ne "low" && $self->{end} ne "high") {
        die "Attempted to call $self\->end() with the value $self\->{end} eq '$self->{end}' - expected 'low' or 'high'";
    }
    return $self->{end};
}

sub dir {
    my $self = shift;
    if (!defined($self->{dir})) {
        die "Attempted to call dir() to get undefined $self\->{dir}!";
    }
    if ($self->{dir} ne "+" && $self->{dir} ne "-") {
        die "Attempted to call $self\->dir() but $self\->{dir} is neither '+' nor '-'!";
    }
    return $self->{dir};
}

sub is_l {
    my $self = shift;
    return $self->dir eq "low";
}

sub is_h {
    my $self = shift;
    return $self->dir eq "high";
}

sub is_fwd {
    my $self = shift;
    if ($self->dir ne '+' && $self->dir ne '-') {
        die "Called $self\->is_fwd() with $self\->{dir} not in ('+', '-')";
    }
    return $self->dir eq '+';
}

sub is_rev {
    my $self = shift;
    if ($self->dir ne '+' && $self->dir ne '-') {
        die "Called $self\->is_fwd() with $self\->{dir} not in ('+', '-')";
    }
    return $self->dir eq '-';
}

sub segment_end {
    my $self = shift;
    if (!defined($self->{segment_end})) {
        die "Attempted to call segment_end() to get undefined $self\->{segment_end}!";
    }
    return $self->{segment_end};
}

sub rg {
    my $self = shift;
    if (!defined($self->{rg})) {
        die "Attempted to call rg() to get undefined $self\->{rg}!";
    }
    return $self->{rg};
}

sub min_reads_pos {
    my $self = shift;
    if (!defined($self->{min_reads_pos})) {
        die "Attempted to call min_reads_pos() to get undefined $self\->{min_reads_pos}!";
    }
    return $self->{min_reads_pos};
}

sub min_reads_clip {
    my $self = shift;
    if (!defined($self->{min_reads_clip})) {
        die "Attempted to call min_reads_clip() to get undefined $self\->{min_reads_clip}!";
    }
    return $self->{min_reads_clip};
}

sub max_reads_pos {
    my $self = shift;
    if (!defined($self->{max_reads_pos})) {
        die "Attempted to call max_reads_pos() to get undefined $self\->{max_reads_pos}!";
    }
    return $self->{max_reads_pos};
}

sub max_reads_clip {
    my $self = shift;
    if (!defined($self->{max_reads_clip})) {
        die "Attempted to call max_reads_clip() to get undefined $self\->{max_reads_clip}!";
    }
    return $self->{max_reads_clip};
}

sub footprint {
    my $self = shift;
    return $self->{footprint};
}
#
# End of getter methods
#


#
# Helper subroutines
#
sub eq {
    my $self = shift;
    my $target = shift;
    return Scalar::Util::refaddr($self) == Scalar::Util::refaddr($target);
}

sub segment {
    my $self = shift;
    return $self->segment_end->segment;
}

sub chr {
    my $self = shift;
    return $self->segment->chr;
}

sub chr_name {
    my $self = shift;
    return $self->chr->name;
}

sub pos {
    my $self = shift;
    return $self->segment_end->pos;
}

sub segment_end_neighbour {
    my $self = shift;
    if ($self->segment_end->end eq "high") {
        return $self->segment->next_seg->low_end;
    }
    elsif ($self->segment_end->end eq "low") {
        return $self->segment->prev_seg->high_end;
    }
    else {
        die;
    }
}

sub id {
    my $self = shift;
    return $self->rg->id;
}

sub mate {
    my $self = shift;
    if ($self->end eq 'low') {
        return $self->rg->high_end;
    }
    else {
        return $self->rg->low_end;
    }
}

sub to_s {
    my $self = shift;
    return $self->chr_name . ":" . $self->pos . ":" . $self->dir;
}

sub print {
    my $self = shift;
    print $self->to_s . "\n";
}

sub is_foldback {
    my $self = shift;
    my %params = @_;
    return $self->rg->is_foldback(%params);
}

sub is_reciprocal_with {
    my $self = shift;
    my $target = shift;
    if (!defined($target) || Scalar::Util::blessed($target) ne 'RearrangementEnd') {
        die "A argument of type 'RearrangementEnd' is needed for $self\->is_reciprocal_with()";
    }
    
    if ($self->chr ne $target->chr) {
        return 0;
    }

    if ($self->is_fwd && $target->is_rev && $self->pos < $target->pos) {
        return 1;
    }

    if ($self->is_rev && $target->is_fwd && $target->pos < $self->pos) {
        return 1;
    }

    return 0;
}

sub is_on_shard {
    my $self = shift;
    my %params = @_;
    return $self->segment->is_shard(%params);
}

sub shard_partner {
    my $self = shift;
    my %params = @_;
    if (!$self->is_on_shard(%params)) {
        die "Attempted to call $self\->shard_partner() on a $self that is not on a shard";
    }

    if ($self->dir eq "+") {
        return $self->segment->low_end->bkpt;
    }
    else {
        return $self->segment->high_end->bkpt;
    }
}

sub traverse_across_shards {
    my $self = shift;
    my %params = @_;
    my $cur_rg_end = $self;
    my $shard_count = 0;
    while ($cur_rg_end->mate->is_on_shard(%params)) {
        $cur_rg_end = $cur_rg_end->mate->shard_partner(%params);
        $shard_count++;

        if ($cur_rg_end == $self) {
            ## If we get here then we got a cyclical shard sequence
            return(undef, $shard_count);
        }
    }

    return($cur_rg_end->mate, $shard_count);
}

sub can_be_in_cis_with {
    my $self = shift;
    my $target = shift;
    if (!defined($target) || Scalar::Util::blessed($target) ne 'RearrangementEnd') {
        die "A argument of type 'RearrangementEnd' is needed for $self\->can_be_in_cis_with()";
    }

    return(
        $self->chr eq $target->chr &&
        (
            ( $self->pos <= $target->pos && $self->is_rev && $target->is_fwd ) ||
            ( $self->pos >= $target->pos && $self->is_fwd && $target->is_rev )
        )
    );
}
#
# End of helper subroutines
#

#
# Worker subroutines
#
sub rg_side_cn {
    my $self = shift;
    my %params = @_;
    if ($self->is_bal_rg_overlap(%params)) {
        if ($self->is_fwd) {
            return $self->segment->prev_seg->cn;
        }
        else {
            return $self->segment->next_seg->cn;
        }
    }
    else {
        return $self->segment->cn;
    }
}

sub non_rg_side_cn {
    my $self = shift;
    my %params = @_;

    # Doesn't matter if the current rg end is an overlap balanced rearrangement or not. 
    if ($self->is_fwd) {
        return $self->segment->next_seg->cn;
    }
    else {
        return $self->segment->prev_seg->cn;
    }
}

sub cn_across_bkpt {
    my $self = shift;
    my %params = @_;
    my $rg_side_cn = $self->rg_side_cn(%params);
    my $non_rg_side_cn = $self->non_rg_side_cn(%params);
    if (!defined($rg_side_cn) || !defined($non_rg_side_cn)) {
        return undef;
    }

    return $rg_side_cn - $non_rg_side_cn;
}

sub relative_rg_side_cn_var {
    my $self = shift;
    my %params = @_;
    if (!$params{ploidy}) {
        die "Parameter 'ploidy' is needed for $self\-cn_change_across_rg()";
    }
    if (!$params{acf}) {
        die "Parameter 'acf' is needed for $self\-cn_change_across_rg()";
    }

    my $rg_side_cn = $self->rg_side_cn(%params);
    if (!defined($rg_side_cn)) {
        die "Attempted to call $self\->relative_cn_diff_stdev() with undefined $self\->rg_side_cn";
        return undef;
    }
    $rg_side_cn = 0 if $rg_side_cn < 0;
    my $adjusted_rg_side_cn = ( 2*(1-$params{acf}) + $rg_side_cn*$params{acf} ) /
                              ( 2*(1-$params{acf}) + $params{ploidy}*$params{acf} );
    
    return $adjusted_rg_side_cn / $self->segment->n_win;
}

sub relative_cn_diff_var {
    my $self = shift;
    my %params = @_;
    if (!$params{ploidy}) {
        die "Parameter 'ploidy' is needed for $self\-cn_change_across_rg()";
    }
    if (!$params{acf}) {
        die "Parameter 'acf' is needed for $self\-cn_change_across_rg()";
    }

    if ($self->is_bal_rg_overlap(%params)) {
        # We can't estimate the copy number because the breakpoint
        # is completely balanced (if not even more).
        return 1e9;
    }

    ## Fold-back rearrangements must be treated a bit differently
    if ($self->rg->is_foldback(%params)) {
        my $rg_side_cn;
        my $non_rg_side_cn;

        if ($self->is_fwd) {
            $rg_side_cn = $self->rg->low_end->rg_side_cn(%params);
        }
        else {
            $rg_side_cn = $self->rg->high_end->rg_side_cn(%params);
        }
        if (!defined($rg_side_cn)) {
            # die "Attempted to call $self\->relative_cn_diff_stdev() with undefined $self\->rg_side_cn";
            return undef;
        }
        $rg_side_cn = 0 if $rg_side_cn < 0;
        if ($self->is_fwd) {
            $non_rg_side_cn = $self->rg->high_end->non_rg_side_cn(%params);
        }
        else {
            $non_rg_side_cn = $self->rg->low_end->non_rg_side_cn(%params);
        }
        if (!defined($non_rg_side_cn)) {
            # die "Attempted to call $self\->relative_cn_diff_stdev() with undefined $self\->non_rg_side_cn";
            return undef;
        }
        $non_rg_side_cn = 0 if $non_rg_side_cn < 0;

        my $adjusted_rg_side_cn = ( 2*(1-$params{acf}) + $rg_side_cn*$params{acf} ) /
                                  ( 2*(1-$params{acf}) + $params{ploidy}*$params{acf} );
        my $adjusted_non_rg_side_cn = ( 2*(1-$params{acf}) + $non_rg_side_cn*$params{acf} ) /
                                      ( 2*(1-$params{acf}) + $params{ploidy}*$params{acf} );

        # Adjusted the expected read count variances by the length of the segments. 
        my $adjusted_cn_diff_var = 0;
        if ($self->is_fwd) {
            $adjusted_cn_diff_var = $adjusted_rg_side_cn / $self->rg->low_end->segment->n_win +
                                    $adjusted_non_rg_side_cn / $self->rg->high_end->segment_end_neighbour->segment->n_win;
        }
        else {
            $adjusted_cn_diff_var = $adjusted_rg_side_cn / $self->rg->high_end->segment->n_win +
                                    $adjusted_non_rg_side_cn / $self->rg->low_end->segment_end_neighbour->segment->n_win;
        }

        return $adjusted_cn_diff_var;
    }
    else {
        ## The adjusted copy numbers are proportional to the expected read counts,
        ## which is proportional to the variance of the expected read counts. 
        my $rg_side_cn = $self->rg_side_cn(%params);
        if (!defined($rg_side_cn)) {
            # die "Attempted to call $self\->relative_cn_diff_stdev() with undefined $self\->rg_side_cn";
            return undef;
        }
        $rg_side_cn = 0 if $rg_side_cn < 0;
        my $adjusted_rg_side_cn = ( 2*(1-$params{acf}) + $rg_side_cn*$params{acf} ) /
                                  ( 2*(1-$params{acf}) + $params{ploidy}*$params{acf} );
        my $non_rg_side_cn = $self->non_rg_side_cn(%params);
        if (!defined($non_rg_side_cn)) {
            # die "Attempted to call $self\->relative_cn_diff_stdev() with undefined $self\->non_rg_side_cn";
            return undef;
        }
        $non_rg_side_cn = 0 if $non_rg_side_cn < 0;
        my $adjusted_non_rg_side_cn = ( 2*(1-$params{acf}) + $non_rg_side_cn*$params{acf} ) /
                                      ( 2*(1-$params{acf}) + $params{ploidy}*$params{acf} );

        # Adjusted the expected read count variances by the length of the segments. 
        my $adjusted_cn_diff_var = $adjusted_rg_side_cn / $self->segment->n_win +
                                   $adjusted_non_rg_side_cn / $self->segment_end_neighbour->segment->n_win;
        return $adjusted_cn_diff_var;
    }
}

sub relative_cn_diff_stdev {
    my $self = shift;
    my %params = @_;
    return sqrt($self->relative_cn_diff_var(%params));
}

sub cn_relative_to_arm {
    my $self = shift;
    my %params = @_;
    if (defined($self->rg_side_cn(%params)) && defined($self->rg_side_arm_cn)) {
        $self->rg_side_cn(%params) - $self->rg_side_arm_cn;
    }   
    else {
        return undef;
    }   
}

sub balanced_bkpt_partner_rg_end {
    # Two rearrangements are defined as balanced if they
    # are at the same copy number breakpoint pointing
    # towards each other. 
    my $self = shift;
    my %params = @_;
    if ($self->is_bal_rg_overlap(%params)) {
        if ($self->is_fwd) {
            return $self->segment->low_end->bkpt;
        }
        else {
            return $self->segment->high_end->bkpt;
        }
    }
    elsif ($self->is_fwd && $self->segment->next_seg->is_bal_rg_gap(%params)) {
        return $self->segment->next_seg->next_seg->low_end->bkpt;
    }
    elsif ($self->is_rev && $self->segment->prev_seg->is_bal_rg_gap(%params)) {
        return $self->segment->prev_seg->prev_seg->high_end->bkpt;
    }
    else {
        return undef;
    }
}

sub is_bal_rg_overlap {
    my $self = shift;
    my %params = @_;
    return $self->segment->is_bal_rg_overlap(%params);
}

sub balanced_bkpt_partner_rg_end_by_footprint {
    # Two rearrangements are defined as balanced if they form a two-breakpoint
    # footprint and point towards each other. 

    my $self = shift;
    my %params = @_;
    if (!$self->footprint->is_balanced_type(%params)) {
        die "Attempted to use $self\->balanced_bkpt_partner_rg_end_by_footprint() "
          . "on a Footprint that is not a balanced breakpoint";
    }

    my($bkpt_1, $bkpt_2) = $self->footprint->rg_ends_array;
    if ($bkpt_1 == $self) {
        return $bkpt_2;
    }
    elsif ($bkpt_2 == $self) {
        return $bkpt_1;
    }
    else {
        die;
    }
}

sub shard_partner_rg_end_by_footprint {
    # Two rearrangements are defined as a shard if they form a two-breakpoint
    # footprint and point towards each other. 

    my $self = shift;
    my %params = @_;
    if (!$self->footprint->is_shard_type(%params)) {
        die "Attempted to use $self\->shard_bkpt_partner_rg_end_by_footprint() "
          . "on a Footprint that is not a shard breakpoint";
    }

    my($bkpt_1, $bkpt_2) = $self->footprint->rg_ends_array;
    if ($bkpt_1 == $self) {
        return $bkpt_2;
    }
    elsif ($bkpt_2 == $self) {
        return $bkpt_1;
    }
    else {
        die;
    }

}

sub has_balanced_bkpt_partner_rg_end {
    my $self = shift;
    my %params = @_;
    return defined($self->balanced_bkpt_partner_rg_end(%params));
}

sub get_us_effector_rg_ends {
    my $self = shift;
    my %params = @_;
    my @effector_rg_ends = ();
    
    if (
        $self->is_fwd &&
        defined($self->segment->low_end->bkpt) &&
        !$self->mate->eq($self->segment->low_end->bkpt)
    ) {
        push @effector_rg_ends, $self->segment->low_end->bkpt;
    }

    my $cur_seg = $self->segment;
    while (defined($cur_seg->prev_seg)) {
        $cur_seg = $cur_seg->prev_seg;
        if (
            defined($cur_seg->low_end->bkpt) &&
            !$self->mate->eq($cur_seg->low_end->bkpt) &&
            !$cur_seg->low_end->bkpt->rg->is_small_td(%params) &&
            !$cur_seg->low_end->bkpt->rg->is_small_del(%params) &&
            !$cur_seg->low_end->bkpt->rg->is_part_of_shard_cycle(%params)
        ) {
            push @effector_rg_ends, $cur_seg->low_end->bkpt;
        }
    }

    return @effector_rg_ends;
}

sub get_ds_effector_rg_ends {
    my $self = shift;
    my %params = @_;
    my @effector_rg_ends = ();
    
    if (
        $self->is_rev &&
        defined($self->segment->high_end->bkpt) &&
        !$self->mate->eq($self->segment->high_end->bkpt)
    ) {
        push @effector_rg_ends, $self->segment->high_end->bkpt;
    }

    my $cur_seg = $self->segment;
    while (defined($cur_seg->next_seg)) {
        $cur_seg = $cur_seg->next_seg;
        if (
            defined($cur_seg->high_end->bkpt) &&
            !$self->mate->eq($cur_seg->high_end->bkpt) &&
            !$cur_seg->high_end->bkpt->rg->is_small_td(%params) &&
            !$cur_seg->high_end->bkpt->rg->is_small_del(%params) &&
            !$cur_seg->high_end->bkpt->rg->is_part_of_shard_cycle(%params)
        ) {
            push @effector_rg_ends, $cur_seg->high_end->bkpt;
        }
    } 

    return @effector_rg_ends;
}

sub neighbours_array {
    # This subroutine looks for shard neighbours
    #

    my $self = shift;
    my %params = @_;
    my @neighbours = ();
    if ($self->is_fwd) {
        # Neighbour downstream of $self
        if (not $self->segment->is_bal_rg_overlap(%params)) {
            if (
                defined($self->segment->next_seg) &&
                !defined($self->segment->next_seg->low_end->bkpt) &&
                $self->segment->next_seg->length <= $params{shard_bypassing_slop} &&
                defined($self->segment->next_seg->high_end->bkpt)
            ) {
                push @neighbours, $self->segment->next_seg->high_end->bkpt;
            }
            elsif (
                defined($self->segment->next_seg) &&
                defined($self->segment->next_seg->next_seg) &&
                $self->segment->next_seg->next_seg->is_bal_rg_overlap(%params) &&
                $self->segment->next_seg->length <= $params{shard_bypassing_slop}
            ) {
                push @neighbours, $self->segment->next_seg->next_seg->high_end->bkpt;
            }
        }

        # Neighbour upstream of $self
        if (
            defined($self->segment->prev_seg) &&
            !defined($self->segment->low_end->bkpt) &&
            $self->segment->length <= $params{shard_bypassing_slop} &&
            defined($self->segment->prev_seg->high_end->bkpt) &&
            not $self->segment->prev_seg->is_bal_rg_overlap(%params)
        ) {
            push @neighbours, $self->segment->prev_seg->high_end->bkpt;
        }
    }
    else {
        # Neighbour downstream of $self
        if (not $self->segment->is_bal_rg_overlap(%params)) {
            if (
                defined($self->segment->prev_seg) &&
                !defined($self->segment->prev_seg->high_end->bkpt) &&
                $self->segment->prev_seg->length <= $params{shard_bypassing_slop} &&
                defined($self->segment->prev_seg->low_end->bkpt)
            ) {
                push @neighbours, $self->segment->prev_seg->low_end->bkpt;
            }
            elsif (
                defined($self->segment->prev_seg) &&
                defined($self->segment->prev_seg->prev_seg) &&
                $self->segment->prev_seg->prev_seg->is_bal_rg_overlap(%params) &&
                $self->segment->prev_seg->length <= $params{shard_bypassing_slop}
            ) {
                push @neighbours, $self->segment->prev_seg->prev_seg->low_end->bkpt;
            }
        }

        # Neighbour upstream of $self
        if (
            defined($self->segment->next_seg) &&
            !defined($self->segment->high_end->bkpt) &&
            $self->segment->length <= $params{shard_bypassing_slop} &&
            defined($self->segment->next_seg->low_end->bkpt) &&
            not $self->segment->next_seg->is_bal_rg_overlap(%params)
        ) {
            push @neighbours, $self->segment->next_seg->low_end->bkpt;
        }
    }
    return @neighbours;
}

sub closest_us_rg_end {
    my $self = shift;
    my %params = @_;
    if (!exists($params{within})) {
        die "Parameter 'within' is required for $self\->closest_us_rg_end()";
    }
    if ($params{within} < 0) {
        $params{within} = 999_999_999_999;
    }
    if (!exists($params{"exclude_rgs"})) {
        die "Parameter 'exclude_rgs' is required for $self\->closest_us_rg_end()";
        # $params{exclude_rgs} = [];
    }

    if (
        $self->is_fwd &&
        defined($self->segment->low_end->bkpt) &&
        !(grep { $self->segment->low_end->bkpt->rg->id eq $_ } @{$params{exclude_rgs}})
    ) {
        if ($self->pos - $self->segment->low_end->pos > $params{"within"}) {
            return undef;
        }
        return $self->segment->low_end->bkpt;
    }

    my $cur_seg = $self->segment;
    my $rg_end;
    while (1) {
        last if !defined($cur_seg->prev_seg);
        $cur_seg = $cur_seg->prev_seg;

        if (defined($cur_seg->high_end->bkpt)) {
            $rg_end = $cur_seg->high_end->bkpt;
            if (!(grep { $rg_end->rg->id eq $_ } @{$params{exclude_rgs}})) {
                if ($self->pos - $rg_end->pos > $params{"within"}) {
                    return undef;
                }
                return $rg_end;
            }
        }
        if (defined($cur_seg->low_end->bkpt)) {
            $rg_end = $cur_seg->low_end->bkpt;
            if (!(grep { $rg_end->rg->id eq $_ } @{$params{exclude_rgs}})) {
                if ($self->pos - $rg_end->pos > $params{"within"}) {
                    return undef;
                }
                return $rg_end;
            }
        }
    }

    return undef;
}

sub closest_ds_rg_end {
    my $self = shift;
    my %params = @_;
    if (!exists($params{"within"})) {
        die "Parameter 'within' is required for $self\->closest_ds_rg_end()";
    }
    if (!exists($params{"exclude_rgs"})) {
        die "Parameter 'exclude_rgs' is required for $self\->closest_ds_rg_end()";
        # $params{exclude_rgs} = [];
    }

    if (
        $self->is_rev &&
        defined($self->segment->high_end->bkpt) &&
        !(grep { $self->segment->high_end->bkpt->rg->id eq $_ } @{$params{exclude_rgs}})
    ) {
        if ($self->segment_end->segment->high_end->pos - $self->pos > $params{"within"}) {
            return undef;
        }
        return $self->segment->high_end->bkpt;
    }

    my $cur_seg = $self->segment;
    my $rg_end;
    while (1) {
        last if !defined($cur_seg->next_seg);
        $cur_seg = $cur_seg->next_seg;

        if (defined($cur_seg->low_end->bkpt)) {
            $rg_end = $cur_seg->low_end->bkpt;
            if (!(grep { $rg_end->rg->id eq $_ } @{$params{exclude_rgs}})) {
                if ($rg_end->pos - $self->pos > $params{"within"}) {
                    return undef;
                }
                return $rg_end;
            }
        }
        if (defined($cur_seg->high_end->bkpt)) {
            $rg_end = $cur_seg->high_end->bkpt;
            if (!(grep { $rg_end->rg->id eq $_ } @{$params{exclude_rgs}})) {
                if ($rg_end->pos - $self->pos > $params{"within"}) {
                    return undef;
                }
                return $rg_end;
            }
        }
    }

    return undef;
}


sub closest_rg_end {
    my $self = shift;
    my %params = @_;
    my $closest_us_rg_end = $self->closest_us_rg_end(%params);
    my $closest_ds_rg_end = $self->closest_ds_rg_end(%params);
    if (!defined($closest_us_rg_end) && !defined($closest_ds_rg_end)) {
        return undef;
    }
    elsif (!defined($closest_us_rg_end)) {
        return $closest_ds_rg_end;
    }
    elsif (!defined($closest_ds_rg_end)) {
        return $closest_us_rg_end;
    }
    else {
        if (
            $self->segment_end->pos - $closest_us_rg_end->segment_end->pos <=
            $closest_ds_rg_end->segment_end->pos - $self->segment_end->pos
        ) {
            return $closest_us_rg_end;
        }
        else {
            return $closest_ds_rg_end;
        }
    }
}

sub is_on_p_arm {
    my $self = shift;
    return $self->segment->is_on_p_arm;
}

sub is_on_q_arm {
    my $self = shift;
    return $self->segment->is_on_q_arm;
}

sub rg_side_arm_cn {
    # The CN of the closest rearrangement side telomere or centromere or the
    # breakpoint.

    my $self = shift;
    if ($self->is_on_p_arm) {
        if ($self->is_fwd) {
            return $self->chr->p_arm_telomere_cn;
        }   
        else {
            return $self->chr->p_centromere_cn;
        }   
    }   
    else {
        if ($self->is_fwd) {
            return $self->chr->q_centromere_cn;
        }   
        else {
            return $self->chr->q_arm_telomere_cn;
        }   
    }
}

sub non_rg_side_arm_cn {
    # The CN of the closest non-rearrangement telomere or centromere or the
    # breakpoint. 

    my $self = shift;
    if ($self->is_on_p_arm) {
        if ($self->is_fwd) {
            return $self->chr->p_centromere_cn;
        }   
        else {
            return $self->chr->p_arm_telomere_cn;
        }   
    }   
    else {
        if ($self->is_fwd) {
            return $self->chr->q_arm_telomere_cn;
        }   
        else {
            return $self->chr->q_centromere_cn;
        }   
    }
}

sub is_inter_chr {
    my $self = shift;
    return $self->chr != $self->mate->chr;
}

sub distance_to_mate {
    # Returns -1 if inter-chromosomal
    my $self = shift;
    if ($self->is_inter_chr) {
        return -1;
    }
    else {
        return abs($self->pos - $self->mate->pos);
    }
}
#
# End of worker subroutines
#


1;
