use warnings FATAL => 'all';
use strict;
use Scalar::Util qw(blessed);
use List::Util qw(sum);

package CopyNumberSegmentArray;

sub new {
    my $class = shift;
    my %params = @_;

    return bless {
        name => $params{name},
        is_sorted => 0,
        is_inverted => ($params{is_inverted} or 0),
        segments => [],
    }, $class;
}

#
# Getter methods
#
sub name {
    my $self = shift;
    if (!defined($self->{name})) {
        die "Attempted to call name() to get undefined $self\->{name}!";
    }
    return $self->{name};
}

sub is_sorted {
    my $self = shift;
    if (!defined($self->{is_sorted})) {
        die "Attempted to call is_sorted() to get undefined $self\->{is_sorted}!";
    }
    return $self->{is_sorted};
}

sub segments {
    my $self = shift;
    if (!defined($self->{segments})) {
        die "Attempted to call segments() to get undefined $self\->{segments}!";
    }
    return $self->{segments};
}

sub segments_array {
    my $self = shift;
    return @{$self->segments};
}

sub size {
    my $self = shift;
    return scalar(@{$self->segments});
}

sub print {
    my $self = shift;
    my $start = ($_[0] or 0);
    my $end   = ($_[1] or $#{$self->segments});
    my $i = $start;
    print "[ Chromosome " . $self->name . "]\n";
    while ($i <= $end) {
        $self->segments->[$i]->print;
        $i++;
    }
    print "\n";
}
#
# End of getter methods
#

#
# Worker subroutines
#
sub add_copy_number_segment {
    my $self = shift;
    my $new_seg = shift;
    if (!defined(Scalar::Util::blessed $new_seg) or Scalar::Util::blessed($new_seg) ne "CopyNumberSegment") {
        die "Attempted to call $self\->add_copy_number_segment() without an argument of class 'CopyNumberSegment'";
    }
    push @{$self->segments}, $new_seg;
}

sub sort_segments {
    my $self = shift;
    
    @{$self->{segments}} = sort {
        $a->low_end->pos <=> $b->low_end->pos ||
        $a->high_end->pos <=> $b->high_end->pos
    } @{$self->segments};

    for (1..$#{$self->{segments}}) {
        $self->{segments}->[$_-1]->{next_seg} = $self->segments->[$_];
        $self->{segments}->[$_]->{prev_seg}   = $self->segments->[$_-1];
    }

    $self->{is_sorted} = 1;
}

sub bkpts_array {
    my $self = shift;
    my($seg, @bkpts);
    for $seg ($self->segments_array) {
        push @bkpts, $seg->low_end->bkpt if defined($seg->low_end->bkpt);
        push @bkpts, $seg->high_end->bkpt if defined($seg->high_end->bkpt);
    }
    return (@bkpts);
}

sub splice_segment {
    my $self = shift;
    my $index = shift;
    my $cur_seg = $self->segments->[$index];

    if (
        $cur_seg->low_end->is_bounded_by_rg ||
        $cur_seg->high_end->is_bounded_by_rg
    ) {
        die "Attempted to splice out a CopyNumberSegment bounded by rearrangements at $self\->splice_segment()";
    }

    if ($index > 0 && $index < $#{$self->segments}) {
        $cur_seg->prev_seg->high_end->{boundary} = $cur_seg->high_end->boundary;
        $cur_seg->next_seg->low_end->{boundary}  = $cur_seg->low_end->boundary;
        $cur_seg->prev_seg->{next_seg} = $cur_seg->next_seg;
        $cur_seg->next_seg->{prev_seg} = $cur_seg->prev_seg;
    }
    elsif ($index > 0) {
        $cur_seg->prev_seg->high_end->{boundary} = $cur_seg->high_end->boundary;
        $cur_seg->prev_seg->{next_seg} = undef;
    }
    elsif ($index < $#{$self->segments}) {
        $cur_seg->next_seg->low_end->{boundary} = $cur_seg->low_end->boundary;
        $cur_seg->next_seg->{prev_seg} = undef;
    }

    splice(@{$self->segments}, $index, 1);
}

sub remove_small_cn_bkpt_segments {
    # Removes small copy number segments that are bounded by copy number
    # segmentation breakpoints alone, regardless of any associated
    # copy number change.
    my $self = shift;
    my %params = @_;
    if (!defined($params{min_cn_bkpt_seg_size})) {
        die "Parameter 'min_cn_bkpt_seg_size' is required in $self\->remove_small_cn_bkpt_segments()";
    }
    if (!$self->is_sorted) {
        die "Attempted to call $self\->remove_small_cn_bkpt_segments() with an unsorted $self";
    }

    my $i;
    my $cur_seg;

    my $changed = 0;
    $i = 0;
    while ($i <= $#{$self->segments}) {
        $cur_seg = $self->segments->[$i];
        if (
            # !$cur_seg->low_end->is_bounded_by_rg &&
            # !$cur_seg->high_end->is_bounded_by_rg &&
            $cur_seg->low_end->is_bounded_by_cn_bkpt &&
            $cur_seg->high_end->is_bounded_by_cn_bkpt &&
            $cur_seg->length < $params{min_cn_bkpt_seg_size}
        ) {
            $self->splice_segment($i);
            $changed = 1;
        }
        else {
            $i++;
        }
    }

    return $changed;
}

sub merge_segment_with_next {
    my $self = shift;
    my $index = shift;

    if (!defined($self->segments->[$index]->next_seg)) {
        die "Attempted to call $self\->merge_segment_with_next() with undefined $self\->segments->[\$index]->{next_seg}";
    }
    if (defined($self->segments->[$index]->high_end->bkpt)) {
        die "Attempted to call $self\->merge_segment_with_next() through a rearrangement breakpoint";
    }

    my($total_len, $total_cn);
    my $cur_seg = $self->segments->[$index];
    if (!defined($cur_seg->cn) && !defined($cur_seg->next_seg->cn)) {
        $cur_seg->{cn} = undef;
    }
    elsif (!defined($cur_seg->cn)) {
        $cur_seg->{cn} = $cur_seg->next_seg->cn;
    }
    elsif (!defined($cur_seg->next_seg->cn)) {
        $cur_seg->{cn} = $cur_seg->cn;
    }
    else {
        $total_len = $cur_seg->length + $cur_seg->next_seg->length;
        $total_cn = $cur_seg->cn * $cur_seg->length +
                    $cur_seg->next_seg->cn * $cur_seg->next_seg->length;
        $cur_seg->{cn} = $total_cn / $total_len;
    }
    $cur_seg->{n_win} += $cur_seg->n_win + $cur_seg->next_seg->n_win - 1;

    $cur_seg->{high_end} = $cur_seg->next_seg->high_end;
    $cur_seg->high_end->{segment} = $cur_seg;
    if (defined($cur_seg->next_seg->high_end->bkpt)) {
        $cur_seg->next_seg->high_end->bkpt->{segment_end} = $cur_seg->high_end;
    }
    if (!$cur_seg->next_seg->high_end->is_bounded_by_tel) {
        $cur_seg->next_seg->next_seg->{prev_seg} = $cur_seg;
        $cur_seg->{next_seg} = $cur_seg->next_seg->next_seg;
    }
    else {
        $cur_seg->{next_seg} = undef;
    }

    splice(@{$self->segments}, $index+1, 1);
}

sub remove_cn_bkpts_without_cn_change {
    # Removes copy number (over-)segmentation boundaries, which don't result in
    # _significant_ copy number level changes. 
    #
    # NOTE: currently the method is implemented naively as a greedy merging
    # method where the segments are merged whenever possible from beginning to
    # end. This may produce inconsistent results with for example if the
    # segments were merged from the end to the beginning. 
    my $self = shift;
    my %params = @_;
    if (!defined($params{min_cn_change})) {
        die "Parameter 'min_cn_change' is required in $self\->remove_cn_bkpts_without_cn_change()";
    }
    if (!$self->is_sorted) {
        die "Attempted to call $self\->remove_cn_bkpts_without_cn_change() with an unsorted $self";
    }

    my($cur_seg, $total_len, $total_cn, $i);
    my $changed = 0;
    $i = 0;
    while ($i < $#{$self->segments}) {
        $cur_seg = $self->segments->[$i];
        if (
            # !$cur_seg->high_end->is_bounded_by_rg &&
            $cur_seg->high_end->is_bounded_by_cn_bkpt &&
            (
                !defined($cur_seg->cn) ||
                !defined($cur_seg->next_seg->cn) ||
                abs($cur_seg->cn - $cur_seg->next_seg->cn) < $params{min_cn_change}
            )
        ) {
            # Sanity check
            if ($cur_seg->next_seg->low_end->is_bounded_by_rg) {
                die;
            }

            # Merge the two segments together
            $self->merge_segment_with_next($i);
            $changed = 1;
        }
        else {
            $i++;
        }
    }

    return $changed;
}

sub normalise_single_bp_segments {
    my $self = shift;
    my %params = @_;
    if (!$self->is_sorted) {
        die "Attempted to call $self\->normalise_single_rg_clusters() with an unsorted $self";
    }

    my $seg = $self->segments->[0];
    my $fixed_segments = 0;
    while (defined($seg->next_seg)) {
        if (
                $seg->low_end->pos == $seg->high_end->pos  and
                !defined($seg->low_end->bkpt)  and
                !defined($seg->high_end->bkpt)  and
                defined($seg->prev_seg->high_end->bkpt)  and
                defined($seg->next_seg->low_end->bkpt)
        ) {
            $seg->low_end->{bkpt}  = $seg->next_seg->low_end->bkpt;
            $seg->low_end->{bkpt}->{segment_end} = $seg->low_end;
            $seg->low_end->{boundary} = "rg";
            $seg->next_seg->low_end->{bkpt} = undef;
            $seg->next_seg->low_end->{boundary} = "rg";
            $seg->high_end->{bkpt} = $seg->prev_seg->high_end->bkpt;
            $seg->high_end->{bkpt}->{segment_end} = $seg->high_end;
            $seg->high_end->{boundary} = "rg";
            $seg->prev_seg->high_end->{bkpt} = undef;
            $seg->prev_seg->high_end->{boundary} = "rg";

            $fixed_segments++;
        }
        else {
            $seg = $seg->next_seg;
        }
    }

    return $fixed_segments;
}

sub print_rg_cns_bedpe {
    my $self = shift;
    for ($self->segments_array) {
        print join(
            "\t",
            $self->name,
            $_->low_end->pos,
            $_->high_end->pos,
            ($_->cn || "NA"),
            $_->low_end->boundary,
            $_->high_end->boundary,
            ($_->n_win || "NA"),
        ) . "\n";
    }
}

# sub print_unique_breakpoints {
sub get_unique_breakpoints {
    my $self = shift;
    my %params = @_;
    if (!$self->is_sorted) {
        die "Attempted to call $self\->print_unique_breakpoints() with an unsorted $self";
    }

    my @out_bkpts = ();

    my $seg = $self->segments->[0];
    while (defined($seg->next_seg)) {
        if (defined($seg->high_end->bkpt)) {
#             ## Ignore small TDs and deletions
#             if (!$seg->high_end->bkpt->rg->is_long_range) {
#                 $seg = $seg->next_seg;
#                 next;
#             }

            if ($seg->high_end->bkpt->is_foldback(%params)) {
                if ($seg->high_end->bkpt->end eq "high") {
                    die;
                }

                # print join(
                push @out_bkpts, join(
                    "\t",
                    $self->name,
                    $seg->next_seg->centre,
                    ($seg->high_end->bkpt->id . "+," . $seg->high_end->bkpt->id . "+")
                ) . "\n";

                if ($seg->high_end->bkpt->mate->segment->is_bal_rg_overlap(%params)) {
                    if (
                        defined($seg->high_end->bkpt->mate->segment->low_end->bkpt) &&
                        $seg->high_end->bkpt->mate->segment->low_end->bkpt->is_foldback(%params)
                    ) {
                        # print join(
                        push @out_bkpts, join(
                            "\t",
                            $self->name,
                            $seg->next_seg->next_seg->next_seg->centre,
                            ($seg->high_end->mate->segment->low_end->bkpt->id . "-," . $seg->high_end->mate->segment->low_end->bkpt->id . "-")
                        ) . "\n";
                        $seg = $seg->next_seg->next_seg->next_seg->next_seg;
                    }
                    elsif (defined($seg->high_end->bkpt->mate->segment->low_end->bkpt)) {
                        # print join(
                        push @out_bkpts, join(
                            "\t",
                            $self->name,
                            $seg->next_seg->next_seg->high_end->pos,
                            $seg->next_seg->next_seg->high_end->bkpt->id . "+" 
                        ) . "\n";
                        $seg = $seg->next_seg->next_seg->next_seg;
                    }
                    else {
                        die;
                    }
                }
                else {
                    $seg = $seg->next_seg->next_seg;
                }
            }
            elsif ($seg->high_end->bkpt->has_balanced_bkpt_partner_rg_end(%params)) {
                if ($seg->is_bal_rg_overlap(%params)) {
                    if ($seg->low_end->bkpt->is_foldback(%params)) {
                        # print join(
                        push @out_bkpts, join(
                            "\t",
                            $self->name,
                            $seg->low_end->pos,
                            $seg->high_end->bkpt->id . "+"
                        ) . "\n";
                        # print join(
                        push @out_bkpts, join(
                            "\t",
                            $self->name,
                            $seg->next_seg->centre,
                            ($seg->low_end->bkpt->id . "-," . $seg->low_end->bkpt->id . "-")
                        ) . "\n";
                        $seg = $seg->next_seg->next_seg;
                    }
                    else {
                        # print join(
                        push @out_bkpts, join(
                            "\t",
                            $self->name,
                            $seg->centre,
                            ($seg->high_end->bkpt->id . "+," . $seg->low_end->bkpt->id . "-")
                        ) . "\n";
                        $seg = $seg->next_seg;
                    }
                }
                elsif ($seg->next_seg->is_bal_rg_gap(%params)) {
                    # print join(
                    push @out_bkpts, join(
                        "\t",
                        $self->name,
                        $seg->next_seg->centre,
                        ($seg->high_end->bkpt->id . "+," . $seg->next_seg->next_seg->low_end->bkpt->id . "-")
                    ) . "\n";
                    $seg = $seg->next_seg;
                }
                else {
                    die;
                }
            }
            else {
                # print join(
                push @out_bkpts, join(
                    "\t",
                    $self->name,
                    $seg->high_end->pos,
                    ($seg->high_end->bkpt->id . "+")
                ) . "\n";
                $seg = $seg->next_seg;
            }
        }
        elsif (defined($seg->next_seg->low_end->bkpt)) {
#             ## Ignore small TDs and deletions
#             if (!$seg->next_seg->low_end->bkpt->rg->is_long_range) {
#                 $seg = $seg->next_seg;
#                 next;
#             }

            if ($seg->next_seg->is_bal_rg_overlap(%params)) {
                # Below will be dealt in the next iteration of the above section
                $seg = $seg->next_seg;
                next;
            }
            elsif ($seg->next_seg->low_end->bkpt->is_foldback(%params)) {
                # print join(
                push @out_bkpts, join(
                    "\t",
                    $self->name,
                    $seg->next_seg->centre,
                    ($seg->next_seg->low_end->bkpt->id . "-," . $seg->next_seg->next_seg->low_end->bkpt->id . "-")
                ) . "\n";
                $seg = $seg->next_seg->next_seg;
            }
            else {
                # print join(
                push @out_bkpts, join(
                    "\t",
                    $self->name,
                    $seg->next_seg->low_end->pos,
                    ($seg->next_seg->low_end->bkpt->id . "-")
                ) . "\n";
                $seg = $seg->next_seg;
            }
        }
        else {
            $seg = $seg->next_seg;
        }
    }

    return @out_bkpts;
}

sub histogram_quantile {
    # First argument is quantile, remaining arguments are and array of
    # hashes with keys 'length' and 'value'.

    my $quantile = shift;
    my @data = sort { $a->{value} <=> $b->{value} } grep({ defined($_->{value}) } @_);
    if (@data == 0) {
        return undef;
    }
    my $total_length = List::Util::sum(map { $_->{length} } @data);
    my $cur_len = 0;
    while ($cur_len + $data[0]->{length} < $total_length * $quantile) {
        $cur_len += $data[0]->{length};
        shift @data;
    }
    return $data[0]->{value};
}

sub p_arm_telomere_cn {
    my $self = shift;
    if (!$self->is_sorted) {
        die "Attempted to call $self\->p_arm_telomere_cn() with an unsorted $self";
    }

    if ($self->segments->[0]->low_end->is_bounded_by_tel) {
        # return $self->segments->[0]->cn;

        my $cur_seg = $self->segments->[0];
        my $cur_len = 0;
        my @data;
        while ($cur_len < 1e6) {
            if ($cur_len + $cur_seg->length > 1e6) {
                push @data, { value => $cur_seg->cn, length => 1e6 - $cur_len };
                last;
            }
            else {
                push @data, { value => $cur_seg->cn, length => $cur_seg->length };
                $cur_len += $cur_seg->length;
                $cur_seg = $cur_seg->next_seg;
            }
        }

        return histogram_quantile(.5, @data);
    }
    else {
        die;
    }
}

sub q_arm_telomere_cn {
    my $self = shift;
    if (!$self->is_sorted) {
        die "Attempted to call $self\->q_arm_telomere_cn() with an unsorted $self";
    }

    if ($self->segments->[$#{$self->segments}]->high_end->is_bounded_by_tel) {
        # return $self->segments->[$#{$self->segments}]->cn;

        my $cur_seg = $self->segments->[$#{$self->segments}];
        my $cur_len = 0;
        my @data;
        while ($cur_len < 1e6) {
            if ($cur_len + $cur_seg->length > 1e6) {
                push @data, { value => $cur_seg->cn, length => 1e6 - $cur_len };
                last;
            }
            else {
                push @data, { value => $cur_seg->cn, length => $cur_seg->length };
                $cur_len += $cur_seg->length;
                $cur_seg = $cur_seg->prev_seg;
            }
        }

        return histogram_quantile(.5, @data);
    }
    else {
        die;
    }
}

sub p_centromere_cn {
    my $self = shift;
    if (!$self->is_sorted) {
        die "Attempted to call $self\->p_centromere_cn() with an unsorted $self";
    }

    my $cur_seg = $self->segments->[0];
    while (!$cur_seg->high_end->is_bounded_by_cen) {
        if (!defined($cur_seg) or !defined($cur_seg->next_seg)) {
            die;
        }
        $cur_seg = $cur_seg->next_seg;
    }
    # return $cur_seg->cn;

    my $cur_len = 0;
    my @data;
    while ($cur_len < 1e6) {
        if ($cur_len + $cur_seg->length > 1e6) {
            push @data, { value => $cur_seg->cn, length => 1e6 - $cur_len };
            last;
        }
        else {
            push @data, { value => $cur_seg->cn, length => $cur_seg->length };
            $cur_len += $cur_seg->length;
            $cur_seg = $cur_seg->prev_seg;
        }
    }

    return histogram_quantile(.5, @data);
}

sub q_centromere_cn {
    my $self = shift;
    if (!$self->is_sorted) {
        die "Attempted to call $self\->q_centromere_cn() with an unsorted $self";
    }

    my $cur_seg = $self->segments->[$#{$self->segments}];
    while (!$cur_seg->low_end->is_bounded_by_cen) {
        if (!defined($cur_seg) or !defined($cur_seg->prev_seg)) {
            die;
        }
        $cur_seg = $cur_seg->prev_seg;
    }
    # return $cur_seg->cn;

    my $cur_len = 0;
    my @data;
    while ($cur_len < 1e6) {
        if ($cur_len + $cur_seg->length > 1e6) {
            push @data, { value => $cur_seg->cn, length => 1e6 - $cur_len };
            last;
        }
        else {
            push @data, { value => $cur_seg->cn, length => $cur_seg->length };
            $cur_len += $cur_seg->length;
            $cur_seg = $cur_seg->next_seg;
        }
    }

    return histogram_quantile(.5, @data);
}

sub has_stable_arm_cn {
    # First argument is arm (p or q).
    # Requirements for stability are:
    # Telomere and centromere CN within 0.5.
    # 0.2 and 0.8 quantiles are within 0.5 of median CN. 
    my $self = shift;
    my $arm = shift;

    if (!$self->is_sorted) {
        die "Attempted to call $self\->has_stable_arm_cn() with an unsorted $self";
    }

    my $arm_cn;
    if ($arm eq 'p') {
        if (
            !defined($self->p_centromere_cn) ||
            !defined($self->p_arm_telomere_cn) ||
            abs($self->p_centromere_cn - $self->p_arm_telomere_cn) > 0.5
        ) {
            return(0);
        }
        $arm_cn = ($self->p_centromere_cn + $self->p_arm_telomere_cn)/2;
    }
    elsif ($arm eq 'q') {
        if (
            !defined($self->q_centromere_cn) ||
            !defined($self->q_arm_telomere_cn) ||
            abs($self->q_centromere_cn - $self->q_arm_telomere_cn) > 0.5
        ) {
            return(0);
        }
        $arm_cn = ($self->q_centromere_cn + $self->q_arm_telomere_cn)/2;
    }
    else {
        die "In $self\->has_stable_arm_cn('\$arm'), \$arm has to be either 'p' or 'q'.";
    }

    my ($low_quantile, $high_quantile, @data, $seg);
    @data = ();
    if ($arm eq 'p') {
        $seg = ($self->segments_array)[0];
        if (defined($seg->cn)) {
            push @data, { length => $seg->length, value => $seg->cn };
        }
        while (!$seg->high_end->is_bounded_by_cen) {
            if (!defined($seg->next_seg)) {
                die;
            }
            $seg = $seg->next_seg;
            if (defined($seg->cn)) {
                push @data, { length => $seg->length, value => $seg->cn };
            }
        }
    }
    else {
        $seg = ($self->segments_array)[-1];
        if (defined($seg->cn)) {
            push @data, { length => $seg->length, value => $seg->cn };
        }
        while (!$seg->low_end->is_bounded_by_tel) {
            if (!defined($seg->prev_seg)) {
                die;
            }
            $seg = $seg->prev_seg;
            if (defined($seg->cn)) {
                push @data, { length => $seg->length, value => $seg->cn };
            }
        }
    }

    $low_quantile  = histogram_quantile(0.2, @data);
    $high_quantile = histogram_quantile(0.8, @data);
    if (
        abs($low_quantile - $arm_cn)  <= 0.5 &&
        abs($high_quantile - $arm_cn) <= 0.5
    ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub mean_cn {
    # Chromosomal CN
    my $self = shift;
    my($sum, $count) = (0, 0);
    for (
        $self->p_arm_telomere_cn,
        $self->p_centromere_cn,
        $self->q_arm_telomere_cn,
        $self->q_centromere_cn
    ) {
        if (defined($_)) {
            $sum += $_;
            $count++;
        }
    }

    die if $count == 0;
    return $sum / $count;
}

# sub remove_balanced_rg_segments {
#     # Removes the small copy number segments between balanced rearrangement
#     # breakpoints.
#     my $self = shift;
#     if (!$self->is_sorted) {
#         die "Attempted to call $self\->remove_balanced_rg_segments() with an unsorted $self";
#     }
#     my %params = @_;
#     if (!exists($params{max_balanced_rg_dist})) {
#         die "Parameter 'max_balanced_rg_dist' is required in $self\->remove_balanced_rg_segments()";
#     }
# 
#     my $i = 0;
#     while ($i <= ($#{$self->segments} - 2)) {
#         if (
#             defined($self->segments->[$i]->high_end->bkpt) &&                    # First segment has a + rearrangement end
#             $self->segments->[$i+1]->length <= $params{max_balanced_rg_dist} &&  # Middle segment is smaller than max_balanced_rg_dist
#             !defined($self->segments->[$i+1]->low_end->bkpt) &&                  # ... and has no associated low end
#             !defined($self->segments->[$i+1]->high_end->bkpt) &&                 # ... or high end rearrangement
#             defined($self->segments->[$i+2]->low_end->bkpt) &&                   # Last segment has a - rearrangement end
#             $self->segments->[$i]->high_end->bkpt->id ne $self->segments->[$i+2]->low_end->bkpt->id  # And it's not a simple deletion
#         ) {
#             $self->segments->[$i]->{next_seg}   = $self->segments[$i+2];
#             $self->segments->[$i+2]->{prev_seg} = $self->segments[$i];
#             splice(@{$self->segments}, $i+1, 1);
#         }
#         else {
#             $i++;
#         }
#     }
# }
# 
# sub remove_foldback_rg_segments {
#     # Removes the small copy number segments between foldback rearrangement
#     # breakpoints.
#     my $self = shift;
#     if (!$self->is_sorted) {
#         die "Attempted to call $self\->remove_foldback_rg_segments() with an unsorted $self";
#     }
#     my %params = @_;
#     if (!exists($params{max_foldback_distance})) {
#         die "Parameter 'max_foldback_distance' is required in $self\->remove_foldback_rg_segments()";
#     }
# 
#     my $i = 0;
#     my $cur_seg;
#     while ($i <= ($#{$self->segments} - 2)) {
#         $cur_seg = $self->segments->[$i];
#             
#         if (
#             defined($cur_seg->low_end->bkpt) &&
#             !defined($cur_seg->high_end->bkpt) &&
#             $cur_seg->length <= $params{max_foldback_distance} &&
#             defined($cur_seg->next_seg->low_end->bkpt) &&
#             $cur_seg->low_end->bkpt->mate == $cur_seg->next_seg->low_end->bkpt
#         ) {
#             # TODO???
#         }
#     }
# }
#
# End of worker subroutines
#


1;
