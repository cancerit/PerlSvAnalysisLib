use strict;
use warnings FATAL => 'all';

use File::Basename;
use lib dirname($0);
use CopyNumberSegmentEnd;

package CopyNumberSegment;

sub new {
    # Required inputs:
    #   rg_of_id
    #   chr
    #   start
    #   end
    #   low_end_bkpt
    #   high_end_bkpt
    #   cn
    #   n_win

    my $class = shift;
    my %params = @_;
    my %rg_of_id = %{$params{rg_of_id}};  # A hash mapping rearrangement IDs to rearrangements
    my $rg_end;

    my $low_end = CopyNumberSegmentEnd->new(
        pos => $params{start},
        end => "low",
    );  
    if ($params{low_end_bkpt} =~ /^(.+):([12])$/) {
        if (!exists($rg_of_id{$1})) {
            die "Rearrangement $1 found in RG_CNS file not found in the RG_BEDPE file!";
        }
        $rg_end = ($2 == 1 ? "low_end" : "high_end");
        $low_end->{boundary} = "rg";
        if ($rg_of_id{$1}->{$rg_end}->dir eq "-") {
            $low_end->{bkpt} = $rg_of_id{$1}->{$rg_end};
            $rg_of_id{$1}->{$rg_end}->{segment_end} = $low_end;
        }
    }   
    else {
        $low_end->{boundary} = $params{low_end_bkpt};
    }   

    my $high_end = CopyNumberSegmentEnd->new(
        pos => $params{end},
        end => "high",
    );  
    if ($params{high_end_bkpt} =~ /^(.+):([12])$/) {
        if (!exists($rg_of_id{$1})) {
            die "Rearrangement $1 found in RG_CNS file not found in the RG_BEDPE file!";
        }
        $rg_end = ($2 == 1 ? "low_end" : "high_end");
        $high_end->{boundary} = "rg";
        if ($rg_of_id{$1}->{$rg_end}->dir eq "+") {
            $high_end->{bkpt} = $rg_of_id{$1}->{$rg_end};
            $rg_of_id{$1}->{$rg_end}->{segment_end} = $high_end;
        }
    }   
    else {
        $high_end->{boundary} = $params{high_end_bkpt};
    }

    $low_end->{other_end} = $high_end;
    $high_end->{other_end} = $low_end;

    my $self = bless {
        chr      => ($params{chr} or undef),  ## Pointer to object of type 'CopyNumberSegmentArray'
        cn       => $params{cn},
        n_win    => $params{n_win},
        low_end  => $low_end,
        high_end => $high_end,
        prev_seg => ($params{prev_seg} or undef),
        next_seg => ($params{next_seg} or undef),
    }, $class;
    $self->{low_end}->{segment} = $self;
    $self->{high_end}->{segment} = $self;

    return $self;
}

#
# Getter methods
#
sub chr {
    my $self = shift;
    if (!defined($self->{chr})) {
        die "Attempted to call $self\->chr_name() to get undefined $self\->{chr}!";
    }
    return $self->{chr};
}

sub chr_name {
    my $self = shift;
    return $self->chr->name;
}

sub cn {
    my $self = shift;
    # if (!defined($self->{cn})) {
    #     die "Attempted to call cn() to get undefined $self\->{cn}!";
    # }
    return $self->{cn};
}

sub cn_s {
    my $self = shift;
    return($self->cn or "NA");
}

sub n_win {
    my $self = shift;
    if (!defined($self->{n_win})) {
        die "Attempted to call n_win() to get undefined $self\->{n_win}!";
    }
    return $self->{n_win};
}

sub low_end {
    my $self = shift;
    if (!defined($self->{low_end})) {
        die "Attempted to call low_end() to get undefined $self\->{low_end}!";
    }
    return $self->{low_end};
}

sub high_end {
    my $self = shift;
    if (!defined($self->{high_end})) {
        die "Attempted to call high_end() to get undefined $self\->{high_end}!";
    }
    return $self->{high_end};
}

sub next_seg {
    my $self = shift;
    # if (!defined($self->{next_seg})) {
    #     die "Attempted to call next_seg() to get undefined $self\->{next_seg} at segment " . $self->segment_position;
    # }
    return $self->{next_seg};
}

sub prev_seg {
    my $self = shift;
    # if (!defined($self->{prev_seg})) {
    #     die "Attempted to call prev_seg() to get undefined $self\->{prev_seg} at segment " . $self->segment_position;
    # }
    return $self->{prev_seg};
}
#
# End of getter methods
#


#
# Helper methods
#
sub centre {
    my $self = shift;
    return int($self->high_end->pos/2 + $self->low_end->pos/2);
}

sub length {
    my $self = shift;
    return($self->high_end->pos - $self->low_end->pos + 1);
}

sub idx_on_chr {
    my $self = shift;
    my $chr = $self->chr;
    if (!($chr->is_sorted)) {
        warn "Call $self\->idx_on_chr() on a CopyNumberSegment on an unsorted CopyNumberSegmentArray.\n";
    }
    my $i = 0;
    for (@{$chr->segments}) {
        if ($_ == $self) {
            return $i;
        }
        else {
            $i++;
        }
    }

    die;
}

sub get_matching_pattern {
    my $pattern = shift;
    my $text = shift;
    if ($text =~ /$pattern/) {
        return $&;
    }
    else {
        return undef;
    }
}

sub is_shard {
    my $self = shift;
    my %params = @_;
    if (!exists($params{max_shard_length})) {
        die "Attempted to call $self\->is_shard() with undefined parameter 'max_shard_length'";
    }
    return(
        $self->length <= $params{max_shard_length} &&
        defined($self->high_end->bkpt) &&
        defined($self->low_end->bkpt) &&
        $self->low_end->bkpt->id ne $self->high_end->bkpt->id
    );
}

sub is_bal_rg_gap {
    my $self = shift;
    my %params = @_;
    return(
        !defined($self->low_end->bkpt) &&
        !defined($self->high_end->bkpt) &&
        $self->length <= $params{max_balanced_rg_dist} &&
        defined($self->prev_seg) &&
        !$self->prev_seg->is_bal_rg_overlap(%params) &&
        defined($self->prev_seg->high_end->bkpt) &&
        defined($self->next_seg) &&
        !$self->next_seg->is_bal_rg_overlap(%params) &&
        defined($self->next_seg->low_end->bkpt) &&
        $self->prev_seg->high_end->bkpt->rg->id ne $self->next_seg->low_end->bkpt->rg->id
    );
}

sub is_bal_rg_overlap {
    my $self = shift;
    my %params = @_;
    if (!exists($params{max_balanced_rg_overlap})) {
        die "$self\->is_bal_rg_overlap called without parameter max_balanced_rg_overlap!";
    }
    return(
        defined($self->low_end->bkpt) &&
        defined($self->high_end->bkpt) &&
        $self->length <= $params{max_balanced_rg_overlap} &&
        $self->low_end->bkpt->min_reads_pos - $self->high_end->bkpt->min_reads_pos >= 5 &&
        $self->low_end->bkpt->max_reads_pos - $self->high_end->bkpt->max_reads_pos >= 5 &&
        $self->low_end->bkpt->rg->id ne $self->high_end->bkpt->rg->id
    );
}

sub is_on_p_arm {
    my $self = shift;
    my $seg = $self;
    while (!$seg->low_end->is_bounded_by_tel && !$seg->low_end->is_bounded_by_cen) {
        if (!defined($seg->prev_seg)) {
            die;
        }
        $seg = $seg->prev_seg;
    }

    return($seg->low_end->is_bounded_by_tel);
}

sub is_on_q_arm {
    my $self = shift;
    return(!$self->is_on_p_arm);
}

sub segment_position {
    my $self = shift;
    return $self->chr_name . ":" . $self->low_end->pos . "-" . $self->high_end->pos;
}

sub to_s {
    my $self = shift;
#     print "Position: " . $self->segment_position . "\n";
#     print "CN: " . $self->cn . "\n";
#     print "n_win: " . $self->n_win . "\n";
#     print "Low end boundary: " . $self->low_end->boundary . "\n";
#     print "Low end bkpt: " . $self->low_end->bkpt . "\n" if defined $self->low_end->bkpt;
#     print "High end boundary: " . $self->high_end->boundary . "\n";
#     print "High end bkpt: " . $self->high_end->bkpt . "\n" if defined $self->high_end->bkpt;
    return join(
        "\t",
        $self->segment_position,
        sprintf("%.1f", ($self->high_end->pos - $self->low_end->pos + 1)/1000),
        (defined($self->cn) ? sprintf("%.2f", $self->cn) : "NA"),
        $self->n_win,
        (
            defined($self->low_end->bkpt) ?
            get_matching_pattern('\d+$', $self->low_end->bkpt->id) . ":" . $self->low_end->bkpt->dir . ":" . substr($self->low_end->bkpt->end, 0, 1) :
            $self->low_end->boundary
        ),
        (
            defined($self->high_end->bkpt) ?
            get_matching_pattern('\d+$', $self->high_end->bkpt->id) . ":" . $self->high_end->bkpt->dir . ":" . substr($self->high_end->bkpt->end, 0, 1) :
            $self->high_end->boundary
        )
    );
}

sub print {
    my $self = shift;
    print($self->to_s, "\n");
}
#
# End of helper methods
#


1;
