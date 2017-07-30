use strict;
use warnings FATAL => 'all';

use lib '/nfs/users/nfs_y/yl3/programs/scripts/pancan_related_scripts/';
use CopyNumberSegmentEnd;

package Shard;

sub new {
    my $class = shift;
    my %params = @_;
    
    if (!exists($params{low_end_segment}) || !defined($params{low_end_segment})) {
        die "Parameter \$params{low_end_segment} is required in $self\->new()";
    }
    if (!exists($params{high_end_segment}) || !defined($params{high_end_segment})) {
        die "Parameter \$params{high_end_segment} is required in $self\->new()";
    }

    ## A few sanity checks. This is why we wanted segments rather than
    ## segment ends. 
    my $les = $params{low_end_segment};
    my $hes = $params{high_end_segment};
    if (
        $les->chr != $hes->chr ||
        !$les->chr->is_sorted ||
        $les->low_end->pos >= $hes->high_end->pos ||
        !defined($les->low_end->bkpt) ||
        !defined($hes->high_end->bkpt)
    ) {
        die;
    }

    return bless {
        low_end => $params{low_end_segment}->low_end,
        high_end => $params{high_end_segment}->high_end,
    }, $class;

    return $self;
}

#
# Getter methods
#
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
#
# End of getter methods
#


#
# Helper methods
#
sub print {
    my $self = shift;
    print $self->low_end->bkpt->rg->to_s . "\n";
    my $cur_seg = $self->low_end->segment->prev_seg;
    while ($cur_seg != $self->high_end->segment) {
        $cur_seg = $cur_seg->next_seg;
        $cur_seg->print;
    }
    print $self->high_end->bkpt->rg->to_s;
}
#
# End of helper methods
#


1;
