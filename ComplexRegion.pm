use strict;
use warnings FATAL => 'all';

use Scalar::Util qw(blessed);

package ComplexRegion;
#
# A linear genomic region of high rearrangement breakpoint density
#

sub new {
    my $class = shift;
    my %params = @_;
    my $region = bless {
        chr_name => $params{chr_name},
        start => $params{start},
        end => $params{end},
    }, $class;

    my($rg_end, @rg_ends);
    my %rgs;
    my %components;
    for my $seg (@{$params{segments}}) {
        $rg_end = $seg->low_end->bkpt;
        if (
            defined($rg_end) &&
            $region->overlaps_rg_end($rg_end)
        ) {
            push @rg_ends, $rg_end;
            $rgs{$rg_end->id} = $rg_end->rg;
            $components{$rg_end->rg->component} = $rg_end->rg->component;
        }

        $rg_end = $seg->high_end->bkpt;
        if (
            defined($rg_end) &&
            $region->overlaps_rg_end($rg_end)
        ) {
            push @rg_ends, $rg_end;
            $rgs{$rg_end->id} = $rg_end->rg;
            $components{$rg_end->rg->component} = $rg_end->rg->component;
        }
    }
    $region->{rg_ends} = \@rg_ends;
    $region->{rgs} = \%rgs;
    $region->{components} = \%components;

    return $region;
}

sub chr_name {
    my $self = shift;
    return $self->{chr_name};
}

sub start {
    my $self = shift;
    return $self->{start};
}

sub end {
    my $self = shift;
    return $self->{end};
}

sub to_s {
    my $self = shift;
    return $self->chr_name . ":" . $self->start . "-" . $self->end;
}

sub rgs {
    my $self = shift;
    return $self->{rgs};
}

sub components {
    my $self = shift;
    return $self->{components};
}

sub overlaps_rg_end {
    my $self = shift;
    my $rg_end = shift;
    if (!defined(Scalar::Util::blessed($rg_end)) || Scalar::Util::blessed($rg_end) ne "RearrangementEnd") {
        die "$self\->overlaps_rg_end() requires an argument of class RearrangementEnd";
    }

    return(
        $rg_end->chr_name eq $self->chr_name &&
        $rg_end->pos >= $self->start &&
        $rg_end->pos <= $self->end
    );
}

sub long_range_rg_count {
    my $self = shift;
    my $long_range_rg_count = 0;
    for (values %{$self->rgs}) {
        if ($_->is_long_range) {
            $long_range_rg_count++;
        }   
    }   
    return $long_range_rg_count;
}

sub is_complex {
    my $self = shift;
    if ($self->long_range_rg_count >= 3) {  ## TODO: put as argument
        return 1;
    }   
    else {
        return 0;
    }   
}

sub is_connected_to {
    ## Is there a rearrangement connecting this ComplexRegion to the target
    ## ComplexRegion?
    my $self = shift;
    my $target = shift;

    if (Scalar::Util::blessed($target) eq "ComplexRegion") {
        for my $rg_id (keys %{$self->rgs}) {
            if (exists($target->rgs->{$rg_id})) {
                return 1;
            }
        }
        return 0;
    }
    elsif (Scalar::Util::blessed($target) eq "ComplexRegionList") {
        for my $region ($target->regions_array) {
            for my $rg_id (keys %{$self->rgs}) {
                if (exists($region->rgs->{$rg_id})) {
                    return 1;
                }
            }
        }
        return 0;
    }
    else {
        die;
    }
}

1;
