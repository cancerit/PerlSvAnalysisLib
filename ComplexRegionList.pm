use strict;
use warnings FATAL => 'all';

use List::Util qw(min max);

use lib '/nfs/users/nfs_y/yl3/programs/scripts/pancan_related_scripts';
use ComplexRegion;

package ComplexRegionList;

sub new {
    my $class = shift;
    my %params = @_;
    return bless {
        regions => $params{regions},
    }, $class;
}

sub new_from_clusters_file {
    my $class = shift;
    my $file = shift;
    my $genome = shift;

    open IN, $file;
    my $cur_chr;
    my @regions_array;
    my @positions;
    while (<IN>) {
        chomp;
        if (/^>(.+)$/) {
            $cur_chr = $1;
        }
        else {
            @positions = split /,/;
            push @regions_array, ComplexRegion->new(
                chr_name => $cur_chr,
                start => List::Util::min(@positions),
                end => List::Util::max(@positions),
                segments => [$genome->chr_by_name($cur_chr)->segments_array],
            );
        }
    }

    return bless {
        regions => \@regions_array,
    }, $class;
}

sub regions_array {
    my $self = shift;
    return @{$self->{regions}};
}

sub is_complex {
    my $self = shift;
    return $self->{is_complex};
}

sub add_region {
    my $self = shift;
    my $target_region = shift;
    push @{$self->{regions}}, $target_region;
}

sub rgs {
    my $self = shift;
    my %rgs = ();
    for my $region ($self->regions_array) {
        %rgs = (%rgs, %{$region->rgs});
    }
    return \%rgs;
}

sub components_array {
    my $self = shift;
    my %components = ();
    for my $region ($self->regions_array) {
        %components = (%components, %{$region->components});
    }
    return(values %components);
}

sub merge_complex_regions_into_events {
    my $all_regions = shift;
    my @events;

    ## Algorithm has three rounds.
    ## 1: Find index events with more than 3 index rearrangements.
    ## 2: Merge index and other regions into events. 
    ## 3: Mark rearrangement ends who are not part of complex events,
    ##    but whose mates fall into complex event regions. 

    my($region, $event);

    COMPLEX_REGION:
    for $region (grep { $_->is_complex } ($all_regions->regions_array)) {
        for $event (@events) {
            if ($region->is_connected_to($event)) {
                $event->add_region($region);
                next COMPLEX_REGION;
            }
        }

        ## If got this far then we didn't join this region to any event
        push @events, ComplexRegionList->new(
            regions => [$region],
        );
    }

    SIMPLE_REGION:
    for $region (grep { !$_->is_complex } ($all_regions->regions_array)) {
        ## Either join an existing event or start a new event.
        for $event (@events) {
            if ($region->is_connected_to($event)) {
                $event->add_region($region);
                next SIMPLE_REGION;
            }
        }
    }

    return @events;
}



1;
