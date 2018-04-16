# FootprintCluster.pm
# A footprint cluster is a group footprints that are connected through
# SVs. Through transitive proximity, footprints and thus rearrangements
# as part of footprints, are assumed to have arisen in the same event
# and thus are annotated as being part of a footprint cluster. 

package FootprintCluster;

use strict;
use Footprint;
use List::Util qw(min max);
use warnings FATAL => 'all';

sub new {
    my $class = shift;
    my %params = @_;
    
    if (!exists($params{id})) {
        die "Parameter \$params{id} is needed in $class\->new(%params)";
    }

    if (!exists($params{genome})) {
        die "Parameter \$params{genome} is needed in $class\->new(%params)";
    }

    my $footprint_cluster = bless {
        id         => $params{id},
        footprints => {},
        genome     => $params{genome},
    }, $class;

    return $footprint_cluster;
}

sub destroy {
    my $self = shift;
    delete $self->{genome}->{rg_clusters}->{$self->id};
}

sub id {
    my $self = shift;
    return $self->{id};
}

sub footprints {
    my $self = shift;
    return $self->{footprints};
}

sub footprints_array {
    my $self = shift;
    return values($self->footprints);
}

sub rgs {
    my $self = shift;
    my %rgs = ();
    for my $fp ($self->footprints_array) {
        $rgs{$_} = $_ for $fp->rgs_array;
    }
    return \%rgs;
}

sub rgs_array {
    my $self = shift;
    return values($self->rgs);
}

sub n_rgs {
    my $self = shift;
    return scalar($self->rgs_array);
}

sub rg_ends_array {
    my $self = shift;
    my @rg_ends = ();
    for ($self->footprints_array) {
        push @rg_ends, $_->rg_ends_array;
    }
    return @rg_ends;
}

sub size {
    my $self = shift;
    return scalar($self->footprints_array);
}

sub add_footprint {
    # Input is footprints

    my $self = shift;
    for (@_) {
        if (!defined(Scalar::Util::blessed $_) or Scalar::Util::blessed($_) ne "Footprint") {
            die "In $self\->add_footprint(), trying to add an object of not type 'Footprint'";
        }
        if (exists($self->footprints->{$_->id})) {
            warn("Warning: in $self\->add_footprint(), the footprint cluster %s already contains footprint %s.\n", $self->id, $_->id);
        }
        elsif ($_->cluster != $self) {
            die sprintf("In $self\->add_footprint(), the footprint %s is already assigned to cluster %s, but is attempted to be assigned to %s", $_->id, $_->cluster->id, $self->id);
        }
        else {
            $self->footprints->{$_->id} = $_;
            $_->set_cluster($self);
        }
    }
}

sub remove_footprint {
    my $self = shift;
    for (@_) {
        if (!exists($self->footprints->{$_->id})) {
            die sprintf("In $self\->remove_footprint(), from cluster %s attempted to remove footprint %s, which is not a member of the cluster", $self->id, $_->id);
        }
        delete $self->footprints->{$_->id};
        $_->{cluster} = undef;
    }
}

sub only_has_deletions {
    # Whether a SV cluster only has deletion-type SVs
    my $self = shift;
    for my $rg ($self->rgs_array) {
        if ($rg->is_del_type) {
            return 0;
        }
    }

    return 1;
}

sub is_intra_chromosomal {
    my $self = shift;
    my %chrs = ();
    for my $f ($self->footprints_array) {
        $chrs{$f->chr_name} = 1;
    }
    return(scalar(keys(%chrs)) == 1);
}


sub cluster_classification {
    # Classification of a SV cluster aka footprint cluster
    #
    # Division of cluster types
    # * A single rearrangement
    # * Clusters with one footprint
    #   * Two rearrangements (returned label)
    #     * Direct inversion (direct_inversion)
    #     * Inversion deletion-gain (inversion_gain_loss)
    #     * Inverted duplication (inverted_duplication)
    #     * Inverted dup-trp-dup (dup_trp_dup)
    #     * Two deletions, two TDs, including two overlapping TDs (complex_unclear)
    #     * Fold-back + (deletion or TD) (complex_unclear)
    #     * Complex (complex_unclear)
    # * Chains and cycles of shards and balanced breakpoints
    #     * shard_(chain|cycle):<length>
    #     * chromoplexy_(chain|cycle):<length>
    #
    # The default type is "complex_unclear"

    my $self = shift;
    my %params = @_;
    my($fp, @rgs);
    my $cluster_type;


    sub interpret_cluster_of_size_1 {
        # Interpret a cluster with exactly one footprint
        my $self = shift;
        my %params = @_;

        if (scalar($self->rgs_array) == 1) {
            # Sanity check = this shouldn't happen here.
            die;
        }

        sub has_only_del_and_td_svs {
            my $self = shift;
            return (scalar(grep {$_->is_td_type or $_->is_del_type} $self->rgs_array)
                    == scalar($self->rgs_array));
        }

        my $footprint = ($self->footprints_array)[0];
        if ($self->n_rgs == 2) {
            # Two rearrangements in a single footprint
            return $footprint->type(%params);
        }
        else {
            # Three or more SVs within a single footprint
            return "complex_unclear";
        }  # End of three or more SVs in a single footprint cluster
    }  # End of sub interpret_cluster_of_size_1()

    # sub two_jump_structure_type {
    #     # This function categorizes the basic two-jump structures
    #     # that are made of two rather than one footprints.
    #     
    #     my $self = shift;
    #     my %params = @_;
    #     
    #     # Basic requirement checks
    #     return "" if $self->size != 2;
    #     my($footprint_1, $footprint_2) = $self->footprints_array;
    #     return "" if $footprint_1->chr_name ne $footprint_2->chr_name;
    #     return "" if $self->n_rgs != 2;
    #     return "" if grep({not $_->is_inv_type} $self->rgs_array);

    #     if ($footprint_2->last_pos < $footprint_1->first_pos) {
    #         ($footprint_1, $footprint_2) = ($footprint_2, $footprint_1);
    #     }
    #     elsif ($footprint_1->last_pos >= $footprint_2->first_pos) {
    #         die;  # Sanity check
    #     }

    #     # my $footprints_distance = min(map {$_->pos} $footprint_2->rg_ends_array)
    #     #          - max(map {$_->pos} $footprint_1->rg_ends_array);
    #     # return "" if $footprints_distance < $params{max_local_two_jump_distance};

    # 
    #     # Different versions of direct inversion
    #     my @footprint_1_rg_ends = $footprint_1->sorted_rg_ends_array(%params);
    #     my @footprint_2_rg_ends = $footprint_2->sorted_rg_ends_array(%params);
    #     if (
    #             $footprint_1->size == 2 and
    #             $footprint_1->n_rgs == 2 and
    #             $footprint_1->is_balanced_type(%params) and
    #             $footprint_2->is_balanced_type(%params)
    #     ) {
    #         return "direct_inversion";
    #     }
    #     if (
    #             $footprint_1->size == 1 and
    #             ($footprint_1_rg_ends[0])->is_fwd and
    #             ($footprint_1_rg_ends[0])->mate == $footprint_2_rg_ends[1] and
    #             ($footprint_2_rg_ends[0])->is_rev
    #     ) {
    #         return "direct_inversion";
    #     }
    #     if (
    #             $footprint_2->size == 1 and
    #             ($footprint_2_rg_ends[0])->is_rev and
    #             ($footprint_2_rg_ends[0])->mate == $footprint_1_rg_ends[1] and
    #             ($footprint_1_rg_ends[2])->is_fwd
    #     ) {
    #         return "direct_inversion";
    #     }

    #     # Inverted gain-losses. 
    #     if (
    #             $footprint_1->size == 1 and
    #             ($footprint_1_rg_ends[0])->is_fwd and
    #             ($footprint_1_rg_ends[0])->mate == $footprint_2_rg_ends[2] and
    #             ($footprint_2_rg_ends[0])->is_rev
    #     ) {
    #         return "inversion_gain_loss";
    #     }
    #     if (
    #             $footprint_2->size == 1 and
    #             ($footprint_2_rg_ends[0])->is_rev and
    #             ($footprint_2_rg_ends[0])->mate == $footprint_1_rg_ends[0] and
    #             ($footprint_1_rg_ends[2])->is_fwd
    #     ) {
    #         return "inversion_gain_loss";
    #     }

    #     # Inverted duplication
    #     if (
    #             $footprint_1->size == 2 and
    #             $footprint_1->n_rgs == 2 and
    #             $footprint_1->is_shard_type(%params) and
    #             $footprint_2->is_shard_type(%params)
    #     ) {
    #         return "inverted_duplication";
    #     }
    #     if (
    #             $footprint_1->size == 1 and
    #             ($footprint_1_rg_ends[0])->is_rev and
    #             ($footprint_1_rg_ends[0])->mate == $footprint_2_rg_ends[1] and
    #             ($footprint_2_rg_ends[0])->is_fwd
    #     ) {
    #         return "inverted_duplication";
    #     }
    #     if (
    #             $footprint_2->size == 1 and
    #             ($footprint_2_rg_ends[0])->is_fwd and
    #             ($footprint_2_rg_ends[0])->mate == $footprint_1_rg_ends[1] and
    #             ($footprint_1_rg_ends[2])->is_rev
    #     ) {
    #         return "inverted_duplication";
    #     }

    #     # Dup-trp-dup
    #     if (
    #             $footprint_1->size == 2 and 
    #             $footprint_1->n_rgs == 2 and
    #             scalar(grep {$_->is_rev} $footprint_1->rg_ends_array) == 2 and
    #             scalar(grep {$_->is_fwd} $footprint_2->rg_ends_array) == 2
    #     ) {
    #         return "dup_trp_dup";
    #     }
    #     if (
    #             $footprint_1->size == 1 and
    #             ($footprint_1_rg_ends[0])->is_rev and
    #             ($footprint_1_rg_ends[0])->mate == $footprint_2_rg_ends[0] and
    #             ($footprint_2_rg_ends[2])->is_fwd
    #     ) {
    #         return "dup_trp_dup";
    #     }
    #     if (
    #             $footprint_2->size == 1 and
    #             ($footprint_2_rg_ends[0])->is_fwd and
    #             ($footprint_2_rg_ends[0])->mate == $footprint_1_rg_ends[2] and
    #             ($footprint_1_rg_ends[0])->is_rev
    #     ) {
    #         return "dup_trp_dup";
    #     }

    #     # Default:
    #     return "";
    # }

    if (scalar($self->rgs_array) == 1) {
        # Cluster with a single rearrangement
        return (($self->rgs_array)[0])->rg_type_s(%params);
    }
    elsif (
            $self->has_only_del_and_td_svs()  and
            $self->size == 1  # To avoid calling two-jumps as simple dels and TDs
    ) {
        # All rearrangements are deletions or tandem duplications
        return "simple_dels_and_tds";
    }
    elsif ($self->size == 1) {
        # Single cluster, but more than one rearrangement
        return $self->interpret_cluster_of_size_1(%params);
    }
    # elsif ($self->two_jump_structure_type(%params)) {
    #     ### Removed to not overcall two-jumps ###
    #     # Two-jump with two footprints
    #     return $self->two_jump_structure_type(%params);
    # }
    elsif (
            $self->size == 2 and
            scalar(grep {$_->is_balanced_type(%params)} $self->footprints_array) == 2 and
            $self->is_intra_chromosomal and
            scalar(grep {$_->is_inv_type} $self->rgs_array) == 2
    ) {
        return "direct_inversion";
    }
    else {
        # Multi-footprint clusters that aren't two-jumps
        #
        my @footprints = $self->footprints_array;

        my @shard_footprints = grep { $_->is_shard_type(%params) } @footprints;
        my @bal_footprints = grep { $_->is_balanced_type(%params) } @footprints;
        if (@shard_footprints) {
            # There are shards - test if the entire event is just a shard
            # chain or cluster.
            my $shard_type = ($shard_footprints[0])->detailed_type(%params);
            if ($shard_type =~ /^shard_(cycle|chain):(\d+)/) {
                # If the event is a pure shard chain, then the number of shards
                # must be the same as the number of footprints + 2. Moreover
                # the end footprints must be singleton footprints. 
                if ($1 eq "cycle" and $2 == $self->size) {
                    # For a shard cycle of two rearrangements, this can also be
                    # unbalanced translocation + TD. 
                    if ($2 == 2) {
                        my @rgs = $shard_footprints[0]->rgs_array;
                        my @cn_changes = (
                            $rgs[0]->weighted_avg_cn_change_across_rg(%params),
                            $rgs[1]->weighted_avg_cn_change_across_rg(%params)
                        );
                        if (grep {!defined($_)} @cn_changes) {
                            return $shard_type
                        }
                        else {
                            @cn_changes = sort(@cn_changes);
                            my $MIN_SEG_SIZE = 200_000;
                            my $bkpt_0_l = ($shard_footprints[0]->sorted_rg_ends_array(%params))[0];
                            my $bkpt_0_h = ($shard_footprints[0]->sorted_rg_ends_array(%params))[1];
                            my $bkpt_1_l = ($shard_footprints[1]->sorted_rg_ends_array(%params))[0];
                            my $bkpt_1_h = ($shard_footprints[1]->sorted_rg_ends_array(%params))[1];

                            if ($cn_changes[0] <= 0) {
                                return $shard_type;
                            }
                            elsif (
                                # $cn_changes[1] - $cn_changes[0] > 0.5 and
                                abs(log($cn_changes[1]/$cn_changes[0])) > log(1.5) and
                                $bkpt_0_l->segment->prev_seg->length >= $MIN_SEG_SIZE and
                                $bkpt_0_h->segment->next_seg->length >= $MIN_SEG_SIZE and
                                $bkpt_1_l->segment->prev_seg->length >= $MIN_SEG_SIZE and
                                $bkpt_1_h->segment->next_seg->length >= $MIN_SEG_SIZE
                            ) {
                                return "TD_after_unbal_transloc";
                            }
                            else {
                                return $shard_type;
                            }
                        }
                    }
                    else {
                        return $shard_type;
                    }
                }  # End of shard cycles
                elsif (
                        $1 eq "chain"
                        and $2 == $self->size - 2
                        and scalar(grep { !($_->is_shard_type(%params)) and $_->size == 1 } @footprints) == 2
                ) {
                    # Two single-breakpoint footprints connected by shards in a chain
                    return $shard_type;
                }
                elsif (
                        $1 eq "chain"
                        and $2 == $self->size - 1
                        and scalar(@bal_footprints) == 1
                ) {
                    # Insertion of a shard (chain)
                    return "shard_chain_insertion";
                }
                elsif (
                        scalar(@bal_footprints) >= 1
                        and scalar(@shard_footprints) >= 1
                ) {
                    # This cluster has both balanced breakpoint and shard footprints. 
                    # Test whether it is a pure balanced breakpoint + shard
                    # cycle or chain. 
                    my $footprint_count = 1;
                    my ($bkpt_1, $bkpt_2) = $bal_footprints[0]->sorted_rg_ends_array(%params);
                    my $starting_bkpt = $bkpt_1;

                    while (
                            $bkpt_1->mate->footprint->is_balanced_type(%params) or
                            $bkpt_1->mate->footprint->is_shard_type(%params)
                    ) {
                        if ($bkpt_1->mate->footprint->is_balanced_type(%params)) {
                            $bkpt_1 = $bkpt_1->mate->balanced_bkpt_partner_rg_end_by_footprint(%params);
                        }
                        else {
                            $bkpt_1 = $bkpt_1->mate->shard_partner_rg_end_by_footprint(%params);
                        }

                        if ($bkpt_1 == $starting_bkpt) {
                            if ($footprint_count == $self->size) {
                                return sprintf(
                                    "chromoplexy_%s_plus_templated_insertions:%s,%s",
                                    (scalar(@bal_footprints) == 1 ? "chain" : "cycle"),  # Single balanced breakpoint cannot be a "cycle"
                                    scalar(@bal_footprints),
                                    scalar(@shard_footprints)
                                );
                            }
                            else {
                                return "complex_unclear";
                            }
                        }    
                        $footprint_count++;
                    }

                    while (
                            $bkpt_2->mate->footprint->is_balanced_type(%params) or
                            $bkpt_2->mate->footprint->is_shard_type(%params)
                    ) {
                        if ($bkpt_2->mate->footprint->is_balanced_type(%params)) {
                            $bkpt_2 = $bkpt_2->mate->balanced_bkpt_partner_rg_end_by_footprint(%params);
                        }
                        else {
                            $bkpt_2 = $bkpt_2->mate->shard_partner_rg_end_by_footprint(%params);
                        }
                        $footprint_count++;
                    }    

                    # If this is a pure chain with only balanced and templated
                    # insertion footprints, then the following must be true.
                    if (
                            $footprint_count + 2 == $self->size  and
                            $footprint_count == scalar(@bal_footprints) + scalar(@shard_footprints)  and
                            $self->n_rgs == $self->size + 1
                    ) {
                        return sprintf(
                            "chromoplexy_chain_plus_templated_insertions:%s,%s",
                            scalar(@bal_footprints),
                            scalar(@shard_footprints)
                        );
                    }

                    # Fallback case, if current footprint cluster is not composed
                    # of just a chain or cycle of balanced and templated
                    # insertion footprints. 
                    return "complex_unclear";
                }
                else {
                    return "complex_unclear";
                }
            }  # End of clusters with shard cycles and/or shard chains
        }  # End of clusters with shard footprints
        elsif (@bal_footprints) {
            # There are balanced breakpoints but no shards - test if the entire
            # event is just balanced breakpoint chain or cluster.
            my $balanced_type = ($bal_footprints[0])->detailed_type(%params);
            if ($balanced_type =~ /chromoplexy_(cycle|chain):(\d+)/) {
                # If the event is a pure balanced breakpoint chain, then the
                # number of balanced breakpoints must be the same as the number
                # of footprints + 2. Moreover the end footprints must be
                # singleton footprints.
                if ($1 eq "cycle" and $2 == $self->size) {
                    return $balanced_type;
                }
                elsif (
                    $1 eq "chain"
                    and $2 == $self->size - 2
                    and scalar(grep { !($_->is_balanced_type(%params)) and $_->size == 1 } @footprints) == 2
                ) {
                    # Two single-breakpoint footprints connected by shards in a chain. 
                    if ($2 == 1) {  # chromoplexy_chain:1
                        return "split_translocation";
                    }
                    else {  # chromoplexy_chain:x, where x > 1
                        return $balanced_type;
                    }
                }
                else {
                    return "complex_unclear";
                }
            }
        }  # End of clusters with balanced breakpoint breakpoint footprints
        elsif (
                $self->size == 2 and
                $self->n_rgs == 2 and
                grep {$_->type(%params) ne "complex_unclear"} @footprints
        ) {
            # Otherwise if there are two SVs and two footprints, one of which is
            # a known type, then that's the even type too. 
            return join(",", grep {$_ ne "single"} map({$_->type(%params)} @footprints));
        }
        else {
            return "complex_unclear";
        }  # End of clusters with neither shards nor balanced breakpoint footprints

        # Default if the cluster has shards and/or balanced breakpoints but is not
        # a pure chain or cycle. 
        return "complex_unclear";
    }  # End of multi-footprint clusters
}


1;
