# Footprint.pm
# A footprint is a genomic region with a high density of rearrangement
# breakpoints. These are defined externally, for example using clustering_index.R.
# A rearrangement cluster is defined as a set of footprints. 

use strict;
use FootprintCluster;
use List::Util qw(min max);
use warnings FATAL => 'all';
use Carp qw(carp);

package Footprint;

sub new {
    my $class = shift;
    my %params = @_;
    if (!exists($params{chrom})) {
        die "\$params{chrom} is required in $class\->new(%params)";
    }
    if (!exists($params{id})) {
        die "\$params{id} is required in $class\->new(%params)";
    }
    if (!exists($params{start})) {
        die "\$params{start} is required in $class\->new(%params)";
    }
    if (!exists($params{end})) {
        die "\$params{end} is required in $class\->new(%params)";
    }
    if (!exists($params{cluster})) {
        die "\$params{cluster} is required in $class\->new(%params)";
    }

    my $footprint = bless {
        chrom    => $params{chrom},
        start    => $params{start},
        end      => $params{end},
        id       => $params{id},       # ID of the footprint
        cluster  => $params{cluster},  # SV cluster of the current footprint
        rg_ends  => {},                # RearrangementEnds of the current footprint
    }, $class;

    return $footprint;
}

sub destroy {
    my $self = shift;

    # Remove cached rg pattern value
    delete $self->{rg_pattern};
    delete $self->cluster->footprints->{$self->id};
    my $cluster = $self->cluster;
    delete $self->{cluster};
    if ($cluster->size == 0) {
        $cluster->destroy;
    }
}

sub chrom {
    my $self = shift;
    return $self->{chrom};
}

sub chr_name {
    my $self = shift;
    return $self->chrom->name;
}

sub start {
    my $self = shift;
    return $self->{start};
}

sub end {
    my $self = shift;
    return $self->{end};
}

sub id {
    my $self = shift;
    return $self->{id};
}

sub cluster {
    my $self = shift;
    return $self->{cluster};
}

sub set_cluster {
    my $self = shift;
    my $new_cluster = shift;
    $self->{cluster} = $new_cluster;
}

sub add_rg_ends {
    my $self = shift;
    for (@_) {
        if ($_->chr_name ne $self->chr_name) {
            die sprintf("In $self\->add_rg_ends(), attempted to add below rearrangement end not matching footprint's chromosome %s:\n%s\n", $_->to_s, $self->chr_name);
        }
        # if (exists($self->rg_ends->{$_->rg->id . ":" . $_->end})) {
        if (defined($_->footprint)) {
            # warn sprintf("Warning: In $self\->add_rg_ends(), rearrangement end %s is also part of footprint %s\n", $_->to_s, $_->footprint->to_s);
            Carp::carp sprintf("Warning: In $self\->add_rg_ends(), rearrangement end %s is also part of footprint %s\n", $_->to_s, $_->footprint->to_s);
        }
        else {
            # $self->rg_ends->{$_->rg->id . ":" . $_->end} = $_;
            $self->rg_ends->{$_} = $_;
            $_->{footprint} = $self;
        }
    }

    # Remove cached rg pattern value
    delete $self->{rg_pattern};
}

sub remove_rg_ends {
    my $self = shift;
    my %params = @_;
    my($update_start, $update_end) = (0, 0);
    my @sorted_bkpts = $self->sorted_rg_ends_array(%params);
    for (@{$params{rg_ends}}) {
        if ($_ == $sorted_bkpts[0]) {
            $update_start = 1;
        }
        if ($_ == $sorted_bkpts[-1]) {
            $update_end = 1;
        }
        $self->{rg_ends}->{$_}->{footprint} = undef;
        delete $self->{rg_ends}->{$_};
    }
    @sorted_bkpts = $self->sorted_rg_ends_array(%params);
    if ($update_start and @sorted_bkpts) {
        $self->{start} = List::Util::min(map {$_->pos} @sorted_bkpts);
    }
    if ($update_end and @sorted_bkpts) {
        $self->{end} = List::Util::max(map {$_->pos} @sorted_bkpts);
    }

    if ($self->n_rg_ends == 0) {
        $self->destroy;
    }

    # Remove cached rg pattern value
    delete $self->{rg_pattern};
}

sub rg_ends {
    my $self = shift;
    return $self->{rg_ends};
}

sub rg_ends_array {
    my $self = shift;
    return values($self->rg_ends);
}

sub n_rg_ends {
    my $self = shift;
    return scalar($self->rg_ends_array);
}

sub sorted_rg_ends_array {
    my $self = shift;
    my %params = @_;
    my @rg_ends_array = sort { $a->pos <=> $b->pos  ||  ($a->is_fwd && $b->is_rev ? 1 : -1) } $self->rg_ends_array;
    my $i = 0;
    while ($i < $#rg_ends_array) {
        if (
                $rg_ends_array[$i]->is_rev
                and $rg_ends_array[$i]->is_bal_rg_overlap(%params)
                and $rg_ends_array[$i]->balanced_bkpt_partner_rg_end(%params) == $rg_ends_array[$i+1]
        ) {
            @rg_ends_array[$i, $i+1] = @rg_ends_array[$i+1, $i];
            $i += 2;
        }
        else {
            $i++;
        }
    }

    return @rg_ends_array;
}

sub rgs {
    my $self = shift;
    my %rgs;
    for ($self->rg_ends_array) {
        $rgs{$_->rg->id} = $_->rg;
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

sub size {
    my $self = shift;
    return $self->n_rg_ends;
}

sub first_pos {
    my $self = shift;
    return List::Util::min(map { $_->pos } $self->rg_ends_array);
}

sub last_pos {
    my $self = shift;
    return List::Util::max(map { $_->pos } $self->rg_ends_array);
}

sub to_s {
    my $self = shift;
    my %params = @_;
    return join(
        " ",
        $self->id,
        $self->chr_name,
        $self->start . "-" . $self->end,
        join(",", map({$_->pos . ($_->is_fwd ? "+" : "-")} $self->sorted_rg_ends_array(%params))),
        sprintf("%d breakpoints", $self->size)
    );
}

sub neighbour_footprints {
    my $self = shift;
    my @rg_ends = $self->rg_ends_array();
    my %neighbours;
    for (@rg_ends) {
        if ($_->mate->footprint eq $self) {
            next;
        }
        $neighbours{$_->mate->footprint} = $_->mate->footprint;
    }
    return values(%neighbours);
}


#
# Interpretation of footprints
#
sub type {
    # Known footprints:
    # "single"          - single breakpoint
    # "balanced"        - balanced breakpoint
    # "shard"           - templated insertion
    # "unphased_size_2" - two breakpoints in same direction
    # "unknown"         - everything else
    # * (chromoplexy|shard)_(chain|cycle):<length>
    # * Complex with gains (complex:<relative_cn>)
    # * Complex with loss (complex:<relative_cn>)
    # * 0+/2-/2+ (bkpt_with_shard:<relative_cn>)
    # * 1-/1+/2+ (shard_before_bkpt:<relative_cn>)
    # * 1-,2-/2+ (foldback_transloc:<relative_cn>)
    # * 0+,2-/2+ (del_shard_transloc:<relative_cn>)
    # * 0+,2+/2- (inv_shard_transloc:<relative_cn>)
    # * Other complex (complex)
    
    my $self = shift;
    my %params = @_;
    
    if ($self->size == 1) {
        return "single";
    }
    elsif (
            $self->cluster->size == 1 and
            $self->size == 2 * $self->n_rgs and
            scalar(grep {$_->is_td_type or $_->is_del_type} $self->rgs_array) == $self->n_rgs
    ) {
        return "simple_dels_and_tds";
    }
    elsif ($self->size == 2) {
        my ($rg_end_l, $rg_end_h) = $self->sorted_rg_ends_array(%params);

        if ($self->n_rgs == 1) {
            my ($rg) = $self->rgs_array;
            if ($rg->is_del_type(%params)) {
                return "del";
            }
            elsif ($rg->is_td_type(%params)) {
                return "td";
            }
        }

        if ($rg_end_l->rg->id eq $rg_end_h->rg->id) {
            return "simple_" . $rg_end_l->rg->rg_type_s(%params);
        }
        elsif (
                $rg_end_l->is_fwd and
                $rg_end_h->is_rev and
                $rg_end_l->segment->is_bal_rg_overlap(%params)
        ) {
            return "balanced_olap";
        }
        elsif ($rg_end_l->is_fwd and $rg_end_h->is_rev) {
            return "balanced";
        }
        elsif ($rg_end_l->is_rev and $rg_end_h->is_fwd) {
            return "shard";
        }
        else {
            return "unphased_size_2";
        }
    }
    elsif ($self->size == 3) {
        my $rg_pattern = $self->rg_pattern(%params);
        if (!defined($rg_pattern)) {
            return "complex_unclear";
        }
        elsif ($rg_pattern->{rg_string} eq "0+/2-/2+") {
            return "bkpt_with_shard";
        }
        elsif ($rg_pattern->{rg_string} eq "1-/1+/2+") {
            return "shard_before_bkpt";
        }
        elsif ($rg_pattern->{rg_string} eq "1-,2-/2+") {
            return "foldback_transloc";
        }
        elsif ($rg_pattern->{rg_string} eq "0+,2-/2+" ) {
            # This is impossible when deletions peeling is enabled. 
            return "del_shard_transloc";
        }
        elsif ($rg_pattern->{rg_string} eq "0+,2+/2-") {
            return "inv_shard_transloc";
        }
        else {
            return "complex_unclear";
        }
    }
    elsif ($self->size == 4 and scalar($self->n_rgs) == 2) {
        # Various types of local two-jumps
        my $rg_pattern = $self->rg_pattern(%params);

        if (!defined($rg_pattern)) {
            return "complex_unclear";
        }   
        elsif ($rg_pattern->{rg_string} eq "0+,2+/2-,4-") {
            return "direct_inversion";
        }   
        elsif ($rg_pattern->{rg_string} eq "0+,3+/2-,3-") {
            return "inversion_gain_loss";
        }   
        elsif ($rg_pattern->{rg_string} eq "1-,3-/1+,3+") {
            return "inverted_duplication";
        }   
        elsif ($rg_pattern->{rg_string} eq "1-,2-/2+,3+") {
            # Check whether copy number pattern is more consistent with 2xBFB
            # or dup-trp-dup.
            my @rgs = $self->rgs_array();
            my $rg_1_cn_change = $rgs[0]->weighted_avg_cn_change_across_rg(%params);
            my $rg_2_cn_change = $rgs[1]->weighted_avg_cn_change_across_rg(%params);
            if (
                    !defined($rg_1_cn_change)  or
                    !defined($rg_2_cn_change)  or
                    $rg_1_cn_change <= 0  or
                    $rg_2_cn_change <= 0
            ) {
                # Play it safe - no CN -> dup-trp-dup which is more frequent
                return "dup_trp_dup";
            }
            elsif (abs(log($rg_1_cn_change/$rg_2_cn_change)) > log(1.5)) {
                # At least 1.5-fold difference in copy number difference.
                return "fb_then_fb";
            }
            else {
                return "dup_trp_dup";
            }
        }   
        else {
            return "complex_unclear";
        }   
    }
    else {
        return "unknown";
    }
}

sub is_del_type {
    my $self = shift;
    my %params = @_;
    return $self->type(%params) eq "del";
}

sub is_td_type {
    my $self = shift;
    my %params = @_;
    return $self->type(%params) eq "td";
}

sub is_balanced_type {
    my $self = shift;
    my %params = @_;
    return($self->type(%params) eq "balanced" or $self->type(%params) eq "balanced_olap");
}

sub is_shard_type {
    my $self = shift;
    my %params = @_;
    return $self->type(%params) eq "shard";
}

sub avg_relative_cn {
    # Whether the rearrangement-side CN of the current footprint is
    # mostly consistent with a gain or a loss. 
    my $self = shift;
    my %params = @_;

    my $chrom = $self->chrom;
    my($p_arm_cn, $q_arm_cn);
    if (!defined($chrom->p_centromere_cn)) {
        $p_arm_cn = $chrom->p_arm_telomere_cn;
    }
    elsif (!defined($chrom->p_arm_telomere_cn)) {
        $p_arm_cn = $chrom->p_centromere_cn;
    }
    else {
        # $p_arm_cn = List::Util::max($chrom->p_centromere_cn, $chrom->p_arm_telomere_cn);
        $p_arm_cn = $chrom->p_centromere_cn/2 + $chrom->p_arm_telomere_cn/2;
    }
    if (!defined($chrom->q_centromere_cn)) {
        $q_arm_cn = $chrom->q_arm_telomere_cn;
    }
    elsif (!defined($chrom->q_arm_telomere_cn)) {
        $q_arm_cn = $chrom->q_centromere_cn;
    }
    else {
        # $q_arm_cn = List::Util::max($chrom->q_centromere_cn, $chrom->q_arm_telomere_cn);
        $q_arm_cn = $chrom->q_centromere_cn/2 + $chrom->q_arm_telomere_cn/2;
    }

    my($cn_diff_sum, $valid_bkpt_count) = (0, 0);
    for my $rg_end (map {($_->low_end, $_->high_end)} ($self->rgs_array)) {
        if (defined($rg_end->rg_side_cn(%params))) {
            if ($rg_end->is_on_p_arm and defined($p_arm_cn)) {
                $cn_diff_sum += $rg_end->rg_side_cn(%params) - $p_arm_cn;
            }
            elsif ($rg_end->is_on_q_arm and defined($q_arm_cn)) {
                $cn_diff_sum += $rg_end->rg_side_cn(%params) - $q_arm_cn;
            }
            $valid_bkpt_count++;
        }
    }

    if ($valid_bkpt_count >= 0) {
        return $cn_diff_sum / $valid_bkpt_count;
    }
    else {
        return undef;
    }
}

sub detailed_type {
    # Detailed type is for example balanced rearrangement vs. chromoplexy
    #

    my $self = shift;
    my %params = @_;
    if ($self->type(%params) eq "single") {
        return "single";
    }
    elsif ($self->is_del_type(%params) or $self->is_td_type(%params)) {
        return $self->type(%params);
    }
    elsif ($self->is_balanced_type(%params)) {
        my($bkpt_1, $bkpt_2) = $self->sorted_rg_ends_array(%params);
        my $starting_bkpt = $bkpt_1;
        my $footprint_count = 1;
        while ($bkpt_1->mate->footprint->is_balanced_type(%params)) {
            $bkpt_1 = $bkpt_1->mate->balanced_bkpt_partner_rg_end_by_footprint(%params);
            if ($bkpt_1 == $starting_bkpt) {
                return "chromoplexy_cycle:$footprint_count";
            }
            $footprint_count++;
        }
        while ($bkpt_2->mate->footprint->is_balanced_type(%params)) {
            $bkpt_2 = $bkpt_2->mate->balanced_bkpt_partner_rg_end_by_footprint(%params);
            $footprint_count++;
        }
        return "chromoplexy_chain:$footprint_count";

        # Otherwise it's just a balanced breakpoint
        return "balanced_breakpoint";
    }  # End of balanced breakpoint type footprint
    elsif ($self->is_shard_type(%params)) {
        my($bkpt_1, $bkpt_2) = $self->sorted_rg_ends_array(%params);
        my $starting_bkpt = $bkpt_1;
        my $footprint_count = 1;
        while ($bkpt_1->mate->footprint->is_shard_type(%params)) {
            $bkpt_1 = $bkpt_1->mate->shard_partner_rg_end_by_footprint(%params);
            if ($bkpt_1 == $starting_bkpt) {
                return "shard_cycle:$footprint_count";
            }
            $footprint_count++;
        }
        while ($bkpt_2->mate->footprint->is_shard_type(%params)) {
            $bkpt_2 = $bkpt_2->mate->shard_partner_rg_end_by_footprint(%params);
            $footprint_count++;
        }
        return "shard_chain:$footprint_count";

        # Otherwise it's just a templated insertion
        return "shard";
    }  # End of shard type footprint
    elsif ($self->size == 3) {
        my $type = $self->type(%params);
        my $footprint_relative_cn = $self->avg_relative_cn(%params);
        $footprint_relative_cn = "NA" if !defined($footprint_relative_cn);
        return sprintf("$type:%.2f", $footprint_relative_cn);
    }  # End of size 3 footprints
    elsif ($self->size == 4 and $self->n_rgs == 2) {
        return $self->type(%params);
    }
    else {
        # Four or more SVs within a single footprint
        # Compute average breakpoint CN compared to (maximum estimate) of arm-level copy number
        # my $avg_relative_cn = $self->avg_relative_cn(%params);
        # $avg_relative_cn = "NA" if !defined($avg_relative_cn);
        # return sprintf("complex_unclear:%.2f", $avg_relative_cn);
        return $self->type(%params);
    }
}

sub rg_pattern {
    # The rearrangement pattern of a pattern ignores all rearrangement
    # breakpoints outside the footprint. Therefore, there is only two ways to
    # represent a rearrangement pattern of a footprint, forward or reverse.
    #
    # Brief description of the algorithm:
    # 1. The footprint is traversed in forward and reverse orientation, yielding
    #    two rearrangement pattern strings. The lexicographically smaller string
    #    is returned.
    # 2. As rearrangements are encountered the first time when the footprint is
    #    traversed, they are pushed to @rg_buckets, so that when the second
    #    breakpoint of the rearrangements are encountered, the second breakpoint
    #    segment location can be stored in the same bucket of its first
    #    breakpoint.
    # 3. As the footprint is traversed, $seg_idx keeps track of the current
    #    segment index. 
    # 4. Copy number changes across breakpoints and their relative variances
    #    are stored along the way. These need to be stored only in one direction
    #    as they can be simply reversed. 
    #
    # Parameter 'symbol_seg_names' sets copy number segments as 'A', 'B' etc.
    # as opposed to '0', '1' etc.

    my $self = shift;
    my %params = @_;
    if (!exists($params{symbol_seg_names})) {
        $params{symbol_seg_names} = 0;
    }

    sub normalise_rg_string {
        my $use_symbol_seg_names = shift;
        if (!$use_symbol_seg_names) {
            return(join(
                "/",
                map( {join(",", @{$_})} @_ )
            ))
        }
        else {
            return(join(
                "/",  # Separator between SVs
                map(
                    {join(
                        "^",
                        map(
                            {chr(substr($_, 0, -1)+65) . substr($_, -1)}
                            @{$_}
                        )
                    )}
                    @_
                )
            ));
        }
    }

    # Set 'no_bal_bkpts' as default - make a gap between all 'balanced' breakpoints
    if (!exists($params{no_bal_bkpts})) {
        $params{no_bal_bkpts} = 1;
    }

    # Don't compute for footprints with more than 20 breakpoints
    if ($self->size >= 20) {
        return undef;
    }

    # Cached value exists?
    if (exists($self->{rg_pattern})) {
        return $self->{rg_pattern}->{
            ($params{symbol_seg_names} ? 'use_symbols' : 'use_numbers')
        };
    }

    my($bkpt, $next_bkpt, $prev_bkpt);
    my($seg, $seg_idx);
    my(%bucket_idx_of_rg, @rg_buckets, $bucket_idx);
    my(@cn_change_across_bkpt, @cn_change_relative_var, @bkpt_cns, @bkpt_pos, @relative_cns);

    sub cache_result {
        # Cache the result and return it.
        my $self = shift;
        $self->{rg_pattern}->{use_symbols} = {
            %{$_[0]},
            rg_string => normalise_rg_string(1, @{$_[0]->{rg_buckets}})
        };
        delete($self->{rg_pattern}->{use_symbols}->{symbol_seg_names});
        $self->{rg_pattern}->{use_numbers} = {
            %{$_[0]},
            rg_string => normalise_rg_string(0, @{$_[0]->{rg_buckets}})
        };
        delete($self->{rg_pattern}->{use_numbers}->{symbol_seg_names});
        if ($_[0]->{symbol_seg_names}) {
            return $self->{rg_pattern}->{use_symbols};
        }
        else {
            return $self->{rg_pattern}->{use_numbers};
        }
    }

    # Special case - if multiple breakpoints have the exact same orientation
    # and breakpoint, then simply just skip. Special case: if two breakpoints
    # come from the same rearrangement, then it's just a fold-back and we can
    # keep it.
    my %pos_and_dir;
    for ($self->rg_ends_array) {
        $pos_and_dir{$_->pos . ($_->is_fwd ? "+" : "-")}->{$_->rg->id} = 1;
    }
    for my $p (values (%pos_and_dir)) {
        if (scalar(keys %{$p}) > 1) {
            return undef;
        }
    } 
    
    my @sorted_bkpts = $self->sorted_rg_ends_array(%params);

    # Sanity check
    if (
        # $sorted_bkpts[0]->is_rev &&
        ($sorted_bkpts[0])->is_bal_rg_overlap(%params)
    ) {
        # This is a balanced overlap breakpoint, so the balanced breakpoint
        # partner better be in the same footprint.
        if (
            @sorted_bkpts < 2 ||
            # $sorted_bkpts[1]->is_rev ||
            !$sorted_bkpts[1]->is_bal_rg_overlap(%params) ||
            $sorted_bkpts[1] != $sorted_bkpts[0]->balanced_bkpt_partner_rg_end(%params)
        ) {
            warn(sprintf(
                "Breakpoint's (%s) balanced overlap mate in a different " .
                    "footprint. Footprint has rgs %s.",
                 $sorted_bkpts[0]->to_s(%params),
                 join(", ", map({$_->id} $self->rgs_array)))
            );
            return undef;
        }
    }
    if (
        # $sorted_bkpts[-1]->is_fwd &&
        $sorted_bkpts[-1]->is_bal_rg_overlap(%params)
    ) {
        # This is a balanced overlap breakpoint, so the balanced breakpoint
        # partner better be in the same footprint.
        if (
            @sorted_bkpts < 2 ||
            # $sorted_bkpts[-2]->is_fwd ||
            !$sorted_bkpts[-2]->is_bal_rg_overlap(%params) ||
            $sorted_bkpts[-2] != $sorted_bkpts[-1]->balanced_bkpt_partner_rg_end(%params)
        ) {
            warn(sprintf(
                "Breakpoint's (%s) balanced overlap mate in a different " .
                    "footprint. Footprint has rgs %s.",
                 $sorted_bkpts[-1]->to_s(%params),
                 join(", ", map({$_->id} $self->rgs_array)))
            );
            return undef;
        }
    }

    # Forward rearrangement pattern
    $seg = $sorted_bkpts[0]->segment;
    $seg_idx = ($sorted_bkpts[0]->is_fwd ? 0 : 1);
    while (1) {
        # Add low end breakpoint if it's part of the footprint
        $bkpt = $seg->low_end->bkpt;
        if (defined($bkpt) && exists($self->rg_ends->{$bkpt})) {
            if ($bkpt->is_bal_rg_overlap(%params)) {
                # Sanity checks
                if (
                    !defined($bkpt->segment->high_end->bkpt) ||
                    $bkpt->balanced_bkpt_partner_rg_end(%params) != $bkpt->segment->high_end->bkpt
                ) {
                    warn(sprintf(
                        "Breakpoint's (%s) balanced overlap mate in a different " .
                            "footprint. Footprint has rgs %s.",
                         $bkpt->to_s(%params),
                         join(", ", map({$_->id} $self->rgs_array)))
                    );
                    return undef;
                }

                # Process the fwd end of the balanced breakpoint first
                $next_bkpt = $bkpt->segment->high_end->bkpt;
                if ($seg_idx > 0) {
                    $seg_idx--;
                }
                push @cn_change_across_bkpt, "NA";
                push @cn_change_relative_var, "NA";
                push @bkpt_pos, $next_bkpt->pos;
                push @bkpt_cns, $next_bkpt->rg_side_cn(%params);
                push @relative_cns, $next_bkpt->cn_relative_to_arm(%params);
                if (exists($bucket_idx_of_rg{$next_bkpt->id})) {
                    $bucket_idx = $bucket_idx_of_rg{$next_bkpt->id};
                    if (scalar(@{$rg_buckets[$bucket_idx]}) != 1) {
                        die;  # Sanity check
                    }
                    push @{$rg_buckets[$bucket_idx]}, "$seg_idx+";
                }
                else {
                    push @rg_buckets, ["$seg_idx+"];
                    $bucket_idx_of_rg{$next_bkpt->id} = $#rg_buckets;
                }

                if (!$params{no_bal_bkpts}) {
                    $seg_idx++;
                }
                else {
                    $seg_idx += 2;
                }

                # Then process the rev end
                push @cn_change_across_bkpt, "NA";
                push @cn_change_relative_var, "NA";
                push @bkpt_pos, $bkpt->pos;
                push @bkpt_cns, $bkpt->rg_side_cn(%params);
                push @relative_cns, $bkpt->cn_relative_to_arm(%params);
                if (exists($bucket_idx_of_rg{$bkpt->id})) {
                    $bucket_idx = $bucket_idx_of_rg{$bkpt->id};
                    if (scalar(@{$rg_buckets[$bucket_idx]}) != 1) {
                        die;  # Sanity check
                    }
                    push @{$rg_buckets[$bucket_idx]}, "$seg_idx-";
                }
                else {
                    push @rg_buckets, ["$seg_idx-"];
                    $bucket_idx_of_rg{$bkpt->id} = $#rg_buckets;
                }
                $seg = $seg->next_seg;

                $prev_bkpt = $bkpt;
            }
            else {
                # We end up here if $bkpt is not a balanced overlap breakpoint

                if (
                    $bkpt->is_bal_rg_overlap(%params) &&
                    !exists($self->rg_ends->{$bkpt->balanced_bkpt_partner_rg_end(%params)})
                ) {
                    warn(sprintf(
                        "Breakpoint's (%s) balanced overlap mate in a different " .
                            "footprint. Footprint has rgs %s.",
                         $bkpt->to_s(%params),
                         join(", ", map({$_->id} $self->rgs_array)))
                    );
                    return undef;
                }

                # These are rearrangements that are not balanced overlap breakpoints
                if (exists($bucket_idx_of_rg{$bkpt->id})) {
                    $bucket_idx = $bucket_idx_of_rg{$bkpt->id};
                    if (scalar(@{$rg_buckets[$bucket_idx]}) != 1) {
                        die;  # Sanity check
                    }
                    push @{$rg_buckets[$bucket_idx]}, ("$seg_idx-");
                }
                else {
                    push @rg_buckets, ["$seg_idx-"];
                    $bucket_idx_of_rg{$bkpt->id} = $#rg_buckets;
                }

                push @cn_change_across_bkpt, ($bkpt->cn_across_bkpt(%params) or "NA");
                push @cn_change_relative_var, ($bkpt->relative_cn_diff_var(%params) or "NA");
                push @bkpt_pos, $bkpt->pos;
                push @bkpt_cns, $bkpt->rg_side_cn(%params);
                push @relative_cns, $bkpt->cn_relative_to_arm(%params);

                $prev_bkpt = $bkpt;
            }
        }

        # Add high end breakpoint if it's part of the footprint
        $bkpt = $seg->high_end->bkpt;
        if (defined($bkpt) && exists($self->rg_ends->{$bkpt})) {
            if (exists($bucket_idx_of_rg{$bkpt->id})) {
                $bucket_idx = $bucket_idx_of_rg{$bkpt->id};
                if (scalar(@{$rg_buckets[$bucket_idx]}) != 1) {
                    die;  # Sanity check
                }
                push @{$rg_buckets[$bucket_idx]}, ("$seg_idx+");
            }
            else {
                push @rg_buckets, ["$seg_idx+"];
                $bucket_idx_of_rg{$bkpt->id} = $#rg_buckets;
            }

            if ($bkpt->has_balanced_bkpt_partner_rg_end(%params)) {
                # Sanity check
                if (!exists($self->rg_ends->{$bkpt->balanced_bkpt_partner_rg_end(%params)})) {
                    warn(sprintf(
                        "Breakpoint's (%s) balanced overlap mate in a different " .
                            "footprint. Footprint has rgs %s.",
                         $bkpt->to_s(%params),
                         join(", ", map({$_->id} $self->rgs_array)))
                    );
                    return undef;
                }
                $seg = $seg->next_seg;
            }

            push @cn_change_across_bkpt, ($bkpt->cn_across_bkpt(%params) or "NA");
            push @cn_change_relative_var, ($bkpt->relative_cn_diff_var(%params) or "NA");
            push @bkpt_pos, $bkpt->pos;
            push @bkpt_cns, $bkpt->rg_side_cn(%params);
            push @relative_cns, $bkpt->cn_relative_to_arm(%params);

            $prev_bkpt = $bkpt;
        }

        # Go to the next segment involving a breakpoint from the footprint
        if (!defined($seg->next_seg)) {
            last;
        }
        $seg = $seg->next_seg;
        while (
            !(defined($seg->low_end->bkpt) && exists($self->rg_ends->{$seg->low_end->bkpt})) &&
            !(defined($seg->high_end->bkpt) && exists($self->rg_ends->{$seg->high_end->bkpt})) &&
            defined($seg->next_seg)
        ) {
            $seg = $seg->next_seg;
        }

        if (
            !(defined($seg->low_end->bkpt) && exists($self->rg_ends->{$seg->low_end->bkpt})) &&
            !(defined($seg->high_end->bkpt) && exists($self->rg_ends->{$seg->high_end->bkpt})) &&
            !defined($seg->next_seg)
        ) {
            last;
        }

        if ($prev_bkpt->is_rev) {
            # Only need to shift if the next breakpoint is rev too
            if (defined($seg->low_end->bkpt) && exists($self->rg_ends->{$seg->low_end->bkpt})) {
                $seg_idx++;
            }
        }
        else {
            if (defined($seg->low_end->bkpt) && exists($self->rg_ends->{$seg->low_end->bkpt})) {
                $seg_idx += 2;
            }
            else {
                $seg_idx++;
            }
        }
    }
    # my $fwd_rg_pattern = join(
    #     "/",
    #     map( {join(",", @{$_})} @rg_buckets )
    # );
    my @fwd_rg_buckets = @rg_buckets;

    # Now generate the reverse rg pattern
    @rg_buckets = %bucket_idx_of_rg = ();
    $seg = $sorted_bkpts[-1]->segment;
    $seg_idx = ($sorted_bkpts[-1]->is_rev ? 0 : 1);
    while (1) {
        # Add high end breakpoint if it's part of the footprint
        $bkpt = $seg->high_end->bkpt;
        if (defined($bkpt) && exists($self->rg_ends->{$bkpt})) {
            if ($bkpt->is_bal_rg_overlap(%params)) {
                # Sanity checks
                if (
                    !defined($bkpt->segment->low_end->bkpt) ||
                    $bkpt->balanced_bkpt_partner_rg_end(%params) != $bkpt->segment->low_end->bkpt
                ) {
                    warn(sprintf(
                        "Breakpoint's (%s) balanced overlap mate in a different " .
                            "footprint. Footprint has rgs %s.",
                         $bkpt->to_s(%params),
                         join(", ", map({$_->id} $self->rgs_array)))
                    );
                    return undef;
                }

                # Process the rev end of the balanced breakpoint first
                $next_bkpt = $bkpt->segment->low_end->bkpt;
                if ($seg_idx > 0) {
                    $seg_idx--;
                }
                if (exists($bucket_idx_of_rg{$next_bkpt->id})) {
                    $bucket_idx = $bucket_idx_of_rg{$next_bkpt->id};
                    if (scalar(@{$rg_buckets[$bucket_idx]}) != 1) {
                        die;  # Sanity check
                    }
                    push @{$rg_buckets[$bucket_idx]}, "$seg_idx+";
                }
                else {
                    push @rg_buckets, ["$seg_idx+"];
                    $bucket_idx_of_rg{$next_bkpt->id} = $#rg_buckets;
                }
                
                if (!$params{no_bal_bkpts}) {
                    $seg_idx++;
                }
                else {
                    $seg_idx += 2;
                }

                # Then process the fwd end
                if (exists($bucket_idx_of_rg{$bkpt->id})) {
                    $bucket_idx = $bucket_idx_of_rg{$bkpt->id};
                    if (scalar(@{$rg_buckets[$bucket_idx]}) != 1) {
                        die;  # Sanity check
                    }
                    push @{$rg_buckets[$bucket_idx]}, "$seg_idx-";
                }
                else {
                    push @rg_buckets, ["$seg_idx-"];
                    $bucket_idx_of_rg{$bkpt->id} = $#rg_buckets;
                }
                $seg = $seg->prev_seg;

                $prev_bkpt = $bkpt;
            }
            else {
                if (
                    $bkpt->is_bal_rg_overlap(%params) &&
                    !exists($self->rg_ends->{$bkpt->balanced_bkpt_partner_rg_end(%params)})
                ) {
                    warn(sprintf(
                        "Breakpoint's (%s) balanced overlap mate in a different " .
                            "footprint. Footprint has rgs %s.",
                         $bkpt->to_s(%params),
                         join(", ", map({$_->id} $self->rgs_array)))
                    );
                    return undef;
                }

                # These are rearrangements that are not balanced overlap breakpoints
                if (exists($bucket_idx_of_rg{$bkpt->id})) {
                    $bucket_idx = $bucket_idx_of_rg{$bkpt->id};
                    if (scalar(@{$rg_buckets[$bucket_idx]}) != 1) {
                        die;  # Sanity check
                    }
                    push @{$rg_buckets[$bucket_idx]}, ("$seg_idx-");
                }
                else {
                    push @rg_buckets, ["$seg_idx-"];
                    $bucket_idx_of_rg{$bkpt->id} = $#rg_buckets;
                }

                $prev_bkpt = $bkpt;
            }
        }

        # Add low end breakpoint if it's part of the footprint
        $bkpt = $seg->low_end->bkpt;
        if (defined($bkpt) && exists($self->rg_ends->{$bkpt})) {
            if (exists($bucket_idx_of_rg{$bkpt->id})) {
                $bucket_idx = $bucket_idx_of_rg{$bkpt->id};
                if (scalar(@{$rg_buckets[$bucket_idx]}) != 1) {
                    die;  # Sanity check
                }
                push @{$rg_buckets[$bucket_idx]}, ("$seg_idx+");
            }
            else {
                push @rg_buckets, ["$seg_idx+"];
                $bucket_idx_of_rg{$bkpt->id} = $#rg_buckets;
            }

            if ($bkpt->has_balanced_bkpt_partner_rg_end(%params)) {
                # Sanity check
                if (!exists($self->rg_ends->{$bkpt->balanced_bkpt_partner_rg_end(%params)})) {
                    warn(sprintf(
                        "Breakpoint's (%s) balanced overlap mate in a different " .
                            "footprint. Footprint has rgs %s.",
                         $bkpt->to_s(%params),
                         join(", ", map({$_->id} $self->rgs_array)))
                    );
                    return undef;
                }
                $seg = $seg->prev_seg;
            }

            $prev_bkpt = $bkpt;
        }

        # Go to the next segment involving a breakpoint from the footprint
        if (!defined($seg->prev_seg)) {
            last;
        }
        $seg = $seg->prev_seg;
        while (
            !(defined($seg->high_end->bkpt) && exists($self->rg_ends->{$seg->high_end->bkpt})) &&
            !(defined($seg->low_end->bkpt) && exists($self->rg_ends->{$seg->low_end->bkpt})) &&
            defined($seg->prev_seg)
        ) {
            $seg = $seg->prev_seg;
        }

        if (
            !(defined($seg->high_end->bkpt) && exists($self->rg_ends->{$seg->high_end->bkpt})) &&
            !(defined($seg->low_end->bkpt) && exists($self->rg_ends->{$seg->low_end->bkpt})) &&
            !defined($seg->prev_seg)
        ) {
            last;
        }

        if ($prev_bkpt->is_fwd) {
            # Only need to shift if the next breakpoint is rev too
            if (defined($seg->high_end->bkpt) && exists($self->rg_ends->{$seg->high_end->bkpt})) {
                $seg_idx++;
            }
        }
        else {
            if (defined($seg->high_end->bkpt) && exists($self->rg_ends->{$seg->high_end->bkpt})) {
                $seg_idx += 2;
            }
            else {
                $seg_idx++;
            }
        }
    }
    # my $rev_rg_pattern = join(
    #     "/",
    #     map( {join(",", @{$_})} @rg_buckets )
    # );
    my @rev_rg_buckets = @rg_buckets;

    sub diff {
        my @values = @_;
        if (@values < 2) {
            return();
        }
        else {
            return map({ abs($values[$_+1] - $values[$_]) } 0..($#_-1));
        }
    }

    # Sanity check
    if (scalar($self->rgs_array) != scalar(@rg_buckets)) {
        die;
    }

    # if ($rev_rg_pattern lt $fwd_rg_pattern) {
    if (normalise_rg_string(0, @rev_rg_buckets) lt normalise_rg_string(0, @fwd_rg_buckets)) {
        return cache_result($self, {
            # rg_string => normalise_rg_string(
            #     $params{symbol_seg_names},
            #     @rev_rg_buckets
            # ),
            symbol_seg_names => $params{symbol_seg_names},
            rg_buckets => [@rev_rg_buckets],
            n_segs => $seg_idx + 1,
            is_inverted => 1,
            cn_changes => [reverse(@cn_change_across_bkpt)],
            cn_change_variances => [reverse(@cn_change_relative_var)],
            distances => [diff(reverse(@bkpt_pos))],
            cns => [reverse(@bkpt_cns)],
            relative_cns => [reverse(@relative_cns)],
        });
    }
    else {
        return cache_result($self, {
            # rg_string => normalise_rg_string(
            #     $params{symbol_seg_names},
            #     @fwd_rg_buckets
            # ),
            symbol_seg_names => $params{symbol_seg_names},
            rg_buckets => [@fwd_rg_buckets],
            n_segs => $seg_idx + 1,
            is_inverted => 0,
            cn_changes => [@cn_change_across_bkpt],
            cn_change_variances => [@cn_change_relative_var],
            distances => [diff(@bkpt_pos)],
            cns => [reverse(@bkpt_cns)],
            relative_cns => [@relative_cns],
        });
    }
}

sub simple_rg_pattern_from_rgs {
    # Given an array of rearrangements on the same chromosome, treat all
    # the breakpoints of the rearrangements as coming from a single
    # footprint, then get the rearrangement pattern of the footprint. 
    my $class = shift;
    my @rgs = @{shift()};
    my %params = @_;

    # Sanity check
    my %chrs = ();
    for (@rgs) {
        $chrs{$_->low_end->segment->chr} = $_->low_end->segment->chr;
        $chrs{$_->high_end->segment->chr} = $_->high_end->segment->chr;
    }
    if (scalar(keys %chrs) > 1) {
        die "Attempted to run Footprint->simple_rg_pattern_from_rgs() with "
            . "rearrangements on multiple chromosomes.";
    }

    my $footprint = Footprint->new(
        chrom => (values %chrs)[0],
        start => "p_tel",
        end => "q_tel",
        id => "afdjga;jga;erjga;oirgja;rigja;gi",  # Totally randomised ID
        cluster => FootprintCluster->new(
            id => "feagrgjar;gjia;reg",
            genome => $rgs[0]->low_end->footprint->cluster->{genome}
        ),
    );

    # Store the original footprint IDs of the rearrangements and set them to
    # the temporary one. 
    my %orig_footprint_of;
    for (@rgs) {
        $orig_footprint_of{$_->low_end} = $_->low_end->footprint;
        $_->low_end->{footprint} = undef;  # To suppress warning
        $orig_footprint_of{$_->high_end} = $_->high_end->footprint;
        $_->high_end->{footprint} = undef;  # To suppress warning
        $footprint->add_rg_ends($_->low_end, $_->high_end);
    }

    my $rg_pattern = $footprint->rg_pattern(%params);

    # Reset the footprint pointers of the rearrangement ends
    for (@rgs) {
        $_->low_end->{footprint} = $orig_footprint_of{$_->low_end};
        $_->high_end->{footprint} = $orig_footprint_of{$_->high_end};
    }

    return $rg_pattern;
}


1;
