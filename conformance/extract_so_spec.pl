#!/usr/bin/env perl
# Extract Ensembl VEP's authoritative consequence specification — the %OVERLAP_CONSEQUENCES table
# in Bio::EnsEMBL::Variation::Utils::Constants — into conformance/data/so_consequences.tsv.
#
# This is THE operational spec our engine ports: each SO term with its SO accession, IMPACT,
# severity rank, tier, target feature, and the VariationEffect.pm `predicate` sub that decides
# whether the term applies. The Sequence Ontology gives the textual definitions; this table gives
# the operational semantics VEP actually computes. We parse the module as text (no Perl module
# loading, so no ensembl-core dependency chain).
#
# Usage: conformance/extract_so_spec.pl [path/to/Constants.pm] > conformance/data/so_consequences.tsv
use strict; use warnings;
my $path = $ARGV[0] || '/root/ensembl-variation/modules/Bio/EnsEMBL/Variation/Utils/Constants.pm';
open my $fh, '<', $path or die "cannot open $path: $!";
local $/; my $src = <$fh>;
my @cols = qw(SO_term SO_accession impact rank tier feature_SO_term predicate);
my @rows;
while ($src =~ /->new_fast\(\{(.*?)\}\s*\)/sg) {
    my $blk = $1; my %f;
    for my $c (@cols) { $f{$c} = ($blk =~ /'$c'\s*=>\s*'([^']*)'/) ? $1 : ''; }
    next unless $f{SO_term};
    $f{predicate} =~ s/.*:://;                       # keep the bare sub name
    push @rows, \%f;
}
print join("\t", @cols), "\n";
for my $r (sort { ($a->{rank} || 999) <=> ($b->{rank} || 999) } @rows) {
    print join("\t", map { $r->{$_} } @cols), "\n";
}
