#!/usr/bin/perl
use strict;

# Sma3s predicted cripsr genes
my %sma;
open in, "./ab_crisprGenes_repeats.tsv";
while (<in>) {
  chomp;

  my ($a1, $a2) = split/\t/;
  $sma{$a1} = $a2;
}
close in;

my @f = `ls ./crispr_abs/ab*/TSV/CRISPR-Cas_summary.tsv`;
foreach my $f (@f) {
  my ($ab) = (split/\//, $f)[6];
  my $n1 = 0;
  my $n2 = 0;

  open in, $f;
  while (<in>) {
    chomp;
  
    next if /^Sequence\(s\)/;
    my ($ar, $gn) = (split/\t/)[2,5];
    $n1 += $gn; 
    $n2 += $ar; 
  }
  close in;
  #print "$ab\t$n1\t$n2\n";
  print "$ab\t$sma{$ab}\t$n2\n";
}
