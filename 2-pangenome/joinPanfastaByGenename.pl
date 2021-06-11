#!/usr/bin/perl
use strict;

# Collapse orthologous clusters from roary which have the same gene name from Sma3s annotation
# INPUT1: pancluster from roary in tsv format
# INPUT2: Sma3s annotation in tsv format
# OUTPUT1: collapsed annotation file in Sma3s format (a 2 is added to the output filename)
# OUTPUT2: collpsed pangenome in the same format than pancluster from roary (a 2 is added to the output filename)
# AJPerez, 2019 (updated June 2021)

my %pf;
my %gn;
my %genes;
my $un = 0;

if (!$ARGV[1]) { die "Please run as: ./joinPanfastaByGenename.pl pan_clusters.tsv pan_uniprot_bacteria_go.tsv\n" }
my $PC_FILE = $ARGV[0];
my $AN_FILE = $ARGV[1];
my ($a1, $a2) = split/\./, $PC_FILE;
my $PC_FILE2 = $a1 . "2." . $a2;
my ($b1, $b2) = split/\./, $AN_FILE;
my $AN_FILE2 = $b1 . "2." . $b2;

# panfasta
open in, $PC_FILE;
while (<in>) {
  chomp;

  my ($pf, @abs) = split/\t/;
  my ($gn, $pf) = split/: /, $pf;
  @{$pf{$pf}} = @abs;
  $gn{$pf} = $gn;
}
close in;

# sma3s
open in, $AN_FILE;
while (<in>) {
  chomp;

  next if /^#/;
  my ($gn) = (split/\t/)[1];
  if (!$gn) {
    $un++;
    $gn = "unknown$un";
  }
  push @{$genes{$gn}}, $_;
}
close in;

# construct panfasta
open ANNOT, ">$AN_FILE2";
foreach my $g (keys %genes) {
  my ($ab, $de) = (split/\t/, $genes{$g}[0])[0, 2];
  print ANNOT "$ab\t$g";
  for (my $x = 2; $x <= 11; $x++) {
    print ANNOT "\t";
    my @total = ();
    for (my $y = 0; $y <= $#{$genes{$g}}; $y++) {
      my (@l) = split/\t/, $genes{$g}[$y];
      my (@an) = split/;/, $l[$x];
      @total = (@total, @an);
    }
    my %hash = map { $_, 1 } @total; @total = keys %hash;
    if ($x == 2) {
      print ANNOT "$total[0]"; # the first description
    } else {
      print ANNOT join ";", @total;
    }
  }
  print ANNOT "\n";
}
close ANNOT;

# cdhit
my $ref;
open PAN, ">$PC_FILE2";
foreach my $g (keys %genes) {
  for (my $y = 0; $y <= $#{$genes{$g}}; $y++) {
    my ($ab) = split/\t/, $genes{$g}[$y];
    if ($y == 0) {
      print PAN "$g: $ab";
    } else {
      print PAN "\t$ab";
    }
    print PAN "\t" if (@{$pf{$ab}});
    print PAN join "\t", @{$pf{$ab}};
  }
  print PAN "\n";
}
close PAN;

print "Created files: $PC_FILE2 and $AN_FILE2\n";

exit;
