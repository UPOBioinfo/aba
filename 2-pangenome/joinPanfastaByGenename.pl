#!/usr/bin/perl
use strict;

my %pf;
my %gn;
my %genes;
my $un = 0;

# panfasta
open in, $ARGV[0];
while (<in>) {
  chomp;

  my ($pf, @abs) = split/\t/;
  my ($gn, $pf) = split/: /, $pf;
  @{$pf{$pf}} = @abs;
  $gn{$pf} = $gn;
}
close in;

# sma3s
open in, $ARGV[1];
while (<in>) {
  chomp;

  my ($gn) = (split/\t/)[1];
  if (!$gn) {
    $un++;
    $gn = "unknown$un";
  }
  push @{$genes{$gn}}, $_;
}
close in;

# construct panfasta
open ANNOT, "./>pangenome_annot2.tsv";
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

# cdhit
my $ref;
foreach my $g (keys %genes) {
  for (my $y = 0; $y <= $#{$genes{$g}}; $y++) {
    my ($ab) = split/\t/, $genes{$g}[$y];
    if ($y == 0) {
      print "$g: $ab";
    } else {
      print "\t$ab";
    }
    print "\t" if (@{$pf{$ab}});
    print join "\t", @{$pf{$ab}};
  }
  print "\n";
}


