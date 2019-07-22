#!/usr/bin/perl
use strict;

my %cl;
my @pl;
my @ab;
my %m;
my %mob;
my %mobref;
my %mobref2;

# MOB-Cluster
open in, "./mobclustering/mob_plasmids_all/clusters.txt"; #cdhit+2ªcolumn of MOB
while (<in>) {
  chomp;

  next if /^id\t/;
  my ($id, $c) = (split/\t/)[0, 1]; # column 0.05 or 0.001
  if ($mobref{$c}) {
    $cl{$id} = $mobref{$c};
  } else {
    $mobref{$c} = $id;
    $cl{$id} = $id;
    push @pl, $id;
  }
}
close in;

# Gather ab identifiers
open ab, "./strains_ab.ab";
while (<ab>) {
    chomp;
    push @ab, $_; 
}
close ab;

# Plasmids Blasts
foreach my $ab (@ab) {
  open in, "./prokka/$ab/$ab.fna_plasmids.tsv";
  while (<in>) {
    chomp;
    
    my ($pl) = (split/\t/)[1];
    $m{$ab}{$pl} = 1;
  }
  close in;
}

# Create matrix
print "\t";
print join "\t", @pl;
print "\n";
foreach my $ab (@ab) {
  print $ab;
  foreach my $pl (@pl) {
    if ($m{$ab}{$pl}) {
      print "\t1";
    } else {
      print "\t";
    }
  }
  print "\n";
}

exit;
