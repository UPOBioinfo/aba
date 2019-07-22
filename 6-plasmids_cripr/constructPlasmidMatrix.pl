#!/usr/bin/perl
use strict;

my %cl;
my @pl;
my @ab;
my %m;

# CD-Hit references
my $ref;
open in, "./cdh_AcinPlasmids.clstr";
while (<in>) {
  chomp;

  if (/>(.+)\.\.\. (.)/) {
    if ($2 eq "*") {
      $ref = $1 if ($2 eq "*");
      push @pl, $ref;
    }
    $cl{$1} = $ref;
  }
}
close in;

# Gather ab identifiers
open ab, "./strains.ab";
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
