#!/usr/bin/perl
use strict;

my @CL;
my %CL;
my @S;
my $file_ab = $ARGV[1] || "./strains.ab";

open in, $file_ab;
while (<in>) {
  chomp;

  $_ =~ s/ab//;
  $_ =~ s/^0+//;
  $S[$_] = 1;
}
close in;

my $FILE = $ARGV[0];
open FILE, $FILE;
while (<FILE>) {
  chomp;

  my ($g, @cl) = split/\t/;
  my ($idc, $ref) = split/: /, $g;
  foreach my $c ($ref, @cl) {
    ($c) = split/_/, $c;
    $c =~ s/ab//;
    $c =~ s/^0+//;
    push @{$CL{$idc}}, $c;
  }
}
close FILE;

# Read clusters
foreach my $c (keys %CL) {
#  print "Cluster $cl $i\n";
#for (my $i = 0; $i <= $#CL; $i++) {
  #print "Cluster $i\n";
  my $cl = $CL{$c};

# remove redundancy
  my %saw; @saw{@{$cl}} = (); @{$cl} = keys %saw;
  my @str;
  my @not;
  foreach my $genes (@{$cl}) {
    $str[$genes] = 1;
  }
  for (my $x = 1; $x <= 2467; $x++) {
    next unless ($S[$x]);
    push @not, $x if !$str[$x];
  }
    
  my $n = $#{$cl} + 1;
  print "Cluster $c\t$n\t";
  print join ",", @{$cl};
  print "\t";
  print join ",", @not;
  print "\n";
}

exit;
