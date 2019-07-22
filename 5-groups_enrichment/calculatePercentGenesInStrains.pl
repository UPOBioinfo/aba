#!/usr/bin/perl
use strict;

my $f1 = $ARGV[0]; # group 1
my $f2 = $ARGV[1]; # group 2
my $pan = $ARGV[2] || "./pan_gn.fasta2";
my %ff1;
my %ff2;
my $tf1 = 0;
my $tf2 = 0;

open in, $f1;
while (<in>) {
  chomp;

  $_ =~ s/ab//;
  $_ =~ s/^0+//;
  $ff1{$_} = 1;
  $tf1++;
}
close in;

open in, $f2;
while (<in>) {
  chomp;

  $_ =~ s/ab//;
  $_ =~ s/^0+//;
  $ff2{$_} = 1;
  $tf2++;
}
close in;

open in, $pan;
while (<in>) {
  chomp;

  my $nf1 = 0;
  my $nf2 = 0;
  my ($id, $ab) = (split/\t/)[0, 2];
  $id =~ s/^Cluster //;
  print "$id\t";
  foreach (split/,/, $ab) {
    $nf1++ if ($ff1{$_});
    $nf2++ if ($ff2{$_});
  }
  printf"%.2f\t%.2f", $nf1/$tf1, $nf2/$tf2;
  print "\n";
}
close in;

exit;
