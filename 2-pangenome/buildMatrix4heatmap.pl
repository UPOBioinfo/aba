#!/usr/bin/perl
use strict;

my %S;
my @par;
my $file_ab = $ARGV[1] || "./strains.ab";

open in, $file_ab;
while (<in>) {
  chomp;

  $_ =~ s/ab//;
  $_ =~ s/^0+//;
  $S{$_} = $_;
}
close in;

open in, $ARGV[0];
while (<in>) {
  chomp;

  my ($ab) = (split/\t/)[2];
  my (@ab) = split/,/, $ab;
  $par[$ab[0]][$ab[0]]++ if $#ab == 0; # one protein by this cluster
  for (my $i = 0; $i <= $#ab-1; $i++) {
    $par[$ab[$i]][$ab[$i]]++;
    for (my $j = $i+1; $j <= $#ab; $j++) {
      $par[$ab[$i]][$ab[$j]]++;
    }
  }
}
close in;

print "#";
for (sort {$a <=> $b} keys %S) { print "\tab$_" }
print "\n";

foreach my $i (sort {$a <=> $b} keys %S) {
  print "ab$i";
  foreach my $j (sort {$a <=> $b} keys %S) {
    $par[$i][$j] = 0 if !$par[$i][$j];
    print "\t$par[$i][$j]";
  }
  print "\n";
}

exit;
