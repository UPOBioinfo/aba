#!/usr/bin/perl
use strict;
# Create a matrix of gene presence/abscence for R
# INPUT1: pan_clusters2_table.tsv from readPangenome.pl
# INPUT2: list of strains id in txt format
# OUTPUT1: matrix in tsv format for R scripts 
# AJPerez, 2019 (updated June 2021)

my %S;
my @par;

if (!$ARGV[0]) { die "Please run as: ./buildMatrix4heatmap.pl pan_clusters2_table.tsv strains.txt\n" }
my $FILE = $ARGV[0];
my ($a1, $a2) = split/\./, $FILE;

my $OUT_FILE = $a1 . ".matrix";

my $file_ab = $ARGV[1] || "./strains.txt";

my $PASS;
my $LETTERS = "ab";
open in, $file_ab;
while (<in>) {
  chomp;

  $_ =~ /([a-z]+)/;
  $LETTERS = $1 if !$PASS;
  $PASS = 1;
  $_ =~ s/$LETTERS//;
  $_ =~ s/^0+//;
  $S{$_} = $_;
}
close in;

open in, $FILE;
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

open OUT, ">$OUT_FILE";
print OUT "#";
for (sort {$a <=> $b} keys %S) { print "\t$LETTERS$_" }
print OUT "\n";

foreach my $i (sort {$a <=> $b} keys %S) {
  print OUT "$LETTERS$i";
  foreach my $j (sort {$a <=> $b} keys %S) {
    $par[$i][$j] = 0 if !$par[$i][$j];
    print OUT "\t$par[$i][$j]";
  }
  print OUT "\n";
}

print "Created file: $OUT_FILE\n";

exit;
