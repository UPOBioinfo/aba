#!/usr/bin/perl
use strict;
# Create a TSV with 4 columns: reference gene name, number of strains with gene, list of strains with and without gene
# INPUT1: pancluster from roary in tsv format (preferently, the collapsed, 2, version)
# INPUT2: list of strains id in txt format
# OUTPUT1: TSV with reference gene, number of strains with this gene, and list of strains with and without the gene 
# AJPerez, 2019 (updated June 2021)

my @CL;
my %CL;
my @S;

if (!$ARGV[0]) { die "Please run as: ./readPangenome.pl pan_clusters2.tsv strains.txt\n" }
my $FILE = $ARGV[0];
my ($a1, $a2) = split/\./, $FILE;
my $OUT_FILE = $a1 . "_table.tsv";

my $file_ab = $ARGV[1] || "./strains.txt";

open in, $file_ab;
while (<in>) {
  chomp;

  $_ =~ s/[a-z]+//;
  $_ =~ s/^0+//;
  $S[$_] = 1;
}
close in;

open FILE, $FILE;
while (<FILE>) {
  chomp;

  my ($g, @cl) = split/\t/;
  my ($idc, $ref) = split/: /, $g;
  foreach my $c ($ref, @cl) {
    ($c) = split/_/, $c;
    $c =~ s/[a-z]+//;
    $c =~ s/^0+//;
    push @{$CL{$idc}}, $c;
  }
}
close FILE;

# Read clusters
open OUT, ">$OUT_FILE";
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
  print OUT "Cluster $c\t$n\t";
  print OUT join ",", @{$cl};
  print OUT "\t";
  print OUT join ",", @not;
  print OUT "\n";
}
close OUT;

print "Created file: $OUT_FILE\n";

exit;
