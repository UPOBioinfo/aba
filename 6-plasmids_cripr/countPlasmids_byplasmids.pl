#!/usr/bin/perl
use strict;
use Statistics::Descriptive;

my $f1 = $ARGV[2];
my $f2 = $ARGV[3];
my %ab;
my @n;
my @str; my @str2;
my %nn; my %nn2;
my @ff; my @ff2;
my $npl = 0;
my $stat = Statistics::Descriptive::Full->new();
my $stat2 = Statistics::Descriptive::Full->new();
my @pl_genes;
my @pl_genes1;
my @pl_genes2;
my $Nstr;
my %C;
my %crispr;
my %crispr2;
my $nnn;
my $nnn2;
my %prokka;

# Gather crispr plasmids
open pl, "./plasmids_predicted_crispr.id";
while (<pl>) {
  chomp;
  $C{$_} = 1;
}
close pl;

# Gather ab identifiers
open ab, $ARGV[1];
while (<ab>) {
  chomp;
  $ab{$_} = 1; 
}
close ab;

# Gather ab identifiers
open ab, "./ab_crisprGenes_repeats_predicted.tsv";
while (<ab>) {
  chomp;
  my ($ab, $n1, $n2) = split/\t/;
  $prokka{$ab} = "$n1\t$n2"; 
}
close ab;

# Go through strains
open in, $ARGV[0];
my $pl = <in>;
chomp $pl;
my @pl = split /\t/, $pl;
while (<in>) {
  chomp;

  $Nstr++;
  my (@l) = split/\t/;
  for (my $x = 1; $x <= $#l; $x++) {
    $n[$x] += $l[$x];
  }

  #next unless $ab{$l[0]}; # only selected strains
  if ($ab{$l[0]}) { # selected strains
    push @str, \@l;
  } else {
    push @str2, \@l;
  }
}
close in;

# percent of strains
my $F1 = $f1 * $Nstr / 100;
my $F2 = $f2 * $Nstr / 100;

# Filter by plasmid frequency
my %abab = ();
for (my $x = 1; $x <= $#n; $x++) {
  #next unless $n[$x] > $F1 && $n[$x] <= $F2 && $n[$x] > 0; # greater than f1 and <= f2
  @{$abab{$pl[$x]}} =();
  $npl++;
  push @pl_genes, $pl[$x];
  foreach my $s (@str) {
    my @l = @{$s};  
    $nn{$pl[$x]} += $l[$x];
    $nnn += $l[$x];
    $crispr{$l[0]}++ if $C{$pl[$x]} == 1 && $l[$x] == 1;
    push @pl_genes1, $pl[$x] if $l[$x] == 1;
    push @{$abab{$pl[$x]}}, $l[0] if $l[$x] == 1;
  }
  foreach my $s (@str2) {
    my @l = @{$s};    
    $nn2{$pl[$x]} += $l[$x];
    $nnn2 += $l[$x];
    $crispr2{$l[0]}++ if $C{$pl[$x]} == 1 && $l[$x] == 1;
    push @pl_genes2, $pl[$x] if $l[$x] == 1;
    push @{$abab{$pl[$x]}}, $l[0] if $l[$x] == 1;
  }
} 

# Freq final
my @out;
foreach (keys %nn) {
  $crispr{$_} = 0 if !$crispr{$_};
  my $nnn = $nn{$_} - $nn2{$_};
  print "$_\t$nn{$_}\t$nn2{$_}\t$nnn";
  print "\t" if @{$abab{$_}};
  print join ";", @{$abab{$_}} if @{$abab{$_}};
  print "\n";
}

exit;
