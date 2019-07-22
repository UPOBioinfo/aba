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
#open pl, "./plasmids_predicted_48_ta.id";
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
#open ab, "./ab_crisprGenes_repeats.tsv";
open ab, "./ab_crisprGenes_repeats_predicted_v3.tsv";
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
for (my $x = 1; $x <= $#n; $x++) {
  next unless $n[$x] > $F1 && $n[$x] <= $F2 && $n[$x] > 0; # greater than f1 and <= f2
  $npl++;
  push @pl_genes, $pl[$x];
  foreach my $s (@str) {
    my @l = @{$s};  
    $nn{$l[0]} += $l[$x];
    $nnn += $l[$x];
    $crispr{$l[0]}++ if $C{$pl[$x]} == 1 && $l[$x] == 1;
    push @pl_genes1, $pl[$x] if $l[$x] == 1;
  }
  foreach my $s (@str2) {
    my @l = @{$s};    
    $nn2{$l[0]} += $l[$x];
    $nnn2 += $l[$x];
    $crispr2{$l[0]}++ if $C{$pl[$x]} == 1 && $l[$x] == 1;
    push @pl_genes2, $pl[$x] if $l[$x] == 1;
  }
} 

# Freq final
my @out;
foreach (keys %nn) {
  $crispr{$_} = 0 if !$crispr{$_};
  push @out, "group 1\t$f1-$f2%\t$crispr{$_}\t$_\t$nn{$_}\t$nnn\t$prokka{$_}";
  push @ff, $nn{$_};
}
foreach (keys %nn2) {
  $crispr2{$_} = 0 if !$crispr2{$_};
  push @out, "group 2\t$f1-$f2%\t$crispr2{$_}\t$_\t$nn2{$_}\t$nnn2\t$prokka{$_}";
  push @ff2, $nn2{$_};
}
$stat->add_data(@ff);
$stat2->add_data(@ff2);
my $means = sprintf "%.02f\t%.03f", $stat->mean(), $stat->standard_deviation();
my $means2 = sprintf "%.02f\t%.03f", $stat2->mean(), $stat2->standard_deviation();
print "set\tpercent\tcrispr\tstrain\tnplasmids\ttotal_plasmids\tprokka_crispr\tprokka_arrays\tmean\tstdev\n";
foreach (@out) {
  my ($gr) = split/\t/;
  print "$_\t";
  if ($gr eq "group 1" ) { 
    print "$means\n";
  } else {
    print "$means2\n";
  }
}

exit;
