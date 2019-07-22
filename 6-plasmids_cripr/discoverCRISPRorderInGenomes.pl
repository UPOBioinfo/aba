#!/usr/bin/perl
use strict;

# crispr genes
my %cri;
#open in, "./crispr_all_gn.ab";
#open in, "./crispr_all_gn_joined.ab";
open in, "./crispr_all_gn_joined_v3.ab";
while (<in>) {
    chomp;
      
    my ($ab, $gn) = split/\t/;
    $cri{$ab} = $gn;
} 
close in;

# pan_gen
my %pan;
open in, "./pan_gn.fasta";
while (<in>) {
  chomp;

  my (@ab) = split/\t/;
  $ab[0] =~ s/.+: //;
  next unless ($cri{$ab[0]});
  
  @{$pan{$ab[0]}} = @ab;
}
close in;

# go through crispr strains
#open in, "./crispr_all_strains.ab";
#open in, "./crispr_all_strains_v3.ab";
open in, "./strains_ab.ab";
while (<in>) {
  chomp;
  my $str = $_;

  my @clu; my %name;
  foreach my $gn (keys %pan) {
    foreach my $g ( grep /$str/, @{$pan{$gn}} ) {
      $name{$g} = $cri{$gn};
      push @clu, $g;
    }
  }

  # Print sorted
  my $prev = 0;
  my $nn = 0;
  my $out = "";
  foreach my $c (sort {$a cmp $b} @clu) {
    $c =~ m/ab[0-9]{5}_0+([0-9]+)/;
    my $current = $1;
    my $diff = $current - $prev;
    if ($diff > 1 && $prev != 0) {
      my @gap;
      if ($diff <= 4) {
        my $n;
        for (my $x = 1; $x < $diff; $x++) {
          $n = sprintf "%05d", $prev + $x;
          push @gap, "$str\_$n";
        }
      }
      $diff--;
      $out .= "\tX($diff";
      $out .= ":" if (@gap);
      $out .= join ";", @gap;
      $out .=  ")";
    }

    $out .= "\t";
    #if ($name{$c} eq "cas9") { my ($a1, $a2) = split/_/, $c; printf "%s\_%05d\-", $a1, $a2-1; }
    $out .= $name{$c};
    #if ($name{$c} eq "cas9") { my ($a1, $a2) = split/_/, $c; printf "-%s\_%05d", $a1, $a2+1; }
    $prev = $current;
    $nn++;
  }
  print "$str\t$nn$out\n";
}
close in;

