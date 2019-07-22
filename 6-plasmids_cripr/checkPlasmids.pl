#!/usr/bin/perl
use strict;

my $s = $ARGV[0];
my $db = $ARGV[1] || "plasmids_ab.fasta";
my $pid = $ARGV[2] || 95;
my %n;
my %p;
my %pl;
my %c;

my %org;
open in, $db;
while (<in>) {
  chomp;

  next unless /^>/;
  my ($id, @org) = split/ /;
  $id =~ s/>//;
  $org{$id} = join " ", @org;
}
close in;

my @blast = `blastn -num_threads 4 -query $s -db $db -perc_identity $pid -outfmt '6 qseqid sseqid sstart send qlen slen qcovs length pident'`;
foreach (@blast) {
  chomp;
  my ($qid, $sid, $p1, $p2, $qlen, $slen, $qcovs, $length, $pident) = split /\t/;

  next unless $slen >= 1000; # plasmids >= 1000nt
  ($p1, $p2) = ($p2, $p1) if $p1 > $p2;
  $pl{$sid} = $slen;

  next unless $qcovs > 95; # Query cov threshold
  push @{$c{$sid}}, "$qid($qlen;$qcovs)";

  for (my $x = $p1; $x <= $p2; $x++) {
    $n{$sid}++ if !$p{$sid}{$x};
    $p{$sid}{$x} = 1;
  }
}

foreach my $pl (keys %pl) {
  my $per = sprintf("%.2f", ($n{$pl} / $pl{$pl}) * 100);
  $per =~ s/\./,/;

  next unless $per > 95;
  #my ($s1, $s2) = split/\|/, $pl;
  #$s2 =~ s/_/ /;
  printf "$s\t$pl\t$org{$pl}\t$pl{$pl}\t$per\t";
  my %saw; @saw{@{$c{$pl}}} = (); @{$c{$pl}} = keys %saw;
  print join ",", @{$c{$pl}};
  print "\n";
}

exit;
