#!/usr/bin/perl

use HyGe;

my $k = $ARGV[0];
my $n = $ARGV[1];
my $K = $ARGV[2];
my $N = $ARGV[3];

print HyGe::ComputeHyperPValue($k, $n, $K, $N) . "\n";

exit;
