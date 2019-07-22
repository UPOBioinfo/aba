#!/usr/bin/perl

# AJPerez, 2007-05-03
# Calculate the frequencies for annotations from gene set
# params: reference_set experiment_set

use lib './';
use HyGe; # Hypergeometric distribution lib

my $file_exp = $ARGV[0];
my $file_ref = $ARGV[1];

my (%goes) = &go_descriptions;

# Experiment file
my $total_exp = 0;
my %id_exp;
open (in, $file_exp) || die "Error! Problem opening $file_exp\n";
while (<in>) {
  chomp $_;
  $total_exp++;

  $id_exp{$_} = 1;
}
close in;

# Run through the reference
my $total_ref = 0;
open (in, $file_ref) || die "Error! Problem opening $file_ref\n";
while (<in>) {
  chomp $_;
  $total_ref++;
  
  my ($a1, $id, @annots, $a2) = split /\t/, $_;
  
  # frequency for each annotation
  my $index = 0;
  foreach my $list (@annots) {
    $index++;
    foreach my $annot (split /\;/, $list) {
      $annot_ref[$index]{$annot}++;
      #print $annot[$index]{$annot} . "\t$index\t$annot\n"; 
    
      next unless ($id_exp{$id});
      $annot_exp[$index]{$annot}++;
      $ids_exp[$index]{$annot} .= "$id,";
    }
  }

}
close in;

# Final annotations
my (@output) = "Annotation\tdescription\tquery_freq\tref_freq\tquery_genes\tref_genes\tHyGe_p-value\tgenes";
for (my $x=1; $x<=$#annot_exp; $x++) {
  foreach (keys %{$annot_exp[$x]}) {
    chop $ids_exp[$x]{$_}; # last comma
    my $link = &includeLink($_);
    #my $link = $_;
    
    my $Hg = sprintf("%.5g", HyGe::ComputeHyperPValue($annot_exp[$x]{$_}, $total_exp, $annot_ref[$x]{$_}, $total_ref));
    #$Hg =~ s/\./\,/;
    
    push (@output, "$link\t$goes{$_}\t$annot_exp[$x]{$_}\t$annot_ref[$x]{$_}\t$total_exp\t$total_ref\t" 
      . $Hg . "\t"
      . $ids_exp[$x]{$_} );
  }
}
print join "\n", map {$_->[0]} sort {$a->[7] <=> $b->[7]} map {[$_,split(/\t/)]} @output;
print "\n";

exit;

# SUBRUTINES
############
sub go_descriptions () {
  my %goes;
  
  # GO
  open (go, "./GO.terms_and_ids");
  while (<go>) {
    chomp;
    next unless (/^GO:/);

    my ($id, $de, $type) = split /\t/;
    $goes{$id} = "($type) $de";
  }
  close go;

  # InterPro
  open (ipro, "./interpro.de");
  while (<ipro>) {
    chomp; chop;
    my ($id, $de) = split /\; /;
    $goes{$id} = $de;
  }
  
  return %goes;
}

###
sub includeLink () {
  my ($id) = @_;

  my %link = (
    ipro => "http://www.ebi.ac.uk/interpro/IEntry?ac=",
    kw   => "http://www.expasy.org/cgi-bin/get-entries?KW=",
    go   => "http://www.ebi.ac.uk/ego/DisplayGoTerm?id=",
    itc  => "http://www.ebi.ac.uk/intact/search/do/search?searchString="
  );

  $id =~ s/(IPR\d+)/=hyperlink("$link{ipro}$1";"$1")/;
  $id =~ s/(GO:\d+)/=hyperlink("$link{go}$1";"$1")/;
  $id =~ s/(EBI-\d+)/=hyperlink("$link{itc}$1";"$1")/;

  return $id;
}
