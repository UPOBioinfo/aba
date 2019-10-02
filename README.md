# ABA Project: Genomics for personalized medicine against *Acinetobacter baumannii*

These are both the computational tools and protocol used to annotate and analyse the pan genome of Acinetobacter baumannii. In this project 2,467 assembled genomes of this bacterium were used , but this protocol is potentially useful to study any other prokaryotic species.

## Genome annotation
FASTA files for all strains are downloaded from NCBI Genome database. Files were used by [Roary](https://sanger-pathogens.github.io/Roary/)  for calculating the core and pan genome, with both -s and -o parameters, and 90% as the identity threshold, which creates a pan_genome file. The reference protein-coding sequences from this file are functionally annotated using [Sma3s v2](http://www.bioinfocabd.upo.es/sma3s/) and the taxonomic division bacteria of UniProt Knowledgebase as annotation source.

Then, the script ***joinPanfastaByGenename.pl*** is used to join genes which has the same gene name.

Now, we can show both the distribution of the number of genes and the average number of shared genes, using ***script_shared_genes.R*** (or ***script_shared_genes_gn.R*** for the final group of strains). These scripts require that the pangenome file were converted by this two Perl scripts: ***readPangenome.pl***, ***buildMatrix4heatmap.pl***.

## Phylogenies
### 16S rRNA-based phylogeny
ssu-align and Infernal are used to find 16S gene in all the strains. Then, they are aligned with MAFFT and trimAL, and the phylogeny is constracted with DNAPARS (PHYLIP). The final tree is shown with R ggtree, and Average Nucleotide Identity (ANI) values are added.

### Core phylogeny (phylogenomics)
Core genes are joined and aligned with MAFFT and trimAL, and the tree is also shown with ggtree. Genome metadata can be added when they are download from NCBI BioSample database.

### Phylogenetic profile (presence/absence matrix)
It is obtained from the Roary output, and it can be shown by Phandango web server.

All phylogenies can be shown by ggtree, and we use two different R scripts for it: ***script_phylo_16S.R*** and ***script_phylo_core***.

## Functional enrichment
When you have a set of gene identifiers (ab00000), you can perform a functional enrichment to find highlighter functions. To do that, you can use the script ***script_core_enrichment.R***, which uses the TopGO R library.

TopGO use GO term, but if you want also to use all the other annotation that give Sma3s, including UniProt Keywords, you can use the following Perl script, which perform a funcional enrichment using the hypergeometric distribution: ***calculateEnrichmentFromSma3s.pl*** (it depends on ***calculatePvalueForAnnotation.pl***).

In this case, to compare the enrichment from two different group you can use ***script_plasmid_enrichment.R***.

To compare two sets of strains, and check the frequency in each set, you can use ***calculatePercentGenesInStrains.pl***.

## Search for plasmids
Plasmid sequences can be downloades from NCBI Nucleotide database, with descriptions "complete genome" or "complete sequence". Then you can search in the genomic sequences by ***checkPlasmids.pl***. It assumes both query and subject coverage of 95%. Then, plasmids can be clustered by psi-CD-HIT to remove redundancy, in addition to use mob_cluster:
>mob_cluster --mode build --infile plasmids_acinetoplasmids_baumannii.fasta --outdir mob_plasmids

The matrix of plasmid presence/absence can be constructed by:
>constructPlasmidMatrix_withMOB.pl > ab_plasm_matrix_2112_130.tsv

And heatmaps can be created using the following R scripts: ***plasmid_matrix.R***, ***plasmid_matrix_ta.R***.

## Search for CRISPR repeats and Cas genes
CRISPR arrays can be found from plasmid and genome sequences, using its GenBank entries:
>grep -E "VERSION|CRISPR Arrays" plasmids_predicted.gb | grep "CRISPR Arrays" -B1 | grep VERSION | cut -c13- > plasmids_predicted_crispr.id

But CRISPRCasFinder for a complete prediction with default parameters. Then, you can collect the results:
>createSummaryFromCrispFinderResults.pl > ab_crisprGenes_repeats_predicted.tsv

To find cas genes, you can search proteins in UniProtKB database by this GO term: GO:0043571; maintenance of CRISPR repeat elements. And then, search for homologs and cut by a threshold with: ***filterBlastp.pl*** and ***searchCRISPR.pl***.

You can also search from the Sma3s annotation or even from the Prokka annotation:
>./searchCRISPRfromGFFProkka.pl > ab_crisprGenes_repeats_abs.tsv

Finally, you can search cas genes in a highly precise mode with RPS-BLAST and PSSM matrices from the CDD database:
>ls *.smp > cdd_crispr.pn
>makeprofiledb -in cdd_crispr.pn -title cdd_crispr -dbtype 'rps'
>rpsblast+ -query pan_gn.faa -db cdd_crispr.pn -evalue 1e-05 -outfmt '6 qseqid sseqid pident qcovs qcovhsp length qlen slen evalue qstart qend sstart send' > pan_gn.rpsblast.tsv

Once you collect all cas genes, you can find CRISPR/Cas cluster in the genomes, using Prokka annotations and ***discoverCRISPRorderInGenomes.pl***.

With all the results, you can create summary tables which can be used to analyse the plasmid distributions, in addition to the CRISPR arrays and cas genes: ***script_plasmids_crispr.R***.

Tables need by the R script can be created by:
> countPlasmids.pl ab_plasm_matrix_2112_84.tsv crispr_all_strains_v3.ab 0 100 > plasmids_crispr_0_100_v3.tsv
> countPlasmids.pl ab_plasm_matrix_2112_84.tsv cas9_v3_strains.ab 0 100 > plasmids_crispr_0_100_cas9_v3.tsv
> countPlasmids.pl ab_plasm_matrix_2112_84.tsv group1.tsv 0 100 > plasmids_group1_0_100_v3.tsv

The variant ***countPlasmids_byplasmids.pl*** can be used to obtains the plasmids linked to each group.

## Software versions
Prokka version 1.13
Sma3s v2
BLAST 2.2.31+ (RPS-BLAST)
ssu-align v0.1.1
Infernal v1.1
MAFFT v7.271
trimAL v1.2
PHYLIP v3.697
RAxML v8.2.9
R ggtree library v1.10.5
pyani v0.2.4
Roary version 3.11.2
FastTree 2.1.8
Phandango
TopGO R package version 2.30.1
CRISPRCasFinder 1.4
psi-CD-HIT
Mob-Suite toolkit v1.4.9.1
R libraries: cowplot, dendextend, dplyr, ggbeeswarm, ggnewscale, ggplot2, ggpubr, ggrepel, ggtree, magrittr, pheatmap, plotly, RColorBrewer, reshape2, tidyverse, topGO, viridis, viridisLite
