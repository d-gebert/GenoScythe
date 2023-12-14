# GenoScythe

The purpose of this program is the generation of guide RNAs (gRNAs) for
CRISPR-Cas9 experiments with a range of different numbers of target
sites in a given genome. For gRNA sequences with multiple target sites
the aim is to spread out target sites across chromosomes as much as 
possible, while only allowing target sites in transposable elements or
gene introns. This ensures that no gene exons or regulatory sequences
will be affected by double strand breaks induced by Cas9.
This script requires at least three parameters. The name of a genome
sequence fasta file, the name of a corresponding geneset GTF file, and 
a repeatmasker output file for the used genome.
