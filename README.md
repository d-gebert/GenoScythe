# gRNA Heuristic Search Tool for CRISPR-Cas9

This repository contains a Perl script for generating and filtering guide RNA (gRNA) sequences for CRISPR-Cas9 experiments, tailored to user-specified genome and annotation files. The tool allows precise selection of gRNAs that meet specific criteria and ensures safe targeting by avoiding gene exons and regulatory regions.

---

## Features
- **gRNA Design:** Generates gRNAs with NGG PAM sequences for CRISPR-Cas9 applications.
- **Customizable Parameters:** Supports user-defined sequence length, mismatch tolerance, chromosome filtering, and more.
- **Off-Target Filtering:** Removes gRNAs with off-target hits in unwanted regions.
- **Integration with Bowtie:** Maps sequences to the genome for hit counting and filtering.
- **Safety Filtering:** Ensures gRNAs target transposable elements or introns only.

---

## Requirements

### Dependencies
- Perl 5 or later
- Modules:
  - `Getopt::Long`
  - `Digest::MD5`
- External tools:
  - `bowtie`
  - `samtools`

Install Perl modules via CPAN:
```bash
cpan install Getopt::Long Digest::MD5
```
Install Bowtie and Samtools via your package manager or from their official websites.

---

## Usage

### Required Inputs
- **Genome fasta file (`-g`)**: The genome sequence in FASTA format.
- **Gene set GTF file (`-f`)**: Genome annotations.
- **RepeatMasker file (`-r`)**: Repeat regions.

### Example Command
```bash
perl gRNA_search.pl -g genome.fa -f annotations.gtf -r repeats.out \
  -seq_len 20 -proximal_len 12 -distal_len 8 -n_distal_mms 2 \
  -max_srep 4 -n_seqs 100000 -max_n_hits 30 -n_sets 1
```

### Options
| Option                  | Description                                                                                  | Default        |
|-------------------------|----------------------------------------------------------------------------------------------|----------------|
| `-genome` or `-g`      | Genome sequence fasta file.                                                                  | Required       |
| `-geneset` or `-f`     | Genome gene set GTF file.                                                                    | Required       |
| `-repeats` or `-r`     | RepeatMasker output file.                                                                    | Required       |
| `-incl_chrs`           | Included chromosomes for hits (comma-separated).                                            | `all`          |
| `-seq_len`             | Total sequence length (excluding PAM).                                                      | 20             |
| `-proximal_len`        | Proximal part length (requires perfect match).                                              | 12             |
| `-distal_len`          | Distal part length (allows mismatches).                                                     | 8              |
| `-n_distal_mms`        | Max allowed mismatches in the distal part.                                                  | 2              |
| `-max_srep`            | Max number of simple repeats in sequence.                                                   | 4              |
| `-n_seqs`              | Number of initially generated sequences.                                                    | 100000         |
| `-max_n_hits`          | Max number of hits per gRNA in the output.                                                  | 30             |
| `-n_sets`              | Number of sets in the output.                                                               | 1              |
| `-filt_mm_seqs`        | Filter gRNAs with off-target hits (mismatches).                                             | No             |
| `-start_g`             | Restrict gRNAs to start with a G at the 5' end.                                             | No             |
| `-help` or `-h`        | Display help message.                                                                       | -              |

---

## Output
- **gRNA FASTA file (`gRNA_seqs.fas`)**: Contains unique gRNA sequences.
- **Summary table (`gRNA_seqs.fas.hits.tbl`)**: Summary of gRNA hits across chromosomes.
- **Detailed hits table (`gRNA_seqs.fas.hits_detail.tbl`)**: Detailed mapping of gRNA hits, including positions and strand.

---

## Notes
- Ensure Bowtie indices are built before running the script if they are not already present.
  ```bash
  bowtie-build genome.fa genome
  ```
- Outputs are overwritten if the script is run multiple times with the same parameters.

---

## License
This project is licensed under the MIT License.

---

## Author
Daniel Gebert  
dg572@cam.ac.uk  

For questions or issues, feel free to reach out or open a GitHub issue.

