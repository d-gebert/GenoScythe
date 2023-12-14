#!/usr/bin/perl
# Heuristic search for gRNA sequences (CRISPR-Cas9) with different numbers of genomic targets (+NGG PAM)

use strict;
use warnings;
# Load modules
use Getopt::Long;
use Digest::MD5 qw(md5_hex);

# Options
my $genome_file = '';               # Genome sequence fasta file
my $gnmgtf_file = '';               # Genome gene set gtf file
my $repmsk_file = '';               # Genome repeatmasker outfile
my $incl_chrs = 'all';              # Included chromosomes (hits only allowed on these)
my $seq_len = 20;                   # Total sequence length (excluding PAM)
my $proximal_len = 12;              # Sequence length of proximal part (requiring perfect match)
my $distal_len = 8;                 # Sequence length of distal part (allowing n mismatches)
my $n_distal_mms = 2;               # Max allowed mismatches in distal part
my $max_srep = 4;                   # Max number of simple repeats in sequences
my $n_seqs = 100000;                # Number of initially generated sequences
my $max_n_hits = 30;                # Max number of hits for final output
my $n_sets = 1;                     # Number of sequence sets to be created in final output
my $filt_mm_seqs = 0;               # Filter sequences that have (off-target) hits with mismatches
my $start_g = 0;                    # Sequence has to start with a G at the 5' end
my $help = 0;                       # Print help message
# Program names/paths (dependencies)
my $bowtie = 'bowtie';
# Global constants
my $id_len = 10;
# Global variables
my @argv = @ARGV;
$|=1; #Autoflush

# Program name
print("\n--- $0 ---\n");
# Test dependencies
my $depends = '';
foreach my $prog ($bowtie) {
	system("$prog -h 2> .progtest_err > .progtest_log");
	open(my $fh, '<', ".progtest_err");
	my $frst_ln = <$fh>;
	if ($frst_ln && $frst_ln =~ /not found/) { $depends .= "$prog\n" }
	unlink('.progtest_err','.progtest_log');
}
if ($depends) { die("\nERROR: Cannot start program.\nFollowing dependencies need to be installed:\n$depends\n"); }

# Collect command line arguments
GetOptions (
	"genome|g=s"	    => \$genome_file,	# Genome sequence fasta file
	"geneset|f=s"	    => \$gnmgtf_file,	# Genome gene set gtf file
	"repeats|r=s"	    => \$repmsk_file,	# Genome repeatmasker outfile
    "incl_chrs=s"	    => \$incl_chrs,		# Included chromosomes (hits only allowed on these)
    "seq_len=i"	        => \$seq_len,		# Total sequence length (excluding PAM)
    "proximal_len=i"	=> \$proximal_len,	# Sequence length of proximal part (requiring perfect match)
    "distal_len=i"	    => \$distal_len,	# Sequence length of distal part (allowing n mismatches)
    "n_distal_mms=i"	=> \$n_distal_mms,	# Max allowed mismatches in distal part
    "max_srep=i"	    => \$max_srep,      # Max number of simple repeats in sequences
    "n_seqs=i"	        => \$n_seqs,        # Number of initially generated sequences
    "max_n_hits=i"	    => \$max_n_hits,    # Max number of hits for final output
    "n_sets=i"	        => \$n_sets,        # Number of sequence sets to be created in final output
    "filt_mm_seqs"	    => \$filt_mm_seqs,  # Filter sequences that have (off-target) hits with mismatches
    "start_g"	        => \$start_g,       # Sequence has to start with a G at the 5' end
	"help|h"		    => \$help)			# Help
or die(usage());

## Handle options / error messaging
# Help message
if ($help) { help() };
# Usage message if run without options
if (!@argv) { die(usage()) }
# Input files
if (!$genome_file || !$gnmgtf_file || !$repmsk_file) { die(usage("\nError! Input files not specified.")) }
unless (-e $genome_file) { die(usage("\nError! File \'$genome_file\' does not exist.")) }
unless (-e $gnmgtf_file) { die(usage("\nError! File \'$gnmgtf_file\' does not exist.")) }
unless (-e $repmsk_file) { die(usage("\nError! File \'$repmsk_file\' does not exist.")) }

# Exaggerated program title
print("
___________________________________________________________

  ###                     ###               #    #          
 #      ##   # #    ##   #      ###  #  #  ####  ###    ##  
 # ##  ####  ## #  #  #   ##   #     #  #   #    #  #  #### 
 #  #  #     #  #  #  #     #  #      ###   #    #  #  #    
  ###   ##   #  #   ##   ###    ###     #    ##  #  #   ##  
                                      ##                    
___________________________________________________________

");

# Process options
# Make sure proximal + distal lengths = total length
if ($proximal_len + $distal_len != $seq_len) {
    die("ERROR! Proximal and distal lengths have to be equal to total sequence length!\n\n");
}
# Included chromosomes
my @incl_chrs = split(',', $incl_chrs);
if ($incl_chrs eq 'all') {
    my $all_chrs = get_fasta_headers($genome_file,1);
    @incl_chrs = @{$all_chrs};
}
my %incl_chrs = map { $_ => 1 } @incl_chrs;
# Get bools for parameters output
my $filt_mm_seqs_bool = $filt_mm_seqs ? 'yes' : 'no';
my $start_g_bool = $start_g ? 'yes' : 'no';

print("Extracting genome...");
# Extract genomic sequences
my $chr_seqs = get_fasta_seqs($genome_file,1);
# Get genome length (only included chromosomes)
my $gnm_len = 0;
foreach my $chr (sort keys %{$chr_seqs}) {
    next unless $incl_chrs{$chr};
    $gnm_len += length($chr_seqs->{$chr});
}
print("ok\n");

print("Extracting GTF data...");
# Extract intron positions from gtf file
my $int_poss = get_intron_positions($gnmgtf_file);
print("ok\n");

print("Extracting repeat data...");
# Extract repeat positions from repeatmasker file
my $rep_poss = get_repeat_positions($repmsk_file);
print("ok\n");

print("Generating random sequences with parameters:\n");
print(
"  Included chromosomes (only hits on these) ...... $incl_chrs\n",
"  Sequence length (excluding PAM) ................ $seq_len\n",
"  Length of proximal part of sequences ........... $proximal_len\n",
"  Length of distal part of sequences ............. $distal_len\n",
"  Max. mismatches in distal part ................. $n_distal_mms\n",
"  Max. number of simple repeats in sequences ..... $max_srep\n",
"  Number of initially generated sequences ........ $n_seqs\n",
"  Max. number of hits for final output ........... $max_n_hits\n",
"  Number of sequence sets in final output ........ $n_sets\n",
"  Filter sequences having hits w/ mismatches ..... $filt_mm_seqs_bool\n",
"  All sequences have to start with a 5'-G ........ $start_g_bool\n",
);
# Seed the random number generator
# time|$$ combines the current time with the current process id
srand(time|$$);

# Get random sequences
my @rand_seqs = ();
# Go through each chromosome
foreach my $chr (sort keys %{$chr_seqs}) {
    # Skip non-specified chromosome
    next unless $incl_chrs{$chr};
    # Get chromosome length
    my $chr_len = length($chr_seqs->{$chr});
    # Skip short contigs
    next if $chr_len < $gnm_len/1000;
    # Determine proportionate number of sequences to be generated
    my $n_seqs_chr = int($n_seqs*($chr_len/$gnm_len));
    # Keep count of iterations
    my $j = 0;
    # Create number of random sequences
    for (my $i = 0; $i < $n_seqs_chr; $i++) {
        # Increment iterations count
        $j++;
        # Leave this chromosome/contig after too many iterations
        last if $j > $n_seqs_chr*100000;
        # Get a random position and the sequence at that position
        my $rand_pos = int(rand($chr_len-$seq_len-1)) + 1;
        my $rpos_end = $rand_pos+$seq_len-1;
        # Sequence position has to be located in intron or repeat
        my $in_int_rep = 0;
        if ($int_poss->{$chr}->{$rand_pos} && $int_poss->{$chr}->{$rpos_end}) {
            $in_int_rep = 1;
        } elsif ($rep_poss->{$chr}->{$rand_pos} && $rep_poss->{$chr}->{$rpos_end}) {
            $in_int_rep = 1;
        }
        unless ($in_int_rep) {
            $i--;
            next;
        }
        # Get sequence at determined random position
        my $rand_seq = substr($chr_seqs->{$chr},($rand_pos-1),$seq_len);
        # Get three bases immediately following the sequence
        my $pam_seq = substr($chr_seqs->{$chr},($rand_pos-1+$seq_len),3);
        # Sequence has to be followed by an NGG PAM
        unless ($pam_seq =~ /GG$/) {
            $i--;
            next;
        }
        # If option 'start sequence with G' is true
        if ($start_g) {
            # Sequence has to have a G at the 5' end
            unless ($rand_seq =~ /^G/) {
                $i--;
                next;
            }
        }
        # Sequence cannot have stretches of simple repeats (e.g. GGGGG or GTGTGTGTGT)
        if ($rand_seq =~ /(.)\1{$max_srep,}/ || $rand_seq =~ /(..)\1{$max_srep,}/) {
            $i--;
            next;
        }
        # Save sequence
        push(@rand_seqs,$rand_seq);
    }
}

# Remove redundant sequences
@rand_seqs = remove_redundant_elements(@rand_seqs);

# Open fasta output file
my $fasfile = find_unique_filename("gRNA_seqs.fas");
my $fas = open_outfile($fasfile);

# IDs hash (for catching collisions)
my %hash_ids = ();
# Save sequence ids
my %rand_seqs = ();
# Go throuch each randomly generated sequence
foreach my $rand_seq (@rand_seqs) {
    # Unique md5 id
    my $id = generate_short_md5_id($rand_seq, $id_len);
    $hash_ids{$id}++;
    # Handle a collision
    $id .= ".$hash_ids{$id}" if $hash_ids{$id} > 1;
    # Save sequence for current id
    $rand_seqs{$id} = $rand_seq;
    # Print sequence to output file
    print($fas ">$id\n$rand_seq\n");
}
close($fas);
print("ok\n");

print("Mapping sequences to genome with $bowtie:\n");
# Build a bowtie index if needed
my $genome_file_idx = $genome_file;
$genome_file_idx =~ s/\.fa.*$//;
unless (-e "$genome_file_idx.rev.1.ebwt") {
    system("$bowtie-build $genome_file $genome_file_idx");
}

# Map random sequences to genome with bowtie
my @cmd = (
    "$bowtie -a",
    "-n $n_distal_mms",
    "-l $distal_len",
    "-x $genome_file_idx",
    "-f $fasfile",
    "-S | samtools view -S -F 4",
    "> $fasfile.$genome_file_idx.sam"
);
system(join(' ', @cmd));
print("ok\n");

print("Filtering sequences according to hits...");
# Extract data from sam file
my $sam_data = get_tab_fields("$fasfile.$genome_file_idx.sam");

# Get sequence hits
my %seq_hits = ();
my %seq_hits_per_chr = ();
my %seq_hit_details = ();
# Parse sam file and get hit counts
foreach my $line (sort keys %{$sam_data}) {
    # Get sequence id and chromosome name
    my $seq_id = $sam_data->{$line}->[0];
    my $chr_id = $sam_data->{$line}->[2];
    # Get bitwise flag and hit position
    my $bwflag = $sam_data->{$line}->[1];
    my $hitpos = $sam_data->{$line}->[3];
    # Get mismatch positions
    my $mmposs = $sam_data->{$line}->[12];
    $mmposs =~ s/MD:Z://;
    # Check perfect matches in proximal part
    my $prox_matches = 0;
    # Plus strand hit
    if ($bwflag == 0) {
        # Number of proximal matches
        ($prox_matches) = ($mmposs =~ /(\d+)$/);
    }
    # Minus strand hit
    elsif ($bwflag == 16) {
        # Number of proximal matches
        ($prox_matches) = ($mmposs =~ /^(\d+)/);
    }
    # Check presence of PAM sequence
    my $pam_seq = '';
    # Plus strand hit
    if ($bwflag == 0) {
        # Get three bases immediately following the sequence (PAM)
        $pam_seq = substr($chr_seqs->{$chr_id},($hitpos-1+$seq_len),3);
    }
    # Minus strand hit
    elsif ($bwflag == 16) {
        # Get three bases immediately following the sequence (PAM)
        $pam_seq = substr($chr_seqs->{$chr_id},($hitpos-1-3),3);
        $pam_seq = rev_com($pam_seq);
    }
    # Save hit count for sequence if PAM present and perfect proximal match
    if ($pam_seq =~ /GG$/ && $prox_matches >= $proximal_len) {
        $seq_hits{$seq_id}++;
        $seq_hits_per_chr{$seq_id}{$chr_id}++;
        # Save hit details
        my $hit_end = $hitpos+$seq_len-1;
        my $hit_str = $bwflag == 0 ? '+' : '-';
        push(@{$seq_hit_details{$seq_id}},[$chr_id,$hitpos,$hit_end,$hit_str,$mmposs]);
    }
}

# Remove sequences with hits outside of the specified chromosomes
foreach my $id (keys %seq_hits) {
    # Get genomic hit count
    my $n_hits_gnm = $seq_hits{$id};
    # Get sum of hit counts for specified chromosomes
    my $n_hits_incl_chrs = 0;
    # Go through each specified chromosome
    foreach my $chr (keys %incl_chrs) {
        # Add up hit counts for specified chromosomes only
        $n_hits_incl_chrs += $seq_hits_per_chr{$id}{$chr} if $seq_hits_per_chr{$id}{$chr};
    }
    # Delete sequence if the hit counts (total vs spec chrs) are not identical
    if ($n_hits_gnm != $n_hits_incl_chrs) {
        delete $seq_hits{$id};
    }
}

# Get hits outside of introns and repeats
my %non_intrep_hits = ();
# Parse sam file and get hits outside of introns and repeats
foreach my $line (sort keys %{$sam_data}) {
    # Get sequence id and chromosome name
    my $seq_id = $sam_data->{$line}->[0];
    my $chr_id = $sam_data->{$line}->[2];
    # Get hit coordinates
    my $hit_pos = $sam_data->{$line}->[3];
    my $pos_end = $hit_pos+$seq_len-1;
    # Save non intron / repeat hit count for this sequence id
    my $in_int_rep = 0;
    if ($int_poss->{$chr_id}->{$hit_pos} && $int_poss->{$chr_id}->{$pos_end}) {
        $in_int_rep = 1;
    } elsif ($rep_poss->{$chr_id}->{$hit_pos} && $rep_poss->{$chr_id}->{$pos_end}) {
        $in_int_rep = 1;
    }
    unless ($in_int_rep) {
        $non_intrep_hits{$seq_id}++;
    }
}

# Remove sequences with hits outside of introns and repeats
foreach my $id (keys %seq_hits) {
    # Delete sequence if it has any hits outside of introns and repeats
    if ($non_intrep_hits{$id}) {
        delete $seq_hits{$id};
    }
}

# Optional filtering of sequences that have hits with mismatches
if ($filt_mm_seqs) {
    # Get hits with mismatches
    my %mm_hits = ();
    # Parse sam file and get hits with mismatches
    foreach my $line (sort keys %{$sam_data}) {
        # Get sequence id
        my $seq_id = $sam_data->{$line}->[0];
        # Get number of mismatches (edit distance)
        my $ed_dis = $sam_data->{$line}->[13];
        # Get number of hits with mismatches
        if ($ed_dis ne "NM:i:0") {
            $mm_hits{$seq_id}++;
        }
    }
    # Remove sequences with mismatch hits
    foreach my $id (keys %seq_hits) {
        # Delete sequences with mismatch hits
        if ($mm_hits{$id}) {
            delete $seq_hits{$id};
        }
    }
}

# Get number of hit chromosomes per sequence
my %hit_chrs_n = ();
# Go through each sequence
foreach my $id (keys %seq_hits) {
    # Save number of hit chromosomes
    my $hit_chrs_n = keys %{$seq_hits_per_chr{$id}};
    $hit_chrs_n{$id} = $hit_chrs_n;
}
print("ok\n");

print("Creating final output...");
# Open final hits output file
my $outfile = "$fasfile.hits.tbl";
my $out = open_outfile($outfile);
my $outfile1 = "$fasfile.hits_detail.tbl";
my $out1 = open_outfile($outfile1);

print($out1 "#seq_id\tchr\tstart\tend\tstrand\tmatch\tguide_RNA_sense\tguide_RNA_antisense\ttarget_with_pam\n");

for (my $i = 0; $i < $n_sets; $i++) {
    # Set number
    my $set_i = $i+1;
    # Output title line
    print($out "#set$set_i\n#seq_id\tguide_RNA_sense\tguide_RNA_antisense\ttotal_hits\thits_per_chr\n");
    # Process and output one sequence for each hit count
    my %done_hits = ();
    # Go through each sequence id, sorted by number of hits (ascending order) + number of hit chromosomes
    foreach my $id (sort {$seq_hits{$a} <=> $seq_hits{$b} || $hit_chrs_n{$b} <=> $hit_chrs_n{$a}} keys %seq_hits) {
        # Skip if a sequence with a particular hit count was already saved
        next if $done_hits{$seq_hits{$id}};
        # Leave loop if max hit count is reached
        last if $seq_hits{$id} > $max_n_hits;
        # Get reverse complement of target sequence
        my $grna_seq = rev_com($rand_seqs{$id});
        # Output sequence and hit count
        print($out "$id\t$rand_seqs{$id}\t$grna_seq\t$seq_hits{$id}\t");
        # Output hits per chromosome
        foreach my $chr (sort keys %{$seq_hits_per_chr{$id}}) {
            print($out "$chr:$seq_hits_per_chr{$id}{$chr}; ");
        }
        print($out "\n");
        # Output hit details
        foreach my $hit (@{$seq_hit_details{$id}}) {
            my $tar_seq = '';
            $tar_seq = substr($chr_seqs->{$hit->[0]},($hit->[1]-1),$seq_len+3) if $hit->[3] eq '+';
            $tar_seq = substr($chr_seqs->{$hit->[0]},($hit->[1]-1-3),$seq_len+3) if $hit->[3] eq '-';
            print($out1 "$id\t$hit->[0]\t$hit->[1]\t$hit->[2]\t$hit->[3]\t");
            print($out1 "$hit->[4]\t$rand_seqs{$id}\t$grna_seq\t$tar_seq\n");
        }
        print($out1 "\n");

        # Mark hit count as already saved
        $done_hits{$seq_hits{$id}} = 1;
        # Delete already printed sequence id
        delete $seq_hits{$id};
    }
    print($out "\n");
}
# Close output file
close($out);
close($out1);
print("ok\n");

print("\nDone! Exiting.\n");

exit;

################################################################################
# Subroutines
################################################################################

sub get_tab_fields {
	# Take name of tab file
	my($infile,$skip_header) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	if ($skip_header) { shift(@in_data); }
	# Global tree hashes for each species
	my %data_fields = ();
	# Set 0 as start index
	my $id_i = 0;
	# Go through file data
	foreach my $line (@in_data) {
		# Get line data
        my @d = split(/\t/,$line);
        # Save data fields
        @{$data_fields{$id_i}} = @d;
		$id_i++;
    }
	# Return data fields
	return \%data_fields;
}

sub remove_redundant_elements {
    # Take array
	my(@array) = @_;
	# Remove redundant elements by creating a hash
	my %hash = map { $_, 1 } @array;
	# Create array without redundancy
	my @unique_array = keys %hash;
    # Return filtered array
	return @unique_array;
}

sub get_exon_positions {
	# Take name of tab file
	my($infile) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Storage variable
	my %intron_positions = ();
    my %exon_positions = ();
    # Gene start flag
    my $start_flag = 0;
    # Intron start and end
    my $intron_beg = 0;
    my $intron_end = 0;
    #my $genes_i = 0;
	# Go through file data
	foreach my $line (@in_data) {
		# Get line data
        my @d = split(/\t/,$line);
        my $chr = $d[0];
        my $typ = $d[2];
        my $beg = $d[3];
        my $end = $d[4];
        # Start of a new gene / transcript
        if ($typ eq 'gene' || $typ eq 'mRNA') {
            $start_flag = 1;
            #$genes_i++ if $typ eq 'gene';
            #if ($genes_i > 1) { exit; }
        }
        # First exon of current gene
        if ($start_flag && $typ eq 'exon') {
            # Start of first intron
            $intron_beg = $end+1;
            # Reset gene start flag
            $start_flag = 0;
        }
        # Following exon of current gene
        elsif ($typ eq 'exon') {
            # End of current intron
            $intron_end = $beg-1;
            #print("$chr\t$intron_beg\t$intron_end\n");
            # Save intron positions
            foreach my $pos ($intron_beg..$intron_end) {
                $intron_positions{$chr}{$pos} = 1;
            }
            # Start of next intron
            $intron_beg = $end+1;
        }
        # Save all exon positions
        if ($typ eq 'exon') {
            foreach my $pos ($beg..$end) {
                $exon_positions{$chr}{$pos} = 1;
            }
        }
    }
    # Remove exonic positions from intron positions
    foreach my $chr (sort keys %intron_positions) {
        foreach my $pos (sort keys %{$intron_positions{$chr}}) {
            if ($exon_positions{$chr}{$pos}) {
                delete($intron_positions{$chr}{$pos});
            }
        }
    }
	# Return data fields
	return \%exon_positions;
}

sub get_intron_positions {
	# Take name of tab file
	my($infile) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Storage variable
	my %intron_positions = ();
    my %exon_positions = ();
    # Gene start flag
    my $start_flag = 0;
    # Intron start and end
    my $intron_beg = 0;
    my $intron_end = 0;
    #my $genes_i = 0;
	# Go through file data
	foreach my $line (@in_data) {
		# Get line data
        my @d = split(/\t/,$line);
        my $chr = $d[0];
        my $typ = $d[2];
        my $beg = $d[3];
        my $end = $d[4];
        # Start of a new gene / transcript
        if ($typ eq 'gene' || $typ eq 'mRNA') {
            $start_flag = 1;
            #$genes_i++ if $typ eq 'gene';
            #if ($genes_i > 1) { exit; }
        }
        # First exon of current gene
        if ($start_flag && $typ eq 'exon') {
            # Start of first intron
            $intron_beg = $end+1;
            # Reset gene start flag
            $start_flag = 0;
        }
        # Following exon of current gene
        elsif ($typ eq 'exon') {
            # End of current intron
            $intron_end = $beg-1;
            #print("$chr\t$intron_beg\t$intron_end\n");
            # Save intron positions
            foreach my $pos ($intron_beg..$intron_end) {
                $intron_positions{$chr}{$pos} = 1;
            }
            # Start of next intron
            $intron_beg = $end+1;
        }
        # Save all exon positions
        if ($typ eq 'exon') {
            foreach my $pos ($beg..$end) {
                $exon_positions{$chr}{$pos} = 1;
            }
        }
    }
    # Remove exonic positions from intron positions
    foreach my $chr (sort keys %intron_positions) {
        foreach my $pos (sort keys %{$intron_positions{$chr}}) {
            if ($exon_positions{$chr}{$pos}) {
                delete($intron_positions{$chr}{$pos});
            }
        }
    }
	# Return data fields
	return \%intron_positions;
}

sub get_repeat_positions {
	# Take repeatmasker file name
	my($repmask_file) = @_;
	# Storage variable
	my %rep_positions = ();
	# Get file data
	my @repmask_data = get_file_data_array($repmask_file);
	# Parse repeatmasker file
	foreach my $line (@repmask_data) {
		# Line starts with sw score
		if ($line =~ /^\s*\d+/) {
			# Remove leading whitespace
			$line =~ s/^\s*//;
			# Get line data
			my @d = split(/\s+/,$line);
			my $chr = $d[4];
            my $beg = $d[5];
            my $end = $d[6];
            my $typ = $d[10];
            # Skip non-interspersed repeats
            next if $typ eq 'Satellite' || $typ eq 'Simple_repeat' || $typ eq 'Low_complexity';
            # Save repeat positions
            foreach my $pos ($beg..$end) {
                $rep_positions{$chr}{$pos} = 1;
            }
		}
	}
	return \%rep_positions;
}


# Save fasta data as hash
# Usage: my $sequences = get_fasta_seqs($infile);
sub get_fasta_seqs {
	# Take fasta file name
	my($file,$short) = @_;
	# Variables
	my $name = '';
	my %sequences = ();
	# Open file
	my $in = open_infile($file);
	# Extract sequence and save with name
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) {
			($name) = ($line =~ />(.*)$/);
			if ($short) { $name =~ s/\s.*// }
		} else {
			$sequences{$name} .= $line;
		}
	}
	return \%sequences;
}

# Open input file
# Usage: my $in = open_infile($infile);
sub open_infile {
	# Take input file name
    my($file) = @_;
    # Open input file
    my $fh;
    if ($file =~ /.gz$/) {
		open($fh, "gunzip -c $file |") or die("Cannot open file '$file': $!\n");
	} else {
    	open($fh, '<', $file) or die("Cannot open file '$file': $!\n");
    }
    # Return filehandle
    return $fh;
}

# Open output file
# Usage: my $out = open_outfile($outfile);
sub open_outfile {
	# Take output file name
    my($file) = @_;
    # Open output file
    open(my $fh, '>', $file) or die("Cannot open file '$file': $!\n");
    # Return filehandle
    return $fh;
}

# Find file name that does not already exist
# Usage: my $filename = find_unique_filename($filename);
sub find_unique_filename {
    # Take file name
    my($filename) = @_;
    # Check if file already exists
    if (-e $filename) {
    	# Filename contains a suffix
    	if ($filename =~ /.+\..+/) {
			my($prefix) = ($filename =~ /(.*)\.[^.]*$/);
			my($suffix) = ($filename =~ /.*(\.[^.]*$)/);
			my $orig_prefix = $prefix;
			# Find a non-existent name
			my $count = 0;
			while (-e $filename) {
				$count++;
				$prefix = $orig_prefix.'_'.$count;
				$filename = $prefix.$suffix;
			}
		}
		# Filename does not contain a suffix
		else {
			my $orig_name = $filename;
			# Find a non-existent name
			my $count = 0;
			while (-e $filename) {
				$count++;
				$filename = $orig_name;
				$filename = $orig_name.'_'.$count;
			}
		}
    }
    # Return unique filename
    return $filename;
}

# Extract file data and save in array
# Usage: my @filedata = get_file_data_array($file);
sub get_file_data_array {
	# Take input file name
    my($file,$ref_opt) = @_;
    my @filedata = ();
    $ref_opt = 0 unless $ref_opt;
	# Open input file
    my $fh = open_infile($file);
	# Extract lines and save in array
    while (my $line = <$fh>) {
    	$line =~ s/\s+$//; #better chomp
    	push(@filedata, $line);
    }
	# Close file
    close($fh) or die("Unable to close: $!\n");
	# Return array containing file data
    if ($ref_opt) {
    	return \@filedata;
    } else {
    	return @filedata;
    }
}

# Save fasta headers as array
# Usage: my $headers = get_fasta_headers($infile);
sub get_fasta_headers {
	# Take fasta file name
	my($file,$short) = @_;
	# Variables
	my $header = '';
	my @headers = ();
	# Open file
	my $in = open_infile($file);
	# Extract header and save in array
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) {
			($header) = ($line =~ />(.*)$/);
			if ($short) { $header =~ s/\s.*// }
			push(@headers,$header);
		}
	}
	return \@headers;
}

# Output reverse complement of input DNA/RNA
sub rev_com {
    # Input sequence
	my($seq) = @_;
    # Reverse sequence
	my $revcom = reverse $seq;
    # Get complement of sequence
	$revcom =~ tr/ACGTUNacgtun/TGCAANtgcaan/;
    # Remove whitespace
	$revcom =~ s/\s//g;
    # Return reverse complement
	return $revcom;
}

# Generate shortened md5 as sequence id
sub generate_short_md5_id {
    # Input: sequence and id length
    my ($dna_sequence, $length) = @_;
    # Create hash value of sequence
    my $hash = md5_hex($dna_sequence);
    # Return shortened hash of given length
    return substr($hash, 0, $length);
}

sub usage {
	# Take optional message
	my($message) = @_;
	# Usage message
	my $USAGE = "perl $0 -g genome.fas -f genes.gtf -r repeats.out [opts]\n";
	# Handle message
	if ($message) { $message .= "\n" unless $message =~ /\n$/ }
	if (!$message) { $message = '' }
	# Print usage to stderr
	return ($message, "\nUsage: $USAGE\n");
}

sub help {
	# Print usage and help message
	die(usage(),
		"Required:\n",
        "  -genome or -g [file]        genome sequence fasta file\n",
        "  -geneset or -f [file]       genome gene set gtf file\n",
		"  -repeats or -r [file]       genome repeatmasker output file\n",
        "\nOptions:\n",
		"  -incl_chrs [chr1,chr2,...]  included chromosomes (hits only allowed on these)\n",
        "                              (chromosome names, delimited by commas)\n",
        "                              [default: all (if nothing specified)]\n",
        "  -seq_len [int]              total sequence length (excluding PAM)\n",
        "                              [default: 20]\n",
        "  -proximal_len [int]         sequence length of proximal part\n",
        "                              (requiring perfect match)\n",
        "                              [default: 12]\n",
        "  -distal_len [int]           sequence length of distal part\n",
        "                              (allowing n distal mismatches; s. next option)\n",
        "                              [default: 8]\n",
        "  -n_distal_mms [int]         max allowed mismatches in distal part\n",
        "                              [default: 2]\n",
        "  -max_srep [int]             max number of simple repeats in sequence stretch\n",
        "                              (e.g. 4 simple repeats: GGGG or ATATATAT, etc.)\n",
        "                              [default: 4]\n",
        "  -n_seqs [int]               number of initially generated sequences\n",
        "                              higher n = longer run + more extensive search\n",
        "                              [default: 100000]\n",
        "  -max_n_hits [int]           max number of hits for final output\n",
        "                              [default: 30]\n",
        "  -n_sets [int]               number of sequence sets in final output\n",
        "                              [default: 1]\n",
        "  -filt_mm_seqs               filter sequences that have hits with mismatches\n",
        "                              [default: no]\n",
        "  -start_g                    sequences have to start with a G at the 5' end\n",
        "                              [default: no]\n",
		"  -help, -h                   show help message\n",
		"\n",
        "  int = integer\n",
        "\n",
	);
}

################################ documentation ################################

=pod

=head1 DESCRIPTION

The purpose of this script is the generation of guide RNAs (gRNAs) for
CRISPR-Cas9 experiments with a range of different numbers of target
sites in a given genome. For gRNA sequences with multiple target sites
the aim is to spread out target sites across chromosomes as much as 
possible, while only allowing target sites in transposable elements or
gene introns. This ensures that no gene exons or regulatory sequences
will be affected by double strand breaks induced by Cas9.
This script requires at least three parameters. The name of a genome
sequence fasta file, the name of a corresponding geneset GTF file, and 
a repeatmasker output file for the used genome.

You can install Perl modules with CPAN. See
http://www.cpan.org/modules/INSTALL.html

=cut
