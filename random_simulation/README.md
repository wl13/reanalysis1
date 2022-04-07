## Code and data used for random simulation 


The PERL scripts require the "MyPerl" module to be included in the PERL5LIB path, or be placed at the same folder with the scripts.

<br />

Run perl script without any parameters will return the usage. 

<br />


### Additional tools used:
(1) Linux awk;

(2) samtools v1.12 (http://www.htslib.org/);

(3) bedtools v2.30.0 (https://github.com/arq5x/bedtools2)

(4) Other required PERL modules: "Getopt::Long", "File::Find::Rule", "List::Util::WeightedChoice", "Data::Random"

<br />


### Detailed procedures

<br />

#### Step1: Run extract_AT_homomeric_regions.pl to get the positions of A/T homomeric regions, i.e., regions with >= 3 As or >= Ts:

		extract_AT_homomeric_regions.pl --fasta TAIR10_chr_all.fasta \
			> AT_homopolymers.csv
    
<br />
   
#### Step2: The above step gives 1-based start position, change it to 0-based bed format

		awk 'BEGIN{OFS="\t";} !/\#/ {print $1,$2-1,$3,$4"x"$5;}' \
			AT_homopolymers.csv > AT_homopolymers.bed

<br />

#### Step3: Obtain non-AT-homomeric regions through subtracting 

		samtools faidx TAIR10_chr_all.fasta
		
		awk 'BEGIN{OFS="\t";} !/\#/ {print $1,0,$2;}' \
			TAIR10_chr_all.fasta.fai > TAIR10_chr_all.bed

		bedtools subtract -nonamecheck -a TAIR10_chr_all.bed \
    		-b AT_homopolymers.bed > non_AT_homopolymers.bed
	
<br />

#### Step4: Mask the non-AT-homomeric regions in the genome fasta file with Ns, this produce a new genome fasta with only AT homomeric sequences

		bedtools maskfasta -fi TAIR10_chr_all.fasta \
		    -fo genomic_AT_homopolymers.fasta \
		    -bed non_AT_homopolymers.bed


<br />

#### Step5: Generate random positions across the new AT-homomeric only genome using "simulate_mutation_pos.pl", this script only simulate mutations in A/T/G/C sites so Ns are ignored. A pipe was used to further reformat the output for plotting purpose. 

				
		simulate_mutation_pos.pl --fasta genomic_AT_homopolymers.fasta \
		    --random-times 1000 --random-size 2000 --exclude mitochondria chloroplast | \
		    perl -ne 'next if (/\#\#/); if (/\#/){print "CHROM,POS,TYPE,src\n"; next;} 
				my @line = (split /\s+/); $line[1] =~ s/Chr//; 
				print "$line[1],$line[2],SNV,set$line[0]\n";' \
		    > 2000000_random_positions_in_AT_homomeric_regions.csv

<br />


