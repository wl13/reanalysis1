#!/usr/bin/perl -w
#
#   simulate_mutation_pos.pl -- Randomly pick up certain positions across query chromosomes.
#
#   Author: Nowind
#   Created: 2013-09-26
#   Updated: 2022-03-07
#   Version: 2.2.1


use strict;

use Getopt::Long;
use File::Find::Rule;
use Data::Random qw(:all);
use List::Util::WeightedChoice qw( choose_weighted );

use MyPerl::FileIO qw(:all);

################### Main #################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.2.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my %options = ();
   $options{rand_times}   = 1;
   $options{rand_size}  = 100;
GetOptions(
            "fasta=s"             => \$options{fasta_file},
            "output=s"            => \$options{output},
            "T|random-times=i"    => \$options{rand_times},
            "S|random-size=i"     => \$options{rand_size},
            "exclude=s{,}"        => \@{$options{exclude_ids}},
           );

unless( $options{fasta_file} ) {
    print <<EOF;

$0  -- Randomly pick up certain positions across query chromosomes

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -f, --fasta  <filename>
        file contain sequences in fasta format, required
    -o, --output <filename>
        output filename, default to STDOUT
        
    -T, --random-times
        random times [default: 1]
    -S, --random-size
        numbers of mutation loci [default: 100]

    -e, --exclude <strings>
        exclude unwanted chromosomes or scaffolds while simulating, all
        chromosomes with ids match strings specified here would be ignored 
EOF

    exit(0);
}

$|++;

print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($options{output}) {
    open (STDOUT, ">", "$options{output}") or die $!;
}



##
## generate random mutation loci in the genome
##
print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
print STDOUT "#RandomTime\tChrom\tPos\tRef\tMut\n";
my $rh_mut_loci = gen_mut_loci($options{fasta_file}, $options{rand_times}, $options{rand_size});

print STDERR "done!\n";



print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################

=head2 gen_mut_loci

    About   : Generate random positions.
    Usage   : gen_mut_loci($file, $rand_size);
    Args    : Reference sequences in fasta format;
              Number of positions needed to pick out.
    Returns : Randomly generated mutations loci.

=cut
sub gen_mut_loci
{
    my ($in, $rand_times, $rand_size) = @_;
    
    ##
    ## read chromosome sequences
    ##
    print STDERR ">> Start reading reference sequences ... ";
    my %chrom_seqs = ();
    parse_fasta_SEQs(\%chrom_seqs, $in);
    print STDERR "done!\n";
    
    my %chrom_lengths = ();
    my $genome_size   = 0;
    for my $chrom (sort keys %chrom_seqs)
    {
        if (@{$options{exclude_ids}} > 0) {
            my $exclude_str = join '|', @{$options{exclude_ids}};
            next if ($chrom =~ /($exclude_str)/);
        }
        
        $chrom_lengths{$chrom} = ($chrom_seqs{$chrom} =~ tr/ATGCatgc/ATGCatgc/); ## use informative length
        
        $genome_size += $chrom_lengths{$chrom};
    }
    

    my @chrom_ids = sort keys %chrom_lengths;
    my $chrom_num = scalar @chrom_ids;
    
    ##
    ## weight each chromosome by its length
    ##
    my @chrom_weights = ();
    for my $chrom (@chrom_ids)
    {
        my $weight = $chrom_lengths{$chrom} / $genome_size;
        push @chrom_weights, $weight;
    }
    
    
    ##
    ## randomly generate mutation loci
    ##
    my %rand_mutations = ();
    my $rand_index = 0;
    while (++$rand_index <= $rand_times)
    {
        print STDERR "\r>> Start generating random mutations .. $rand_index \/" . $rand_times;
        for (my $i=0; $i<$rand_size; $i++)
        {
            my $chrom = choose_weighted([@chrom_ids], [@chrom_weights]);
            my $pos   = int(rand($chrom_lengths{$chrom}+1));
            
            my $ref_base  = uc(substr($chrom_seqs{$chrom}, $pos-1, 1));
            
            if ($ref_base !~ /A|T|G|C/) { ## skip non-ATGC characters
                $i--;
                next;
            }
            
            my @alt_bases = grep { $_ ne $ref_base } qw(A T G C);
            
            my $mut_base  = $alt_bases[int(rand(3))];
            
            print STDOUT "$rand_index\t$chrom\t$pos\t$ref_base\t$mut_base\n";
        }
    }

    print STDERR "\tdone!\n";
}

