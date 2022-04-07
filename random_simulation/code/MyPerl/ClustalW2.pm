#!/usr/bin/perl -w
#
#   ClustalW2.pm -- Run ClustalW2 to draw phylogenetic tree
#   Author: Nowind
#   Created: 2010-11-14
#   Last Modified: 2010-12-01

use strict;

package MyPerl::ClustalW2;

use Bio::SeqIO;

our (@ISA, @EXPORT);

@ISA    = qw(Exporter);

@EXPORT = qw(Draw_NJ_Tree);

######################### Main #########################

sub Draw_NJ_Tree
{
    my ($input) = @_;
    
    #
    # Translate nucleotide sequences to protein sequences
    #
    
    print "Translating to protein sequences...";
    
    $input =~ /^(.*)\./;
    my $protein = $1 . '.pro.fasta';
    
    my %nt_seq;
    
    my $seqio = Bio::SeqIO->new( -format => 'fasta', -file => $input );
    my $seqot = Bio::SeqIO->new( -format => 'fasta', -file => "> $protein" );
    
    while ( my $seq = $seqio->next_seq() )
    {
        $nt_seq{$seq->display_id} = \($seq->seq);
        
        my $trans = $seq->translate();
        $seqot->write_seq($trans);
    }
    
    print "done!\n";
    
    
    #
    # Align with the protein sequences
    #
    
    my $aln_pro = $protein . '.aln';
    Align_AA($protein, $aln_pro);
    
    #
    # ReTranslate aligned protein sequences
    # to aligned nucleotide sequences
    #
    
    print "Retranslate to nucleotide sequences...";
    
    my $aa_io = Bio::SeqIO->new( -format => 'fasta', -file => $aln_pro );
    my $nt_ot = $aln_pro . '.nt';
    
    open (W, "> $nt_ot") or die $!;
    while ( my $aln_aa = $aa_io->next_seq() )
    {
        my $raln_nt_seq = ReTrans( $nt_seq{$aln_aa->display_id}, \($aln_aa->seq) );
        
        print W ">" . $aln_aa->display_id . "\n";
        print W $$raln_nt_seq . "\n";
    }
    close W;
    
    print "done!\n";
    
    #
    # Generate phylogenetic trees use
    # aligned nucleotide sequences
    #
    
    Draw_NT($nt_ot);
}


sub Align_AA
{    
    my ($infile, $outfile) = @_;

    print "Start alignment...\n";
    
    my $aln = <<CMD;
    clustalw2 -INFILE=$infile
              -ALIGN
              -TYPE=PROTEIN
              -OUTFILE=$outfile
              -OUTPUT=FASTA
              -QUIET
CMD
    
    system $aln;
    
    print "done!\n";
}

sub Draw_NT
{
    my ($infile) = @_;
    
    print "Start drawing trees...\n";
    
    my $draw = <<CMD;
    clustalw2 -INFILE=$infile
              -OUTPUTTREE=phylip
              -BOOTSTRAP=1000
              -SEED=111
              -BOOTLABELS=branch
              -CLUSTERING=NJ
              -QUIET
CMD

    system $draw;
    
    print "done!\n";
}

sub ReTrans
{
    my ($rnt, $raln_aa) = @_;
    
    my @aa = split //, $$raln_aa;
    
    my $aln_nt = '';
    my $pos    = 0;
    
    for (my $i=0; $i<= $#aa; $i++)
    {
        if ( $aa[$i] eq '-' ) {
            $aln_nt .= '-' x 3;
        }
        else {
            my $codon = substr($$rnt, $pos, 3);
            $aln_nt .= $codon;
            
            $pos += 3;
        }
    }
    
    return \$aln_nt;
}

sub usage
{
    my $msg = <<EOF;
    
    Run_ClustalW.pl -- Run clustalw2 to align sequences by codons, and
                       generate phylogenetic trees use NJ method.
                       
    Usage:
    perl Run_ClustalW.pl -in <fasta file>
    
EOF

    print $msg;
    exit(1);
}

1;