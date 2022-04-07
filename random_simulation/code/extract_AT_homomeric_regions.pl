#!/usr/bin/perl -w
#
#   extract_AT_homomeric_regions.pl
#
#
#   Author: wanglong@nju.edu.cn
#   Created: 2022-04-04
#   Updated: 2022-04-04
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 22/04/04: The initial version.


=head1 NAME

extract_AT_homomeric_regions.pl


=head1 SYNOPSIS

  extract_AT_homomeric_regions.pl --help/?

=head1 DESCRIPTION

Extract positions of A/T homomeric runs, e.g., AAA.. or TTT... (size >= 3bp).

=cut


use strict;

use Getopt::Long;
use File::Find::Rule;

use MyPerl::FileIO qw(:all);

################### Main #################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my $max_threads = 1;
my ($fasta_file, $out_locations, @out_triplets, $output);
GetOptions(
            "fasta=s"                => \$fasta_file,
            "O|output=s"             => \$output,
            
            "positions"              => \$out_locations,
            "triplets=s{,}"          => \@out_triplets,
           );

unless( $fasta_file ) {
    print <<EOF;

$0  -- Extract positions of A/T homomeric runs, e.g., AAA.. or TTT... (size >= 3bp).

Version: $VERSION

Usage:   perl $0 [--fasta FILE | STDIN] [Options]

Options:

    -f, --fasta     <filename>
        sequence file in fasta format, required
    -O, --output    <filename>
        output file, default to STDOUT
    
EOF

    exit(0);
}

$|++;

print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";


if ($output) {
    open (STDOUT, "> $output") || die $!;
}


search_homopolymers($fasta_file);



print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################


=head2 search_homopolymers

    Usage   : search_homopolymers($fasta_file);
    Args    : Fasta file.
    Returns : Null

=cut
sub search_homopolymers
{
    my ($in) = @_;
    
    print STDERR ">> Start reading $fasta_file ... ";
    my @SEQs = ();
    my @IDs  = parse_fasta_SEQs(\@SEQs, $fasta_file);
    print STDERR "done!\n";
    
    for (my $i=0; $i<@IDs; $i++)
    {
        next if ($IDs[$i] =~ /(mitochondria|chloroplast|scaffold)/);  ## skip mitochondria, chloroplast and etc.
        
        print STDERR "\r>> Start counting ... $IDs[$i]";
        
        next if (length($SEQs[$i]) < 3);
        
        my $seq = uc ($SEQs[$i]);
        
        while($seq =~ /(A{3,}|T{3,})/g)
        {
            my $hp_length = length($1);
            my $hp_end    = pos($seq);
            my $hp_start  = pos($seq) - $hp_length + 1;
            my $hp_base   = substr($1, 0, 1);
            
            print "$IDs[$i]\t$hp_start\t$hp_end\t$hp_base\t$hp_length\n";
        }
    }
    print STDERR "\tdone!\n";
}

