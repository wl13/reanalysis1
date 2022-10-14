#!/usr/bin/perl -w
#
#   MegaNGKaKs.pm -- Use to calculate dN dS values
#   Last Modified By: Nowind
#   Last Modified At: 2010.12.20

package MyPerl::MegaNGKaKs;

@ISA    = qw(Exporter);

@EXPORT = qw(Calc_KaKs_ratio Check_divergence get_pn_ps);


my %genetic_code;
my (%codon2syn_site, %codon2non_site);
my %Codon1_Codon2_syn_non_num_info;

sub initialize
{
    my ($pkg) = shift;
    
    get_genetic_code();
    get_codon2syn_site_AND_codon2non_site();
    get_Codon1_Codon2_syn_non_num_info();
    
    return 1;
}

#Calculate Ka Ks Ka/Ks and Pi Value
sub Calc_KaKs_ratio
{
    my ($pN_num, $pN_site, $pS_num, $pS_site) = @_;
    
    my $dN = 0 - 3/4 * log(1 - 4/3 * $pN_num/$pN_site);
    my $dS = 0 - 3/4 * log(1 - 4/3 * $pS_num/$pS_site);
    
    my $ratio = $dS ? $dN/$dS : 99;
       $ratio = 0 if (!$dN && !$dS);
    
    return ($dN, $dS, $ratio);
}

#Check if num/site > 3/4
sub Check_divergence
{
    my ($num, $site) = @_;
    
    if ($num/$site < 0.75) {
        return 1;
    } else {
        return 0;
    }
}

#Calc the pN pS Value
sub get_pn_ps
{
    my ($seq1, $seq2) = @_;
    
    my ($pN_num, $pN_site1, $pN_site2, $pS_num, $pS_site1, $pS_site2) = (0,0,0,0);
    #print ">Seq1\n", $seq1, "\n", ">Seq2\n", $seq2, "\n\n";
    
    $seq1 = uc $seq1;
    $seq2 = uc $seq2;
    
    my $align_len_seq1 = length $seq1;
    #print "\rLength: $align_len_seq1";
    my (@base1, @base2) = ();
    foreach my $index ( 0 .. $align_len_seq1 - 1){
        $base1[$index] = substr($seq1,$index,1);
        $base2[$index] = substr($seq2,$index,1);
    }
    
    my $triple_switch = 1;
    my ($codon1, $codon2) = ($base1[0], $base2[0]);#OMG, So When U Got --ATG... , Then U May
                                                   #Got some surprises(= =), Just Be Careful, It
                                                   #Dose Make Difference!!!
    #print "Codon1: $codon1; Codon2: $codon2\n";
    foreach my $index ( 1 .. $align_len_seq1 - 1){
        my $b1 = $base1[$index];
        my $b2 = $base2[$index];
        if ($b1 ne '-'){
            if ($triple_switch < 3){
                #print $codon1, '.=', $b1, "\n"; sleep 1;
                $codon1 .= $b1;
                $codon2 .= $b2;
                $triple_switch ++;
            }
            if ($triple_switch == 3){
                if ($genetic_code{$codon1} and $genetic_code{$codon2}
                    and $genetic_code{$codon1} ne '_'
                    and $genetic_code{$codon2} ne '_'){
                    my $codon1_non = $codon2non_site{$codon1};
                    my $codon2_non = $codon2non_site{$codon2};
                    my $codon1_syn = $codon2syn_site{$codon1};
                    my $codon2_syn = $codon2syn_site{$codon2};
                    $pN_site1 += $codon1_non;
                    $pN_site2 += $codon2_non;
                    $pS_site1 += $codon1_syn;
                    $pS_site2 += $codon2_syn;

                    my $pS_num_this_codon = 0;
                    my $pN_num_this_codon = 0;
                    if (get_diff_base_number($codon1,$codon2) == 0){ # 0 base diff in codon
                    }
                    elsif (get_diff_base_number($codon1,$codon2) >= 1){
                        $pS_num += $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2}->[0];          
                        $pN_num += $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2}->[1];
                    }   
                }
                $codon1 = '';
                $codon2 = '';
                $triple_switch = 0;
            }
        }
    }
    my $pN_site = ($pN_site1 + $pN_site2) / 2;
    my $pS_site = ($pS_site1 + $pS_site2) / 2;
    
    return ($pN_num, $pN_site, $pS_num, $pS_site);
}


sub get_intermediate_codon {
    my ($codon1, $codon2) = @_;
    my $diff_num = 0;
    my %hash = ();
    for my $index (0..2){
        my $base1 = substr($codon1,$index,1);
        my $base2 = substr($codon2,$index,1);
        $hash{$index} = [$base1, $base2] if $base1 ne $base2;
    }
    
    my @intermediate_codons = ();
    foreach my $ind (keys %hash){
        my $new_codon1 = $codon1;
        substr($new_codon1,$ind,1,$hash{$ind}->[1]);
        if ($new_codon1 ne 'TAA' and
            $new_codon1 ne 'TAG' and
            $new_codon1 ne 'TGA'){
            push @intermediate_codons,$new_codon1;
        }
    }
    return @intermediate_codons;
}

sub get_diff_base_number {
    my ($codon1, $codon2) = @_;
    my $diff_num = 0;
    for (0..2){
        $diff_num++ if substr($codon1,$_,1) ne substr($codon2,$_,1);
    }
    return $diff_num;
}


sub get_genetic_code {
    %genetic_code  = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

}

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my ( %genetic_code ) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if( exists $genetic_code{$codon} ) {
        return $genetic_code{$codon};
    }
    else{
        print STDERR "Bad codon \"$codon\"!!\n";
        exit;
    }
}


sub get_Codon1_Codon2_syn_non_num_info
{     
    foreach my $codon1 (grep !/TAA/, grep !/TAG/, grep !/TGA/, keys %genetic_code)
    {   # Stop Codon is not count, followed the rules of Mega4
        foreach my $codon2 (grep !/TAA/, grep !/TAG/, grep !/TGA/, keys %genetic_code)
        {   #Stop Codon is not count, followed the rules of Mega4
            if ($codon1 ne $codon2) {
                if (get_diff_base_number($codon1, $codon2) == 1) {  # 2 Codon Diff 1
                    if (codon2aa($codon1) eq codon2aa($codon2)) {
                        $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2} = [1,0];
                    }
                    else {
                        $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2} = [0,1];
                    }
                }
                elsif (get_diff_base_number($codon1, $codon2) == 2) {  # 2 Codon Diff 1
                    #print "$codon1 $codon2\n"; 
                    my @intermediate_codons = get_intermediate_codon($codon1,$codon2);
                    my $syn_way_number = 0;
                    my $non_way_number = 0;
                    foreach my $intermediate_codon (@intermediate_codons)
                    {
                        #next if codon2aa($intermediate_codon) eq '_';
                        if (codon2aa($codon1) eq codon2aa($intermediate_codon)){ 
                            $syn_way_number += 1;
                        }
                        elsif (codon2aa($codon1) ne codon2aa($intermediate_codon)){
                            $non_way_number += 1;
                        }
                        if (codon2aa($intermediate_codon) eq codon2aa($codon2)){ 
                            $syn_way_number += 1;
                        }
                        elsif (codon2aa($intermediate_codon) ne codon2aa($codon2)){
                            $non_way_number += 1;
                        }                        
                    }
                    my $syn_way_number_final = 0;
                    my $non_way_number_final = 0;
                    if (($syn_way_number + $non_way_number) != 0) {
                        $syn_way_number_final = 2 * $syn_way_number / ($syn_way_number + $non_way_number);
                        $non_way_number_final = 2 * $non_way_number / ($syn_way_number + $non_way_number);
                    }
                    $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2} = [$syn_way_number_final,$non_way_number_final];
                }
                elsif (get_diff_base_number($codon1, $codon2) == 3) {  # 2 Codon Diff 1
                    my @intermediate_codons = get_intermediate_codon($codon1,$codon2);
                    my $syn_way_number = 0;
                    my $non_way_number = 0;
                    my $empty_syn_and_non_number = 0;
                    foreach my $intermediate_codon (@intermediate_codons)
                    {   # AAA - 'BBB' - CCC - DDD
                        # next if codon2aa($intermediate_codon) eq '_';                   
                        my @intermediate_codons_further = get_intermediate_codon($intermediate_codon,$codon2);
                        foreach my $intermediate_codon_further (@intermediate_codons_further)
                        {   # AAA - BBB - 'CCC' - DDD
                            # next if codon2aa($intermediate_codon_further) eq '_';
                            if (codon2aa($codon1) eq codon2aa($intermediate_codon)) { 
                                $syn_way_number += 1;
                            }
                            elsif (codon2aa($codon1) ne codon2aa($intermediate_codon)) {
                                $non_way_number += 1;
                                #print "+\n";
                            }     
                            
                            #print "$codon1 ", codon2aa($codon1), " $intermediate_codon ", codon2aa($intermediate_codon), " $intermediate_codon_further ", codon2aa($intermediate_codon_further), " $codon2 ", codon2aa($codon2), "\n";
                            if (codon2aa($intermediate_codon) eq codon2aa($intermediate_codon_further)) { 
                                $syn_way_number += 1;
                            }
                            elsif (codon2aa($intermediate_codon) ne codon2aa($intermediate_codon_further)) {
                                $non_way_number += 1;
                            }
                            if (codon2aa($intermediate_codon_further) eq codon2aa($codon2)) { 
                                $syn_way_number += 1;
                            }
                            elsif (codon2aa($intermediate_codon_further) ne codon2aa($codon2)) {
                                $non_way_number += 1;
                            }
                        }
                    }
                    my $syn_way_number_final = 0;
                    my $non_way_number_final = 0;
                    #print "Syn: $syn_way_number\n", "Non: $non_way_number\n\n";
                    if (($syn_way_number + $non_way_number) != 0) {
                        $syn_way_number_final = 3 * $syn_way_number / ($syn_way_number + $non_way_number);
                        $non_way_number_final = 3 * $non_way_number / ($syn_way_number + $non_way_number);
                    }
                    $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2} = [$syn_way_number_final,$non_way_number_final];
                }
            }
        }
    }
}

sub get_codon2syn_site_AND_codon2non_site {

    %codon2syn_site = (
        'AAA' => 1/3,
        'AAC' => 1/3,
        'AAG' => 1/3,
        'AAT' => 1/3,
        'ACA' => 1,
        'ACC' => 1,
        'ACG' => 1,
        'ACT' => 1,
        'AGA' => 5/6,
        'AGC' => 1/3,
        'AGG' => 2/3,
        'AGT' => 1/3,
        'ATA' => 2/3,
        'ATC' => 2/3,
        'ATG' => 0,
        'ATT' => 2/3,
        'CAA' => 1/3,
        'CAC' => 1/3,
        'CAG' => 1/3,
        'CAT' => 1/3,
        'CCA' => 1,
        'CCC' => 1,
        'CCG' => 1,
        'CCT' => 1,
        'CGA' => 3/2,
        'CGC' => 1,
        'CGG' => 4/3,
        'CGT' => 1,
        'CTA' => 4/3,
        'CTC' => 1,
        'CTG' => 4/3,
        'CTT' => 1,
        'GAA' => 1/3,
        'GAC' => 1/3,
        'GAG' => 1/3,
        'GAT' => 1/3,
        'GCA' => 1,
        'GCC' => 1,
        'GCG' => 1,
        'GCT' => 1,
        'GGA' => 1,
        'GGC' => 1,
        'GGG' => 1,
        'GGT' => 1,
        'GTA' => 1,
        'GTC' => 1,
        'GTG' => 1,
        'GTT' => 1,
        'TAC' => 1,
        'TAT' => 1,
        'TCA' => 1,
        'TCC' => 1,
        'TCG' => 1,
        'TCT' => 1,
        'TGC' => 1/2,
        'TGG' => 0,
        'TGT' => 1/2,
        'TTA' => 2/3,
        'TTC' => 1/3,
        'TTG' => 2/3,
        'TTT' => 1/3
    );
    foreach my $codon (keys %codon2syn_site){
        $codon2non_site{$codon} = 3 - $codon2syn_site{$codon};
    }
}

1;

__END__

#sub get_codon2syn_site_AND_codon2non_site {     #use decimal instead of fraction
#    my %codon2syn_site = ();
#    my %codon2non_site = ();
#
#    %codon2syn_site = (
#        'AAA' => 0.333,
#        'AAC' => 0.333,
#        'AAG' => 0.333,
#        'AAT' => 0.333,
#        'ACA' => 1,
#        'ACC' => 1,
#        'ACG' => 1,
#        'ACT' => 1,
#        'AGA' => 0.833,
#        'AGC' => 0.333,
#        'AGG' => 0.667,
#        'AGT' => 0.333,
#        'ATA' => 0.667,
#        'ATC' => 0.667,
#        'ATG' => 0,
#        'ATT' => 0.667,
#        'CAA' => 0.333,
#        'CAC' => 0.333,
#        'CAG' => 0.333,
#        'CAT' => 0.333,
#        'CCA' => 1,
#        'CCC' => 1,
#        'CCG' => 1,
#        'CCT' => 1,
#        'CGA' => 1.5,
#        'CGC' => 1,
#        'CGG' => 1.333,
#        'CGT' => 1,
#        'CTA' => 1.333,
#        'CTC' => 1,
#        'CTG' => 1.333,
#        'CTT' => 1,
#        'GAA' => 0.333,
#        'GAC' => 0.333,
#        'GAG' => 0.333,
#        'GAT' => 0.333,
#        'GCA' => 1,
#        'GCC' => 1,
#        'GCG' => 1,
#        'GCT' => 1,
#        'GGA' => 1,
#        'GGC' => 1,
#        'GGG' => 1,
#        'GGT' => 1,
#        'GTA' => 1,
#        'GTC' => 1,
#        'GTG' => 1,
#        'GTT' => 1,
#        'TAC' => 1,
#        'TAT' => 1,
#        'TCA' => 1,
#        'TCC' => 1,
#        'TCG' => 1,
#        'TCT' => 1,
#        'TGC' => 0.5,
#        'TGG' => 0,
#        'TGT' => 0.5,
#        'TTA' => 0.667,
#        'TTC' => 0.333,
#        'TTG' => 0.667,
#        'TTT' => 0.333
#    );
#    foreach my $codon (keys %codon2syn_site){
#        $codon2non_site{$codon} = 3 - $codon2syn_site{$codon};
#    }
#    return (\%codon2syn_site, \%codon2non_site);
#}
#
