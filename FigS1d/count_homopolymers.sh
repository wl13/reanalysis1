

perl -ne 'chomp; next unless(/SNV/); my @line = (split /\,/); print "Chr$line[0]\t$line[1]\t$line[3]\n";' \
    raw_variants.csv | sort -k1,1 -k2,2n \
    > raw_variants.SNVs.pos.csv


awk 'BEGIN{OFS="\t";} !/\#/ && /MA_training/ {print $1,$2-10,$2+10,$3;}' \
    raw_variants.SNVs.pos.csv | \
    fasta_process.pl --query - --fasta TAIR10_chr_all.fasta \
    --rows 0 1 2 3 --subset 1 2 --out-format tabular | sed 's/\_/\t/g' \
    > raw_variants.MA_training_SNVs.ex10.csv

perl -ne 'my @line = (split /\s+/); my @hps = ($line[5] =~ m/(A+)|(G+)|(C+)|(T+)/g);
    my $longest_hp = (sort { length($a) <=> length($b) } (@hps))[-1]; my $llen = length($longest_hp);
    my $pos = $line[1]+10; print "$line[0]\t$pos\t$line[5]\t$longest_hp\t$llen\n";' \
    raw_variants.MA_training_SNVs.ex10.csv \
    > raw_variants.MA_training_SNVs.ex10.hps.csv
