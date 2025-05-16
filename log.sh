
#awk -F "\t" '{print $1"\t"$4"\t"$3}' index_primer.txt|sed 's/ Mogolia_CAUNMG_seq04//g' |grep -v '^order' > index.tsv
#head -96 index.tsv|cut -f 1,2 |sed 's/^0//g' > cell.index.tsv
#cut -f 3 index.tsv|uniq |awk '{print "P"NR"\t"$1}' >plate.index.tsv
python3 /Users/yangchentao/Desktop/Work/Cyclone/bin/CycFqFilter.py -q 7 -l 700 -L 770 -g 0.2 -G 0.6 -o test.clean test.fastq.gz
seqkit fq2fa -w 0 -o test.clean.fa test.clean.fq.gz
python pcr_demultiplex.py -p primer.txt --plate-index plate.index.tsv --well-index cell.index.tsv -f test.clean.fa  -o output
