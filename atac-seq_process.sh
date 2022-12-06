1. nohup fastq-dump --gzip --split-3 -A ID -O . &
2. nohup trimmomatic PE -phred33 ID_1.fastq.gz ID_2.fastq.gz -baseout ID.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads 50 &
3. nohup fastqc /ID _1.fastq.gz /ID_2.fastq.gz -t 5 ./ &
4. nohup bowtie2 -p 70 -x /index/hg19 -1 /ID_1P.fq.gz -2 ID_2P.fq.gz | samtools sort -O bam -o ID.bam &
5. gatk MarkDuplicates -I ID.bam -O ID.rmdup.bam --REMOVE _DUPLICATES true -M test.log
6. nohup samtools view -h -f 2 -q 30 ID.rmdup.bam | grep -v chrM | samtools sort -O bam -@ 5-o - > ID.last.bam &
7. bedtools bamtobed -i ID.last.bam  > ID.bed
8. macs2 callpeak -t ID.bed -g hs --nomodel --shift -100 --extsize 200 -n  ID --outdir ./peaks/
9. zcat  gencode.v40.annotation.gtf.gz | grep   protein_coding |perl -alne '{next unless $F[2] eq "gene" ;/gene_name \"(.*?)\";/; print "$F[0]\t$F[3]\t$F[4]\t$1" }' >protein_coding.hg38.position
10. perl -p -i -e 's/ /\t/g' ID_summits.bed
11. bedtools intersect -a ID_summits.bed  -b ./protein_coding.hg38.position  -wa -wb >gene.tsv
~
 
~

