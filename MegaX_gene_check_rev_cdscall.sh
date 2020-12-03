# input $1 = position - $2 = gene name - $3 = strand specifity (1 => reverse , all other is forward [0])

# get the reference sequence for the postion of interest - if "-" use the reverse complement strand
samtools faidx /media/michael/Daten/Reference_genomes/Barley_30052019/Hordeum_vulgare.IBSC_v2.dna.toplevel.fa $1 > ref_$2_gene.fa

# create a one liner file 
../FASTX_Toolkit/fasta_formatter -i ref_$2_gene.fa -o ref_$2_gene_2.fa -w 0 

# get the reverse strand
if [ $3 -eq 1 ]
then
../FASTX_Toolkit/fastx_reverse_complement -i ref_$2_gene_2.fa -o ref_$2_gene_2.fa
fi  

##################### golf ################################

# get the reads aligned at the position of interest
# EXAMPLE ~/Dokumente/software/samtools-1.8/samtools view output_sorted.bam chr5H:632099213-632101201 -o phosphate_gene.bam
samtools view ../A_Golf/output_sorted.bam $1 -o golf.$2.bam

# call the bases for each base locus along the gene
bcftools mpileup -Ou -f Hordeum_vulgare.IBSC_v2.dna.toplevel.fa golf.$2.bam | bcftools call -c -Oz -o golf.$2.vcf.gz

# index
bcftools index golf.$2.vcf.gz

# get the sequence 
cat ref_$2_gene.fa | bcftools consensus golf.$2.vcf.gz > golf.consensus$2.fa

# create a one liner
../FASTX_Toolkit/fasta_formatter -i golf.consensus$2.fa -o golf.consensus$2_2.fa -w 0 

# get the reverse complement strand if asked for
if [ $3 -eq 1 ]
then
../FASTX_Toolkit/fastx_reverse_complement -i golf.consensus$2_2.fa -o golf.consensus$2_2.fa
fi

######################## isr 42-8 ############################

# get the reads aligned at the position of interest
samtools view ../A_ISR_42-8/output_sorted.bam $1 -o isr.$2.bam

# call the bases for each base locus along the gene
bcftools mpileup -Ou -f Hordeum_vulgare.IBSC_v2.dna.toplevel.fa isr.$2.bam | bcftools call -c -Oz -o isr.$2.vcf.gz

# index
bcftools index isr.$2.vcf.gz

# get the sequence 
cat ref_$2_gene.fa | bcftools consensus isr.$2.vcf.gz > isr.consensus$2.fa

# create a one liner
../FASTX_Toolkit/fasta_formatter -i isr.consensus$2.fa -o isr.consensus$2_2.fa -w 0 

# get the reverse complement strand if asked for
if [ $3 -eq 1 ]
then
../FASTX_Toolkit/fastx_reverse_complement -i isr.consensus$2_2.fa -o isr.consensus$2_2.fa
fi


#######################################################################
rm ref_$2_gene.fa ref_$2_gene_2.fa golf.$2.bam isr.$2.bam golf.$2.vcf.gz isr.$2.vcf.gz golf.$2.vcf.gz.csi isr.$2.vcf.gz.csi isr.consensus$2.fa  golf.consensus$2.fa
