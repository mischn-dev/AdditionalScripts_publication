 ~/Dokumente/software/samtools-1.8/samtools mpileup -uvf Hordeum_vulgare.IBSC_v2.dna.toplevel.fa -I -q 25 -Q 30 -t AD --positions parents_reference_SNP_positions_in_Global_ref.txt \
../A_Golf/filtered_reads.bam \
../A_ISR_42-8/filtered_reads.bam \
../F3P1-1/filtered_reads_sorted.bam \
../F3P1-2/filtered_reads_sorted.bam  \
../F3P1-1_2/filtered_reads_sorted.bam  \
../F3P1-2_2/filtered_reads_sorted.bam  \
../F12P1K1/filtered_reads_sorted.bam  \
../F12P1k3/filtered_reads_sorted.bam  \
../F12P1Ö1/filtered_reads_sorted.bam  \
../F12P1Ö2/filtered_reads_sorted.bam \
../F16P1K1/filtered_reads_sorted.bam  \
../F16P1K2/filtered_reads_sorted.bam  \
../F16P1Ö1/filtered_reads_sorted.bam  \
../F16P1Ö2/filtered_reads_sorted.bam  \
../F22P1K1/filtered_reads_sorted.bam  \
../F22P1KEP/filtered_reads_sorted.bam  \
../F22P1Ö1/filtered_reads_sorted.bam  \
../F22P1Ö2/filtered_reads_sorted.bam  \
../F23P1K1/filtered_reads_sorted.bam  \
../F23P1K2/filtered_reads_sorted.bam   \
../F23P1Ö1/filtered_reads_sorted.bam  \
../F23P1Ö2/filtered_reads_sorted.bam  \
 | ~/Dokumente/software/bcftools-1.8/bcftools call -vm -Ov --threads 30  > ../snpcalling_sambamba_unique_known_pos_all_dup_rm.vcf

