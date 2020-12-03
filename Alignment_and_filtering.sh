for i in $(find  WGS_Population_1/F* -maxdepth 0 -type d)
do
cd $i 
echo $i
echo "align, sort and filter - not removal!!!"

# two data sets have to be matched together
echo "merge the datasets togeter"
cat  *_1.fq.gz > R1.fq.gz
cat  *_2.fq.gz > R2.fq.gz

# align the reads and filter directly, save in bam format to save disc space
echo "align and sort the files"
bwa mem -t 30 /media/michael/Daten/Reference_genomes/Barley_30052019/bwa_index/balrey R1.fq.gz R2.fq.gz | samtools sort -t 30 -o output_sorted.bam -

#echo "filter the alignments"
sambamba_v0.6.6 view -t 30 -h -f bam -p -F "mapping_quality >= 30 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" -o filtered_reads.bam output_sorted.bam

# remove the duplicates and sort
echo "remove the duplicates"
sambamba_v0.6.6 markdup -r -t 30 -p filtered_reads.bam sambamba_filtered_markdup.bam

sambamba_v0.6.6 sort -t 30 -p -o filtered_reads_sorted.bam sambamba_filtered_markdup.bam

done

