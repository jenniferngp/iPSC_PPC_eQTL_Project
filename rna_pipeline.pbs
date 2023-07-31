#$ -pe smp 16
#$ -V
#$ -cwd
#$ -e logs
#$ -o logs

hostname >& 2
export PATH=/software/biobambam2-2.0.95/bin:/software/STAR-2.7.3a/bin/Linux_x86_64:/software/samtools-1.9:/software/sambamba_v0.6.7:/software/rsem-1.2.20:$PATH

if [ $# -lt 2 ]; then
        echo "$0 OUT_DIR FASTQ1 [FASTQ2]" >& 2;
        exit 1;
fi
out_dir=$1; shift
in_files=( "$@" )
reference_star=/reference/private/STAR/hg19_gencode34
reference_rsem=/reference/private/RSEM/hg19_gencode34/hg19_gencode34

echo $in_files >& 2
echo $dir >& 2

# 1. Copy reference STAR
date >& 2
reference=$reference_star

if [ ! -d $out_dir ]; then mkdir $out_dir; fi

# 2. Run STAR Alignment
cmd="STAR --runThreadN 16 --genomeDir ${reference} --genomeLoad NoSharedMemory \
--readFilesCommand zcat --readFilesIn ${in_files[@]} --outSAMattributes All \
--outSAMunmapped Within \
--outSAMattrRGline ID:1 PL:ILLUMINA PU:CARDIPS LB:${name} SM:${name} \
--outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --alignIntronMin 20 \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMtype BAM Unsorted \
--outFileNamePrefix ${out_dir}/ --quantMode TranscriptomeSAM"
echo $cmd >& 2; $cmd

# 3. Sort bam 
date >& 2
#cmd="sambamba sort -m 32GB -t 16 --tmpdir ${dir} ${dir}/Aligned.out.bam"
cmd="samtools sort -m 2G -n -o ${out_dir}/Aligned.out.namesorted.bam -@ 16 ${out_dir}/Aligned.out.bam"
echo $cmd >& 2; $cmd

# 4. Fill in mate coordinates
date >& 2
cmd="samtools fixmate -@ 16 -m ${out_dir}/Aligned.out.namesorted.bam ${out_dir}/Aligned.out.namesorted.fixmate.bam"
echo $cmd >& 2; $cmd

# 5. Sort bam
date >& 2
cmd="samtools sort -m 2G -o ${out_dir}/Aligned.out.sorted.bam -@ 16 ${out_dir}/Aligned.out.namesorted.fixmate.bam"
echo $cmd >& 2; $cmd

# 6. Index bam
date >& 2
#cmd="sambamba index -t 16 ${dir}/Aligned.out.sorted.bam"
cmd="samtools index -@ 16 ${out_dir}/Aligned.out.sorted.bam"
echo $cmd >& 2; $cmd

# 7. Mark duplicates
date >& 2
#cmd="bammarkduplicates markthreads=16 tmpfile=${dir} \
#I=${dir}/Aligned.out.sorted.bam O=${dir}/Aligned.out.sorted.mdup.bam \
#M=${dir}/Aligned.out.sorted.mdup.txt"
cmd="samtools markdup -@ 16 -s -T ${out_dir}/temp \
${out_dir}/Aligned.out.sorted.bam ${out_dir}/Aligned.out.sorted.mdup.bam"
echo $cmd >& 2; $cmd

# 8. Index bam
date >& 2
cmd="samtools index -@ 16 ${out_dir}/Aligned.out.sorted.mdup.bam"
#echo $cmd >& 2; $cmd

# 9. Run flagstat
date >& 2
cmd="samtools flagstat -@ 16 ${out_dir}/Aligned.out.sorted.mdup.bam > ${out_dir}/Aligned.out.sorted.mdup.bam.flagstat"
echo $cmd >& 2; eval $cmd

# 10. Run idxstats
date >& 2
cmd="samtools idxstats -@ 16 ${out_dir}/Aligned.out.sorted.mdup.bam > ${out_dir}/Aligned.out.sorted.mdup.bam.idxstats"
echo $cmd >& 2; eval $cmd

# 10. Run samtools stats
date >& 2
cmd="samtools stats -@ 16 ${out_dir}/Aligned.out.sorted.mdup.bam > ${out_dir}/Aligned.out.sorted.mdup.bam.stats"
echo $cmd >& 2; eval $cmd

# 11. Reformat header
date >& 2
cmd="samtools view -H ${out_dir}/Aligned.out.sorted.mdup.bam | sed 's,^@RG.*,@RG\tID:None\tSM:None\tLB:None\tPL:Illumina,g' | samtools reheader - ${out_dir}/Aligned.out.sorted.mdup.bam > ${out_dir}/Aligned.out.sorted.mdup.reheader.bam"
echo $cmd >& 2; eval $cmd

# 12. Copy Gencode Reflat
date >& 2
refFlat=/reference/private/Gencode.v34lift37/refFlat.txt.gz

# 13. Run Picard's CollectRnaSeqMetrics
cmd="java -jar /software/picard-2.20.1/picard.jar CollectRnaSeqMetrics \
I=${out_dir}/Aligned.out.sorted.mdup.reheader.bam \
O=${out_dir}/picard.RNA_Metrics \
VALIDATION_STRINGENCY=SILENT \
REF_FLAT=$refFlat \
STRAND=NONE"
echo $cmd >& 2; $cmd

# 14. Copy reference for RSEM
reference=$reference_rsem

# 15. Run RSEM
date >& 2
cmd="rsem-calculate-expression --bam --num-threads 16 --no-bam-output --seed 3272015 \
--estimate-rspd --forward-prob 0"
if [ ${#in_files[@]} -gt 1 ]; then
        cmd=$cmd" --paired-end"
fi
cmd=$cmd" ${out_dir}/Aligned.toTranscriptome.out.bam ${reference} ${out_dir}/rsem"
echo $cmd >& 2; $cmd

date >& 2
