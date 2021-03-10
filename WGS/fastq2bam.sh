#!/bin bash

sample_id=""
threads=""
mem_gb=""
ref_hg=""
outdir=""
R1=""
R2=""

function print_usage() {
/bin/cat << EOF

Usage :
    ${PROGRAM_NAME} [-p prefix] [-t threads] [-m mem_gb] [-r reference] [-o outdir] [-1 Read1] [-2 Read2]
Option :
    -p, prefix    : sample_id for creating output directory and prefix of file name
    -t, threads
    -m, mem_gb
    -r, reference
    -o, outdir    : output abpath
    -1, read1     : fastq R1 abpath
    -2, read2     : fastq R2 abpath

EOF
}

while getopts "p:t:m:r:o:1:2:h" option
do
  case $option in
    p)
      sample_id="$OPTARG"
      ;;
    t)
      threads="$OPTARG"
      ;;
    m)
      mem_gb="$OPTARG"
      ;;
    r)
      ref_hg="$OPTARG"
      ;;
    o)
      outdir="$OPTARG"
      ;;
    1)
      R1="$OPTARG"
      ;;
    2)
      R2="$OPTARG"
      ;;
    h)
      print_usage
      exit 0
      ;;
  esac
done

#echo "smaple_id threads mem_gb ref_hg outdir R1 R2 $smaple_id $threads $mem_gb $ref_hg $outdir $R1 $R2"

mkdir -p ${outdir}/${sample_id}/FASTQC
mkdir -p ${outdir}/${sample_id}/BAM



echo "[INFO] FASTQC"
echo "[INFO] command : python3 /opt/geninus-packages/fastqc.py -r1 ${R1} -r2 ${R2} -t ${threads}  -o ${outdir}/${sample_id}/FASTQC"
python3 /opt/geninus-packages/fastqc.py -r1 ${R1} -r2 ${R2} -t ${threads} -o ${outdir}/${sample_id}/FASTQC



echo "[INFO] trim, cut-adaptor"
echo "[INFO] command : /opt/bbmap/bbduk.sh in=${R1} in2=${R2} out=${outdir}/${sample_id}/FASTQC/${sample_id}.R1.bbduk.r.Q10.fastq.gz out2=${outdir}/${sample_id}/FASTQC/${sample_id}.R2.bbduk.r.Q10.fastq.gz qtrim=r trimq=10"
/opt/bbmap/bbduk.sh in=${R1} in2=${R2} out=${outdir}/${sample_id}/FASTQC/${sample_id}.R1.bbduk.r.Q10.fastq.gz out2=${outdir}/${sample_id}/FASTQC/${sample_id}.R2.bbduk.r.Q10.fastq.gz qtrim=r trimq=10
R1=${outdir}/${sample_id}/FASTQC/${sample_id}.R1.bbduk.r.Q10.fastq.gz
R2=${outdir}/${sample_id}/FASTQC/${sample_id}.R2.bbduk.r.Q10.fastq.gz



echo "[INFO] create bam"
echo "[INFO] command : bwa mem -R \"@RG\\tID:${sampleID}\\tSM:${sampleID}\\tPL:illumina\" -M -t ${threads} ${ref_hg} {R1} {R2} | samtools view -@ ${threads} -Suh -F 12 -f 2 -q 20 - -o ${outdir}/${sample_id}/BAM/${sample_id}.bam"
bwa mem -R "@RG\\tID:${sampleID}\\tSM:${sampleID}\\tPL:illumina" -M -t ${threads} ${ref_hg} {R1} {R2} | samtools view -@ ${threads} -Suh -F 12 -f 2 -q 20 - -o ${outdir}/${sample_id}/BAM/${sample_id}.bam



echo "[INFO] samtools sort bam"
echo "[INFO] command : samtools sort -@ ${threads} ${outdir}/${sample_id}/BAM/${sample_id}.bam -o ${outdir}/${sample_id}/BAM/${smaple_id}.sorted.bam"
samtools sort -@ ${threads} ${outdir}/${sample_id}/BAM/${sample_id}.bam -o ${outdir}/${sample_id}/BAM/${smaple_id}.sorted.bam



echo "[INFO] deduplication"
echo '''[INFO] command :
java8 -Xmx${mem_gb}g \
      -jar /opt/data/picard-2.20.1/picard.jar \
      MarkDuplicates \
      REMOVE_SEQUENCING_DUPLICATES=true \
      INPUT=${outdir}/${sample_id}/BAM/${smaple_id}.sorted.bam \
      OUTPUT=${threads} ${outdir}/${sample_id}/BAM/${sample_id}.dedup.sorted.bam \
      METRICS_FILE=${threads} ${outdir}/${sample_id}/BAM/${sample_id}.dedup.sorted.bam.metrics \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT
'''
java8 -Xmx${mem_gb}g \
      -jar /opt/data/picard-2.20.1/picard.jar \
      MarkDuplicates \
      REMOVE_SEQUENCING_DUPLICATES=true \
      INPUT=${outdir}/${sample_id}/BAM/${smaple_id}.sorted.bam \
      OUTPUT=${threads} ${outdir}/${sample_id}/BAM/${sample_id}.dedup.sorted.bam \
      METRICS_FILE=${threads} ${outdir}/${sample_id}/BAM/${sample_id}.dedup.sorted.bam.metrics \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT
