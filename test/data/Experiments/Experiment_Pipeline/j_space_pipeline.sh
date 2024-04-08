echo Script name: $0
echo $# arguments
if [ "$#" -ne 3 ];
then
    echo "j_space_vcf_pipeline: illegal number of parameters"
    echo "SYNOPSIS"
    echo "  j_space_vcf_pipeline fasta_file fastq_dir output_dir"
fi
genomeFa=$1
fastqDir=$2
basedir=$3

# PIPELINE FROM FASTQ to ANNOVAR
#basedir="/dati-raid/BimiB/share/aguidi/JOG-Space/ExperimentPipe/"


# REQUIRED EXTERNAL TOOLS:
# STAR 2.7.3a
# trimmomatic 0.39
# samtools 1.9
# gatk3 3.8-1-0-gf15c1c3ef
# picard 2.21.3
# annovar


# SETTING DIRECTORIES
#fastqDir=${basedir}FastaQ/
fastqTrimDir=${fastqDir}
singleSamDir=${basedir}singleSAM/
singleBamDir=${basedir}singleBAM/
picardDir=${basedir}picard/
vcfDir=${basedir}vcf/
#annovarDir=${basedir}annotated/

# SETTING REFERENCE FILES
#genomeFa=${basedir}Reference/reference.fasta
#annovar_db=${basedir}reference/annovar/RefGene_hg38/

###????###
##list of known driver in vcf format
#dbSNP_vcf=${basedir}reference/1000GENOMES-phase_3.vcf
echo $genomeFa
echo $fastqDir
echo $basedir 

# IMPORT FASTQs
fastqlist="$(find "$fastqDir" -type f  -exec basename {} \; | sort | cut -f 1 -d '.' | tr '\n' ' ')"

jobs=20


#===B0===
echo "trimming fastq"

mkdir -p $fastqTrimDir
if [ ! -d "$fastqTrimDir" ]; then
echo "Error mkdir"
exit 1
fi

# Default parameters. No illumina adapters
#parallel -j $jobs trimmomatic SE -phred33 -threads 3 -quiet ${fastqDir}{}.fastq.gz ${fastqTrimDir}{}.trim.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -summary ${fastqTrimDir}{}.trim.summary ::: $fastqlist

#===E0===



#===B1===
echo "bwa-mem2 -- generate genome indexing file"

bwa-mem2 index $genomeFa

#===E1===



#===B2===
echo "bwa-mem2 -- align fastq"

filelist="$(find "$fastqTrimDir" -type f -name "*.fq"| sort | tr '\n' ','| sed 's/,$//g')"
filelist="$(find "$fastqTrimDir" -type f -name "*.fq"| sort | tr '\n' ' '| sed 's/,$//g')"

mkdir -p $singleSamDir
if [ ! -d "$singleSamDir" ]; then
echo "Error mkdir"
exit 1
fi

cd $singleSamDir

#-R "@RG\t@ID:Lane1\tLB:Tumor\tPL:art\tPU:1\tSM:{1.}L001"
parallel  "bwa-mem2 mem -t 2 -M -R '@RG\tID:Lane1\tLB:Tumor\tPL:art\tSM:cell' -L 50 $genomeFa {1} {2} > ${singleSamDir}{1/.}.sam" ::: ${fastqDir}*_1.fq :::+ ${fastqDir}*_2.fq

cd $basedir

#===E2===



#===B=samtools===
echo "indexing fasta"

#indexing and dictionary creation of the fasta file inside 'reference' folder
samtools faidx $genomeFa
picard CreateSequenceDictionary R=$genomeFa

echo "The fasta file is indexed."
#===E=samtools===



#===B2=samtools===

#sort, mark duplicates, and create index
mkdir -p $singleBamDir
if [ ! -d "$singleBamDir" ]; then
echo "Error mkdir"
exit 1
fi

filelist="$(find "$singleSamDir" -type f -name "*.sam" -exec basename {} \; | sort | cut -f 1 -d '.' | tr '\n' ' ')"

#sam to bam
parallel -j 40 picard SamFormatConverter I=${singleSamDir}/{}.sam O=${singleBamDir}/{}Aligned.bam ::: $filelist

#gatk ClipReads   -I input_reads.bam   -O output_reads.bam

#sort
parallel -j 40 picard SortSam I=${singleBamDir}/{}Aligned.bam O="${singleBamDir}/{}Aligned.sortedByCoord.out.bam" SO=coordinate ::: $filelist

#getting bam index of the raw bam
parallel -j 40 picard BuildBamIndex I="${singleBamDir}/{}Aligned.sortedByCoord.out.bam" O="${singleBamDir}/{}Aligned.sortedByCoord.out.bai" ::: $filelist

#===E2=samtools===




#===B=Picard===

#Add read groups and mark duplicates
#The above step produces a SAM file, which we then put through the usual Picard processing steps: adding read group information, sorting, marking duplicates and indexing.

echo "add or replace read groups"

mkdir -p $picardDir
if [ ! -d "$picardDir" ]; then
echo "Error mkdir"
exit 1
fi

# GATK BEST PRACTICE
parallel -j 40 picard AddOrReplaceReadGroups I="${singleBamDir}{}Aligned.sortedByCoord.out.bam" O="${picardDir}{}.arrg_s1.bam" SO=coordinate RGID=0 RGPL=art RGLB=Tumor RGPU=group RGSM=cell ::: $filelist

parallel -j 20 picard MarkDuplicates I="${picardDir}{}.arrg_s1.bam" O="${picardDir}{}.md_s2.bam" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M="${picardDir}output_metrics_{}.txt" ::: $filelist

parallel -j 50 samtools index -b "${picardDir}{}.md_s2.bam" ::: $filelist

#parallel -j 20 gatk3 -Xmx32G -T BaseRecalibrator -R ${genomeFa} -I "${picardDir}{}.snc_s3.bam" -OQ -o "${picardDir}{}.table" -knownSites ${dbSNP_vcf} ::: $fastqlist

#parallel -j 40 gatk3 -Xmx3G -T PrintReads -kpr -R ${genomeFa} -I "${picardDir}{}.snc_s3.bam" -OQ -BQSR "${picardDir}{}.table" -o "${picardDir}{}.bbq_s4.bam" ::: $fastqlist

echo "PICARD DONE."

#===E=Picard===



#===B=HaplotypeCaller===
echo "HaplotypeCaller and Filtration"

mkdir -p $vcfDir
if [ ! -d "$vcfDir" ]; then
echo "Error mkdir"
exit 1
fi

#parallel -j 40 gatk3 -T HaplotypeCaller -R ${genomeFa} -I "${picardDir}{}.bbq_s4.bam" -o "${vcfDir}{}.raw.snps.indels.vcf" -dontUseSoftClippedBases -stand_call_conf 20 --dbsnp ${dbSNP_vcf} ::: $failelist

parallel -j 40 gatk3 -T HaplotypeCaller -R ${genomeFa} -I "${picardDir}{}.md_s2.bam" -o "${vcfDir}{}.raw.snps.indels.vcf" -dontUseSoftClippedBases -stand_call_conf 20 ::: $filelist

parallel -j 40 gatk3 -T VariantFiltration -R ${genomeFa} -V "${vcfDir}{}.raw.snps.indels.vcf" -window 35 -cluster 3 -filterName "FS" -filter \"FS '>' 30.0\" -filterName "QD" -filter \"QD '<' 2.0\" -o "${vcfDir}{}.filtered.vcf" ::: $filelist

#parallel -j 40 gatk3 -T VariantFiltration -R ${genomeFa} -V "${vcfDir}{}.raw.snps.indels.vcf" -filterName "BS" -filter \"FS '>' 60.0\" -o "${vcfDir}{}.filtered.vcf" ::: $filelist


#===E=HaplotypeCaller===

#===B=Annovar===

#echo "Annotation with annovar"
#
#mkdir -p $annovarDir
#if [ ! -d "$annovarDir" ]; then
# echo "Error mkdir"
# exit 1
#fi
#
#parallel -j 40 "${annovar_c2a_cmd}" -format vcf4 "${vcfDir}{}.filtered.vcf" --outfile "${annovarDir}{}.anninput" --includeinfo ::: $fastqlist
#parallel -j 40 "${annovar_av_cmd}" -out ${annovarDir}{} -exonicsplicing -build hg38 "${annovarDir}{}.anninput" "${annovar_db}" ::: $fastqlist
#
#===E=Annovar===


