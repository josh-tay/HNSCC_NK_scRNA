
# runs STAR aligner for RNAseq fastqs, maps to hg19 ercc92 concatenated genomes first

import sys  # import the necessary modules
import os
import os.path
import glob
import re
import time

sample_fastq = sys.argv[1]  # get file with sample names and fastq

datetime = time.strftime("%Y%m%d-%H%M")

STAR_executable = "/media/data/Josh/packages/STAR-2.7.0f/bin/Linux_x86_64/STAR"
STAR_hg19genomeDir = "/media/data/Josh/NPC_star_indices/hg19ercc92_starindex"
OUTPUT_Root = "/media/data/Josh/NK_scRNA"
OUTPUT_Dir = OUTPUT_Root + "/" + datetime + "_bam"

make_outputdir_command = "mkdir " + OUTPUT_Dir
os.system(make_outputdir_command)

# load genome references for STAR

print("Loading STAR genome references")
load_genome_command = STAR_executable + " --genomeLoad LoadAndExit --genomeDir " + STAR_hg19genomeDir
print(load_genome_command)
os.system(load_genome_command)


with open(sample_fastq, "r") as sample_fastq_file:

    for line in sample_fastq_file:
        line = line.strip()
        line_list = line.split('\t')
        sample_no = line_list[0]
        fastq_1 = line_list[1]
        fastq_2 = line_list[2]

        print(sample_no, fastq_1, fastq_2)

        STAR_cmd_str = STAR_executable + \
            " --runThreadN 36" + \
            " --genomeLoad LoadAndKeep" + \
            " --genomeDir " + STAR_hg19genomeDir + \
            " --limitBAMsortRAM 40147483648" + \
            " --readFilesCommand zcat" + \
            " --readFilesIn " + fastq_1 + " " + fastq_2 + \
            " --outFileNamePrefix " + OUTPUT_Dir + "/" + sample_no + \
            " --outSAMtype BAM SortedByCoordinate" + \
            " --quantMode GeneCounts" + \
            " --outFilterType BySJout" + \
            " --outFilterMultimapNmax 20" + \
            " --alignSJoverhangMin 8" + \
            " --alignSJDBoverhangMin 1" + \
            " --outFilterMismatchNmax 999" + \
            " --alignIntronMin 20" + \
            " --alignIntronMax 1000000" + \
            " --alignMatesGapMax 1000000" + \
            " --outSAMattributes NH HI AS NM MD"

            # " --outFilterMismatchNoverReadLmax 0.04" + \
            # " --outFilterScoreMinOverLread 0.85" + \
            # " --outFilterIntronMotifs RemoveNoncanonicalUnannotated" + \
            # " --outReadsUnmapped Fastx" + \
            # " --outSAMunmapped Within" + \
            # no clipping performed, apart from clip3pAdapterMMp default: 0.1

        print(STAR_cmd_str)
        os.system(STAR_cmd_str)

print("Unloading STAR genome references")
unload_genome_command = STAR_executable + " --genomeLoad Remove --genomeDir " + STAR_hg19genomeDir
os.system(unload_genome_command)

print("Creating BAM indexes")
bam_index_command = "ls " + OUTPUT_Dir + "/*.bam | parallel samtools index '{}'"
print(bam_index_command)
os.system(bam_index_command)







