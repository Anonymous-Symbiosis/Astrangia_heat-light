# Tag-based RNA-seq reads processing pipeline

# More info on the packages
# fastx_toolkit: http://hannonlab.cshl.edu/fastx_toolkit/download.html
# bowtie2: http://bowtie-bio.sourceforge.net/index.shtml 
# pipeline came from https://github.com/z0on/tag-based_RNAseq

#-----Grab data and scripts-----
# extract and copy  fastq.gz files from multiple subdirectory to a single folder
# in my case, from multiple folders within Astrangia_MPCC2022 to the working directory
find /projectnb/coral/user/Astrangia_MPCC2022/BaseSpace/ -name "*.fastq.gz" -exec cp {} /projectnb/coral/user/Astrangia_MPCC2022/job_May2025/ \;

# unzip gz files
gunzip *.gz

# check data
cd /projectnb/coral/user/Astrangia_MPCC2022/job_May2025
ls -l

# count number of files
ls -1| wc -l

# concatenate host and symbiont genomes (from Dan)
cat ../reference_files/apoculata.genome.fasta ../reference_files/algae_reftranscriptome.fasta > holobiont.fasta

# count the number of contigs in your fasta file
grep '>' holobiont.fasta | wc -l 

# concatenate host and symbiont annotation files (from Dan)
cat ../reference_files/apoculata.gff3 algae_ref.gff > holobiont.gff3

#-----Trimming-----
# adaptor trimming, deduplicating, and quality filtering
# load module fastx-toolkit
module load fastx-toolkit

# create and launch the cleaning process for all files in the same time:
tagseq_trim_launch.pl '\.fastq$' > clean 

# check the folder clean
nano clean

# create a batch job and submit to execute all commands written to file ‘clean’
scc6_qsub_launcher.py -N trim -P coral -jobsfile clean 
qsub trim_array.qsub
## these commands (tagseq_clipper.pl) are clipping poly A tails, Illumina adapters, short reads (<20bp) and low quality (to keep 90% of a read must be quality score >20)

# check status
qstat -u user

# when the job is done, have a look in the trim.e* file
cat trim.(tab complete)
##this has all of the info for trimming. You'll see many sequences are PCR duplicates because this is TagSeq data and remember that we incorporated the degenerate bases into the cDNA synthesis. 

mv trim* parallel_outputs

#-----Mapping-----
# create bowtie2 index for reference transcriptome
module load bowtie2
bowtie2-build holobiont.fasta holobiont.fasta 

# create and launch the mapping process for all files in the same time:
for file in *.trim
	do echo "bowtie2 -x holobiont.fasta -U $file --local -p 4 -S ${file/.trim/}.sam">> maps
done
nano maps

# create a batch job and submit to execute all commands written to file ‘maps’
scc6_qsub_launcher.py -N maps -P coral -jobsfile maps
qsub maps_array.qsub

nano maps.o(tab complete)
mv maps* parallel_outputs

#-----Sort the sam alignment files and convert to bams-----
module load samtools

for file in *.sam
	do echo "samtools sort -n -O bam -o ${file/.sam/}.bam $file" >> sortConvert
done 

scc6_qsub_launcher.py -N sortconverting -P coral -jobsfile sortConvert
qsub sortconverting_array.qsub

mv sortconverting* parallel_outputs

#-----Counting-----

module load subread

# Choose the GFF
MY_GFF="holobiont.gff3"; GENE_ID="ID"
featureCounts -a $MY_GFF -t gene -g $GENE_ID -o feature_counts_out.txt -T 32 *.bam

# get alignment counts
for file in *.bam
	do echo "samtools flagstat $file > ${file/.bam/}_flagStats.txt" >> getInitialAlignment
done

module load samtools
scc6_qsub_launcher.py -N samtools_alignment_counts -P coral -jobsfile getInitialAlignment
qsub samtools_alignment_counts_array.qsub

mv samtools* parallel_outputs

#-----Making raw read count table-----
module purge
module load python2
module load htseq/0.11.0

MY_GFF="holobiont.gff3"; GENE_ID="ID"

for file in *.sam
	do echo "htseq-count -t gene -i ID -m intersection-nonempty --stranded=no $file $MY_GFF > ${file/sam/counts.txt}" >> doCounts
done

scc6_qsub_launcher.py -N counts -P coral -jobsfile doCounts
qsub counts_array.qsub 

#-----Compile to make big expression table-----

expression_compiler.pl *.counts.txt > holobiont_counts.txt


#-----Making a table compiling the number of counts in each step-----

>mapped_count.tsv
for file in *_flagStats.txt
do pp=$(grep "mapped" $file | head -n 1)
 echo -e "$file\t$pp" |\
 awk '{split($1, a, "_flagStats.txt")
 print a[1]"\t"$2"\tmapped"}' >> mapped_count.tsv
 done

# raw
wc -l *.fastq |\
	awk '{split($2, a, ".fastq")
	print a[1]"\t"$1/4"\rawCounts"}' |\
grep -v total > raw_read_counts.tsv 

# trimmed
wc -l *.trim |\
	awk '{split($2, a, ".trim")
	print a[1]"\t"$1/4"\ttrimmedCounts"}' |\
grep -v total > trimmed_read_counts.tsv 

# Compile pipeline counts for supps
cat *.tsv > final_pipeline_counts.txt
