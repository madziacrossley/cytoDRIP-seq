#cytoDRIP-seq alignment from Crossley et al., 2024, Nature PMID: 36544021
#Raw fastq files available from GEO under accession number GSE178841.

#create sample_names.txt
#CON_1_S9p6
#CON_2_S9p6
#SETX_1_S9p6
#SETX_2_S9p6
#CON_1_IGG
#CON_2_IGG
#SETX_1_IGG
#SETX_2_IGG

#give everyone read only permission on raw gzip 
#raw_files/chmod 444 .fq.gz*

##quality control check with fastqc

mkdir fastqc
find *.fastq.gz | grep -v 'Und' | xargs fastqc -o ../fastqc/ --threads 10

#Look at html files generated:
#Samples 1-4 (i.e. S9.6 IPs) look good, good sequence quality and other scores. As expected for R2s see G tail. For IgG samples (5-8) per base sequence content look lower and seems to be more duplication. Overall IgGs have lower diversity libraries than IPs.
#Note from 1S kit on FASTQC: https://swiftbiosci.com/wp-content/uploads/2019/02/16-0853-Tail-Trim-Final-442019.pdf
#"Please note: Quality control software, such as FastQC FastQC (Babraham Bioinformatics) may raise “Per base sequence content” or “Per base GC content” flags at the beginning of R2"
#Therefore, trim off G-rich tail prior to alignment

#trimming adapters
#trim off adapters and secondly G-rich tails (order is important, G-rich tails are added with the ssDNA library kit. This improves alignment)

mkdir adapter_tail_trim
#Using P5 TruSeq LT adapter (universal sequence)
#Actual sequences
#5'-3' AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
#5'-3' GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNCTCGTATGCCGTCTTCTGCTTG
#with Ns representing the unique barcode

ADAPTER_FWD=GATCGGAAGAGCACACGTCTGAACTCCAGTCACN{8}
ADAPTER_REV=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

cat sample_names.txt | xargs -I {} -P10 sh -c 'cutadapt -j 10 --minimum-length 30 -a $ADAPTER_FWD -A ADAPTER_REV -o a.{}.R1.fq.gz -p a.{}.R2.fq.gz {}.R1.fq.gz {}.R2.fq.gz 1> a.{}.txt'
cat sample_names.txt | xargs -I {} -P10 sh -c 'cutadapt -j 10 --minimum-length 30 -a "C{30}" -U 15 -o tail_trim.{}.R1.fq.gz -p tail_trim.{}.R2.fq.gz a.{}.R1.fq.gz a.{}.R2.fq.gz 1> tail.{}.txt'
#check cutadapt reports .txt files to make sure trimming has worked. Can also repeat fastqc to see if per base GC content flag at beginning of R2 is now gone

##Update 2023: since original experiments were done IDT acquired Swift Biosciences. The adapter types are now different
#Using i5 and i7 TruSeq HT adapters
#Actual adapter sequences
#Index 1 (i7) Adapters
#5´–GATCGGAAGAGCACACGTCTGAACTCCAGTCACXXXXXXXXATCTCGTATGCCGTCTTCTGCTTG–3´
#Index 2 (i5) Adapters
#5´–AATGATACGGCGACCACCGAGATCTACACYYYYYYYYACACTCTTTCCCTACACGACGCTCTTCCGATCT–3´
#With X and Y standing for the unique barcodes

#Use Ns in place of 8 base unique barcode for each i5 and i7
#ADAPTER_FWD = GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
#ADAPTER_REV = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT
##End of update

#remove adapter trim files to keep adapter and tail trimmed fq files
rm a.*.R1.fq.gz
rm a.*.R2.fq.gz

#bowtie2 alignment
#cytoDRIP-seq - aligned to hg38 (end to end, very sensitive settings), 
mkdir -p ../bowtie_e2e_alignment

INDEX="/Genomics/indexes/hg38/hg38"
cat sample_names.txt | xargs -I {} -P10 sh -c 'bowtie2 -p10 --no-unal --no-mixed --no-discordant -x $INDEX -1 adapter_tail_trim/tail_trim.{}.R1.fq.gz -2 adapter_tail_trim/tail_trim.{}.R2.fq.gz | samtools view -bSh > bowtie_e2e_alignment/{}.bam'
#run script
sh ./cytodrip_align.sh 2> bowtie_e2e_alignment/bowtie_report.txt
#alignment rates: for S9.6 IPs 84-86%, for IgG 42-55%


#Split reads into plus and minus strands
#CytoDRIP-seq is stranded, since ssDNA of RNA-DNA hybrid is captured in library prep
OUTDIR="aligned_split_reads"
SPLIT_FORWARD='(NF < 15 || (NF > 15 && ($2 == 99 || $2 == 147))){print}'
SPLIT_REVERSE='(NF < 15 || (NF > 15 && ($2 == 83 || $2 == 163))){print}'

mkdir -p $OUTDIR

for FILENAME in *.bam; do
    echo "echo 'Writing ${OUTDIR}/p_${FILENAME}...'"
    echo "samtools view -h $FILENAME | awk '${SPLIT_FORWARD}' | samtools view -bSh > ${OUTDIR}/p_$FILENAME"
    echo "echo 'Writing ${OUTDIR}/m_${FILENAME}...'"
    echo "samtools view -h $FILENAME | awk '${SPLIT_REVERSE}' | samtools view -bSh > ${OUTDIR}/m_$FILENAME"
done

#run in aligned_split_reads to generate the script file 
sh ../../bowtie_e2e_alignment/cytodrip_align.sh > run_read_split.sh
sh ./run_read_split.sh

#Use bedtools utilities to generate genome tracks from split m and p bam files

GENOMEFILE="/shared_files/util/hg38.chrom.sizes"
TOBIGWIG="/shared_files/util/bedGraphToBigWig"
BEDTOOLSPATH="bedpe"
CUTPATH_HUMAN="bed"
UNIQ="rmdup"
TRACKS="tracks"
BIGWIG="bigwig"

mkdir -p $BEDTOOLSPATH $CUTPATH_HUMAN $SORTED $UNIQ $TRACKS $BIGWIG

for FILENAME in *.bam; do
    echo "echo \"Processing file $FILENAME ...\""
    echo "bedtools bamtobed -bedpe -i $FILENAME > ${BEDTOOLSPATH}/${FILENAME}.bedpe" # Convert from BAM to BEDPE
    echo "cut -f 1,2,6 ${BEDTOOLSPATH}/${FILENAME}.bedpe > ${CUTPATH_HUMAN}/${FILENAME}.bed" # Take columns 1 (chr R1), 2 (start R1), and 6 (End R2)
    echo "gsort --parallel=10 -u -k1,1 -k2,2n -k3,3n ${CUTPATH_HUMAN}/${FILENAME}.bed > ${UNIQ}/${FILENAME}.bed" # Sort according to chromosome, then start, then end, then emove duplicate fragments from the file, parallelized sort over 10 cores
    echo "bedtools genomecov -bg -g $GENOMEFILE -i ${UNIQ}/${FILENAME}.bed > ${TRACKS}/${FILENAME}.bedGraph" #Make bedGraph
    echo "$TOBIGWIG ${TRACKS}/${FILENAME}.bedGraph $GENOMEFILE ${BIGWIG}/${FILENAME}.bw"  # Convert to bigWig format, note this bw is not read count normalized
done

#run in /split_aligned_reads
sh ../bowtie_e2e_alignment/cytodrip_align.sh > run_bedtools.sh
sh ./run_bedtools.sh

#Count the reads for each p and m file
rm -ri /bed #don't need these files, leave bed files in /rmdup
#in /rmdup
wc -l *.bed > read_count.txt
#Add up read counts for m_ and p_ bed files for each sample, then divide by 10e3 to get unique reads per million mapped

#Read-count normalization (RPM)
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,-$4/8.495319)}' m.CON_1_IGG.bedGraph > CON_1_IGG.m.bedGraph.norm
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/8.495319)}' p.CON_1_IGG.bedGraph > CON_1_IGG.p.bedGraph.norm
#repeat for all files, using RPM normalization factors 
#Normalization factors
#CON_1_S9p6 = 28.867118
#CON_2_S9p6 = 34.340676
#SETX_1_S9p6 = 33.874164
#SETX_2_S9p6 = 41.024399
#CON_1_IGG = 8.495319
#CON_2_IGG = 9.398247
#SETX_1_IGG = 2.537287
#SETX_2_IGG = 4.623568

#Make RPM normalized bigWig files
find *.norm | sed 's/\.bedGraph.norm//' | xargs -I {} -P 10 /shared_files/util/bedGraphToBigWig {}.bedGraph.norm /shared_files/util/hg38.chrom.sizes {}.bw

