# Commands and parameters used for RNAseq data processing.
# Input/Output file names are omitted where possible for simplicity
# Basic workflow:
# 1) Illumina adapter removal
# 2) Map to T. spiralis genome
# 3) Find candidate read pairs for TSL screening
# 4) Screen candidate read pairs for TSL sequences
# 5) Map TSL read pairs to genome and count reads against genes/exons
# 6) Assemble de novo gene annotations from non-TSL read pairs

#
# 1) Illumina adapter removal
#
trim_galore --illumina --paired --retain_unpaired --stringency 3 --length 30 -q 20

#
# 2) Map to T. spiralis genome
#
# since the spliced leader sequence at the 5' end of a read will not map to 
# the genome, we enforce end-to-end mapping. Unmapped reads with mapped mates
# are then candidate reads for TSL screening.
hisat2 -x genome/T.spiralis -p 12 --no-softclip --no-discordant --dta \
	| samtools sort -O bam -@ 12 -o out.bam -

#
# 3) extract candidate read pairs for TSL sequence screening
#
# extract unmapped reads with mapped mates
# these reads need to be screened for TSL sequence
samtools view -b -h -f 5 -F 8 out.bam | \
	samtools sort -n -O bam -o out.bam.TSL.candidates -
samtools fastq -N -s out.bam.TSL.candidates.fastq.gz out.bam.TSL.candidates

# extract mapped reads with unmapped mates. These will be the mates for above.
samtools view -b -h -f 9 -F 260 out.bam | \
	samtools sort -n -O bam -o out.bam.TSL.candidates.mates -
samtools fastq -N -s out.bam.TSL.candidates.mates.fastq.gz out.bam.TSL.candidates.mates

#
# 4) Screen candidate read pairs for TSL sequences
#
# either with minimum 8bp or 10bp hits (-O) and error rates of 0.11 or 0.09 (-e) respectively
cutadapt -g file:TSL_leader_seqs_FWD.fasta -m 20 \
    -O 8 -e 0.11 \ 
    -o {name}.R1.fq.gz \
    -p {name}.R2.fq.gz \
    --untrimmed-output notags.R1.fq.gz \
    --untrimmed-paired-output notags.R2.fq.gz \
    out.bam.TSL.candidates.fastq.gz \
    out.bam.TSL.candidates.mates.fastq.gz
#					   
# 5) Map TSL read pairs back to genome and count reads against genes/exons
#
hisat2 -x genome/T.spiralis --no-softclip 
	| samtools sort -O bam -@ 8 -o TSL.bam -

# since we are only interested in the location of the 5' trans-spliced read, 
# we need to extract R1 reads from the BAM files.
# counting PE fragments would lead to too much ambiguity due to multiple overlapping.
samtools view -h -O bam -f 64 TSL.bam > TSL.bam.R1

# quantify genes
featureCounts -a genome/trichinella_spiralis.PRJNA12603.WBPS10.annotations.gff3.gtf \
				-t exon -g gene_id -s 0 --largestOverlap 
				-o TSL.bam.R1.FC.gene TSL.bam.R1
# quantify exons
featureCounts -a genome/trichinella_spiralis.PRJNA12603.WBPS10.annotations.gff3  \
				-f -t exon -g ID -s 0 --largestOverlap  \
				-o TSL.bam.R1.FC.exon TSL.bam.R1
#	
# 6) Assemble de novo gene annotations from non-TSL read pairs
#
# Trinity assembly
Trinity --genome_guided_bam out.bam \
		--genome_guided_max_intron 1000000 \
		--no_normalize_reads --full_cleanup
TransDecoder.LongOrfs -m 50 -t Trinity-GG.fasta --gene_trans_map Trinity-GG.fasta.gene_trans_map
TransDecoder.Predict --single_best_only -t Trinity-GG.fasta
cd-hit -c 0.99 -i Trinity-GG.fasta.transdecoder.pep
cd-hit-est -c 0.99 -i Trinity-GG.fasta.transdecoder.cds

# BRAKER annotation
RepeatMasker -species nematoda -xmall 
braker.pl --genome=genome.fa.masked --softmasking=on \
	--bam=out.bam --prot_seq=Trinity-GG.fasta.transdecoder.pep.cdhit99 \
	--gth2traingenes --gff3
braker.pl --genome=genome.fa.masked --softmasking=on \
	--bam=out.bam --gff3
	
# STRINGTIE annotation
stringtie out.bam -l Tsp_denovo -m 50 
# problem: mono-exonic genes have no strand information
# can infer strand by predicting ORFs from these transcripts
grep -e "\stranscript\s.*\s\.\s\.\s" stringtie.gtf \
	| gffread - -g genome.fa -w stringtie.nostrand.fasta
TransDecoder.LongOrfs -m 1 -t stringtie.nostrand.fasta
TransDecoder.Predict -t stringtie.nostrand.fasta --single_best_only --retain_long_orfs_mode strict --retain_long_orfs_length 10
sedpatterns=""
while read gene strand
do
	gene=${gene//\./\\.}
	sedpatterns+="-e '/${gene}/ s/\.\s\./${strand}\t\./' "
done < <(grep "gene" transdecoder.gff3 | cut -f 1,7)
eval "sed $sedpatterns stringtie.gtf > stringtie.gtf.strands"

# generate a merged annotation set across two BRAKER runs and STRINGTIE
gffcompare -T -C -X \
	-o tmp.merged \
	-s genome.fa.masked \
	stringtie.gtf.strands \
	braker1/augustus.hints.gtf.fixed \
	braker2/augustus.hints.gtf.fixed

# tidy up all GFF/GTF files for use with featureCounts
# also extract transcript sequences from genome based on annotations
gffcompare -T -C -X -s genome.fa.masked 
gffread -g genome.fa.masked \
	-o gtf.clean \
	-w gtf.clean.transcripts.fasta \
	-M -Q -T
bedtools sort -i gtf.clean > gtf.clean.sorted
	
# extract unique exons for featureCounts 
grep "\sexon\s" gtf.clean.sorted \
	| bedtools sort -i stdin \
	| bedtools merge -s -c 4,5,9 -o min,max,first -i stdin \
	| awk -v prog="BRAKER" -v OFS='\t' -F'\t' '{print $1, prog, "exon", $5, $6, ".", $4, ".", $7}' \
	> gtf.clean.sorted.unique_exons

# The final GTF files can then be used with featureCounts as above, using stranded mode (-s 1)