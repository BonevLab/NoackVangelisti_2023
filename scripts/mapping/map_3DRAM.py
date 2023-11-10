#!/usr/bin/env python

#### Modified from TAURUS Lee et al., Nat Methods 2019 (https://github.com/dixonlab/Taurus-MH)

import os
import sys
com=sys.argv
if '-1' not in com or "-2" not in com:
	print ("USAGE: python map_3DRAM.py -r <genome_folder> -1 <G to A converted mate> -2 <C to T converted mate> -Trim1 <Trim off first Xbp of R1> -Trim2 <Trim off first Xbp of R2> -split1 <Use first Xbp of unmapped read for 1st split reads> -split2 <Use last Xbp of unmapped read for 2nd split reads> -Threads <number of threads> -cutoff <CpG coverage cutoff to be considered for analysis> -genome <genome to use, either mm10 or hg38>")
picard="/home/hpc/bonev/software/Taurus-MH"
bismark="/home/hpc/common/bin"
bowtie="/home/hpc/common/bin"
python="/home/hpc/bonev/software/miniconda3/envs/taurus/bin"
TAURUS_loc="scripts/mapping"
juicer_loc="/home/hpc/bonev/juicer"
trim_loc="/home/hpc/bonev/software/miniconda3/envs/taurus/bin"

R1_pre=sys.argv[com.index('-1')+1]
R2_pre=sys.argv[com.index('-2')+1]
if '-r' in com:
  REF=sys.argv[com.index('-r')+1]
if '-r' not in com:
  REF='/home/hpc/bonev/annotations/mm10/'
if '-Trim1' in com:
	T1=sys.argv[com.index('-Trim1')+1]
if '-Trim1' not in com:
	T1='1'
if '-Trim2' in com:
	T2=sys.argv[com.index('-Trim2')+1]
if '-Trim2' not in com:
	T2='15'
if '-split1' in com:
	S1=sys.argv[com.index('-split1')+1]
if '-split1' not in com:
	S1='40'
if '-split2' in com:
	S2=sys.argv[com.index('-split2')+1]
if '-split2' not in com:
	S2='40'
if '-Threads' in com:
	Threads=sys.argv[com.index('-Threads')+1]
if '-Threads' not in com:
	Threads='4'
if '-cutoff' in com:
	cutoff=sys.argv[com.index('-cutoff')+1]
if '-cutoff' not in com:
	cutoff='6'
if '-genome' in com:
	Genome=sys.argv[com.index('-genome')+1]
if '-genome' not in com:
  Genome='mm10'
cutsite='GATC'
	
R1=R1_pre.split('/')[-1]+'_trimmed.fastq'
R2=R2_pre.split('/')[-1]+'_trimmed.fastq'
preTrim_R1=R1.split('.f')[0]+'_val_1.fq'
preTrim_R2=R2.split('.f')[0]+'_val_2.fq'

if 'gz' in R1.split('.')[-1]:
	R1_mod='.'.join(R1.split('.')[:-2]).split('/')[-1]
if 'gz' not in R1.split('.')[-1]:
	R1_mod='.'.join(R1.split('.')[:-1]).split('/')[-1]
if 'gz' in R2.split('.')[-1]:
        R2_mod='.'.join(R2.split('.')[:-2]).split('/')[-1]
if 'gz' not in R2.split('.')[-1]:
        R2_mod='.'.join(R2.split('.')[:-1]).split('/')[-1]

rfh=open('Run_TAURUS-MH_'+R1.split('/')[-1]+'.sh','w')

rfh.write("bismark="+bismark\
+"\nbowtie="+bowtie\
+"\npicard="+picard\
+"\nTAURUS_loc="+TAURUS_loc\
+"\npython="+python\
+"\nREF="+REF\
+"\nR1_pre="+R1_pre\
+"\nR2_pre="+R2_pre\
+"\nR1="+R1\
+"\nR2="+R2\
+"\nR1_mod="+R1_mod\
+"\nR2_mod="+R2_mod\
+"\n"
+"\nmkdir ./tmp"
+"\n${trim_loc}/trim_galore --nextseq 30 --dont_gzip -j "+Threads+" --clip_R1 "+T1+" --clip_R2 "+T2+" --length 20 --paired "+R1_pre+" "+R2_pre+" &"\
+"\nwait"\
+"\nmv "+preTrim_R1+" "+R1\
+"\nmv "+preTrim_R2+" "+R2\
+"\n"
+"\n${bismark}/bismark --parallel "+Threads+" -un --bowtie2 ${REF} ${R1} & "
+"\n${bismark}/bismark --parallel "+Threads+" -un --pbat --bowtie2 ${REF} ${R2} &"
+"\nwait"\
+"\n"\
+"\n${python}/python ${TAURUS_loc}/ram/3piece_read_split_R1.py ${R1}_unmapped_reads.fq.gz "+S1+" "+S2+" &"\
+"\n${python}/python ${TAURUS_loc}/ram/3piece_read_split_R2.py ${R2}_unmapped_reads.fq.gz "+S1+" "+S2+" &"\
+"\nwait"\
+"\n"
+"\n${bismark}/bismark --parallel "+Threads+" --bowtie2 ${REF} ${R1}_unmapped_reads.fq.gz_r1.fq &"\
+"\n${bismark}/bismark --parallel "+Threads+" --pbat --bowtie2 ${REF} ${R2}_unmapped_reads.fq.gz_r1.fq &"\
+"\nwait"\
+"\n"
+"\nrm ${R1}_unmapped_reads.fq.gz ${R2}_unmapped_reads.fq.gz"\
+"\nrm ${R1}_unmapped_reads.fq.gz_r1.fq ${R2}_unmapped_reads.fq.gz_r1.fq"\
+"\nrm ${R1} ${R2}"\
+"\n"
+"\njava -jar -Xmx10g ${picard}/picard.jar MergeSamFiles VERBOSITY=ERROR TMP_DIR=./tmp MAX_RECORDS_IN_RAM=2000000 SO=queryname \\"
+"\nI=${R1_mod}_bismark_bt2.bam \\"
+"\nI=${R2_mod}_bismark_bt2.bam \\"
+"\nI=${R1}_unmapped_reads.fq.gz_r1_bismark_bt2.bam \\"
+"\nI=${R2}_unmapped_reads.fq.gz_r1_bismark_bt2.bam \\"
+"\nO=${R1}_all_merged_3split.bam "\
+"\n"
+"\nrm ${R1_mod}_bismark_bt2.bam ${R2_mod}_bismark_bt2.bam ${R1}_unmapped_reads.fq.gz_r1_bismark_bt2.bam ${R2}_unmapped_reads.fq.gz_r1_bismark_bt2.bam"\
+"\n"
+"\n${python}/python ${TAURUS_loc}/ram/Bam_to_multi_contact2.py ${R1}_all_merged_3split.bam"\
+"\n${python}/python ${TAURUS_loc}/ram/Deduplicate_multi_contact.py ${R1}_all_merged_3split.bam_multi_split_aligned.txt"\
+"\n${python}/python ${TAURUS_loc}/ram/Multi_contact_to_two_contact_stat.py ${R1}_all_merged_3split.bam_multi_split_aligned.txt_man_dedupped.txt"\
+"\n${python}/python ${TAURUS_loc}/ram/Deduplicate_bam_with_contact_alt.py ${R1}_all_merged_3split.bam"\
+"\n${python}/python ${TAURUS_loc}/ram/Stat_summary.py ${R1_pre}"\
+"\nperl ${TAURUS_loc}/ram/fragment_4dnpairs.pl ${R1}_all_merged_3split.bam_multi_split_aligned.txt_man_dedupped.txt_2_contacts.txt "+preTrim_R1.split('_R1')[0]+"_frag.txt /home/hpc/bonev/juicer/restriction_sites/"+Genome+"_DpnII.txt"\
+"\nsed -i s/\-/16/g "+preTrim_R1.split('_R1')[0]+"_frag.txt"\
+"\nsed -i s/\+/0/g "+preTrim_R1.split('_R1')[0]+"_frag.txt"\
+"\nperl ${juicer_loc}/scripts/common/statistics.pl -s ${juicer_loc}/restriction_sites/"+Genome+"_DpnII.txt -o summary2.txt "+preTrim_R1.split('_R1')[0]+"_frag.txt"\
+'''\nawk 'NR==2 {printf"Total sequenced:\\t%s\\nMapped reads:\\t%s (%s %)\\nUnique:\\t%s (%s %)\\nCpG Methylation:\\t%s\\n",$2,$4,$4/$2*100,$5,$5/$4,$12}' summary.txt > summary1.txt'''\
+"\ncat summary1.txt summary2.txt > summary.txt"\
+"\nrm summary1.txt summary2.txt"\
+"\nrm *_SE_report.txt"\
+"\nrm ${R1}_all_merged_3split.bam"\
+"\nbismark_methylation_extractor --parallel "+Threads+" --ample_memory -s --bedGraph --ignore 6 --CX --genome_folder ${REF} ${R1}_all_merged_3split_dedupped.bam"\
+"\ncoverage2cytosine --genome_folder ${REF} --nome-seq -o ${R1} ${R1}_all_merged_3split_dedupped.bismark.cov.gz"\
+'''\nawk '($5+$6)>=1{printf"%s\\t%s\\t%s\\t%.0f\\n",$1,$2-1,$2,$4}' ${R1}.NOMe.CpG.cov | sort -k1,1 -k2,2n > ${R1}_CpG_cov1x.bedGraph'''\
+'''\nawk '($5+$6)>=1{printf"%s\\t%s\\t%s\\t%.0f\\n",$1,$2-1,$2,$4}' ${R1}.NOMe.GpC.cov | sort -k1,1 -k2,2n > ${R1}_GpC_cov1x.bedGraph'''\
+"\nbedGraphToBigWig ${R1}_CpG_cov1x.bedGraph ${REF}/"+Genome+".chrom.sizes "+preTrim_R1.split('_R1')[0]+"_CpG_cov1x.bw"\
+"\nbedGraphToBigWig ${R1}_GpC_cov1x.bedGraph ${REF}/"+Genome+".chrom.sizes "+preTrim_R1.split('_R1')[0]+"_GpC_cov1x.bw"
+"\nrm *_all_merged_3split_dedupped.txt *_man_*"
)

rfh.close()
#os.system("nohup sh Run_TAURUS-MH_"+R1.split("/")[-1]+".sh 2>log.txt &")
