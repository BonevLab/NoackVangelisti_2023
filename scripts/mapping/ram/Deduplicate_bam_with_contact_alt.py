#### Modified from TAURUS Lee et al., Nat Methods 2019 (https://github.com/dixonlab/Taurus-MH)

import pysam  
import sys
import os
os.system('sort -k1n '+sys.argv[1]+'_multi_split_aligned.txt_man_dedupped.txt > '+sys.argv[1]+'_multi_split_aligned.txt_man_dedupped_queryname_sorted.txt')
bam=sys.argv[1]
dfh=pysam.AlignmentFile(bam, "rb")
dfh2=open(sys.argv[1]+'_multi_split_aligned.txt_man_dedupped_queryname_sorted.txt','r')
rfh=pysam.AlignmentFile(sys.argv[1].replace('.bam','_dedupped.bam'), "wb", template=dfh)

ids=[]

for i in dfh2:
    line=i.split()
    id=':'.join(line[1].split(':')[4:7])
    ids.append(id)
ids=set(ids)
for a in dfh:
    read=str(a).split()
    id=read[0].split('_')[0]
    if ':'.join(id.split(':')[4:7]) in ids:
		rfh.write(a)
