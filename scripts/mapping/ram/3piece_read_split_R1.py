#!/usr/bin/python

#"usage: 3piece_read_split.py Input_file First_piece_size(bp) Last_piece_size(bp) upstream_to_trim(bp) downstream_to_trim(bp)"
import gzip
import sys
fn1=sys.argv[1]

if 'gz' in fn1:
        dfh1=gzip.open(fn1,'r')

if 'gz' not in fn1:
	dfh1=open(fn1,'r')

ID1=dfh1.readline().split()
line1=dfh1.readline().rstrip()
plus1=dfh1.readline().split()
QS1=dfh1.readline().rstrip()
rfhs=[]

size1=int(sys.argv[2])
size2=int(sys.argv[3])
trim1=3
trim2=5

for i in range(1,2):
	rfhs.append(open(fn1+'_r'+str(i)+'.fq','w'))


while line1:
  a = ['GATTGATT', 'GATTGATC']
  size1=40
  c=filter(lambda x:  x in line1, a)
  if c:
    size1=line1.find(c[0])
  c2=filter(lambda x:  x in line1[trim1+size1:len(line1)], a)
  size2=len(line1)
  if c2:
    size2=line1.find(c2[0])
  if len(line1[0:size1+trim1])>=30:
	  rfhs[0].write(ID1[0].split('_')[0]+'_'+ID1[0].split('_')[1].split(':')[0]+'-1:'+':'.join(ID1[0].split('_')[1].split(':')[1:])+'\n'+line1[0:size1+trim1]+'\n'+plus1[0]+'\n'+QS1[0:size1+trim1]+'\n')
  if len(line1[trim1+1+size1:size2+trim1])>=30:
	  rfhs[0].write(ID1[0].split('_')[0]+'_'+ID1[0].split('_')[1].split(':')[0]+'-2:'+':'.join(ID1[0].split('_')[1].split(':')[1:])+'\n'+'GATN'+line1[8+size1:size2+trim1]+'\n'+plus1[0]+'\n'+QS1[4+size1:size2+trim1]+'\n')
  if len(line1[4+size2:len(line1)])>=30:
	  rfhs[0].write(ID1[0].split('_')[0]+'_'+ID1[0].split('_')[1].split(':')[0]+'-3:'+':'.join(ID1[0].split('_')[1].split(':')[1:])+'\n'+'GATN'+line1[size2+8:len(line1)]+'\n'+plus1[0]+'\n'+QS1[size2+4:len(line1)]+'\n')
  ID1=dfh1.readline().split()
  line1=dfh1.readline().rstrip()
  plus1=dfh1.readline().split()
  QS1=dfh1.readline().rstrip()

map(lambda x:x.close(),rfhs)
