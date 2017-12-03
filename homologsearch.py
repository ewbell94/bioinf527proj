from copy import deepcopy
import subprocess
from sys import argv
import random
from multiprocessing.dummy import Pool as ThreadPool
import time

Nprots=1000
mats=["randblos.mat","blosum62.mat","new99.mat"]

f=open("astral-scopedom-seqres-gd-sel-gs-bib-40-2.06.fa")
prots=[]
seq=""
currprot=["",""]

for line in f:
    if line[0]==">":
        if seq=="":
            scopid=line.split()[1]
            parts=scopid.split(".")
            currprot[0]=".".join(parts[:3])
        else:
            currprot[1]=seq
            prots.append(deepcopy(currprot))
            seq=""
            scopid=line.split()[1]
            parts=scopid.split(".")
            currprot[0]=".".join(parts[:3])
    else:
        seq+=line.strip().upper()

currprot[1]=seq
prots.append(currprot)

#Nprots=len(prots)
currmat="Default"

def align(i):
    conttab=[[0,0],[0,0]]
    for j in range(i):
        if i!=j:
            if len(prots[i][1])>len(prots[j][1]):
                seq1=prots[i][1]
                seq2=prots[j][1]
            else:
                seq1=prots[j][1]
                seq2=prots[i][1]
            seq1=prots[i][1]
            seq2=prots[j][1]
            process=subprocess.Popen(["./align","-p",seq1,"-q",seq2,"-m",currmat,"-Q"],stdout=subprocess.PIPE)
            seqid,err=process.communicate()
            seqid=float(seqid)
            if seqid>0.25:
                if prots[i][0]==prots[j][0]:
                    conttab[0][0]+=1
                else:
                    conttab[0][1]+=1
            else:
                if prots[i][0]==prots[j][0]:
                    conttab[1][0]+=1
                else:
                    conttab[1][1]+=1
    return conttab

if len(argv) > 1:
    random.seed(argv[1])
random.shuffle(prots)

for mat in mats:
    print (time.strftime("%m/%d/%Y %H:%M:%S"))
    currmat=mat
    pool=ThreadPool(8)
    eyes=[i for i in range(Nprots)]
    res=pool.map(align,eyes)
    pool.close()
    pool.join()
    conttab=[[0,0],[0,0]]
    for i in res:
        conttab[0][0]+=i[0][0]
        conttab[0][1]+=i[0][1]
        conttab[1][0]+=i[1][0]
        conttab[1][1]+=i[1][1]
    print(conttab)
    sumall=float(sum([conttab[0][0],conttab[0][1],conttab[1][0],conttab[1][1]]))
    print(sum([conttab[0][0],conttab[1][1]])/sumall)
    print(float(conttab[0][0])/sum([conttab[0][0],conttab[0][1]]))
    print(float(conttab[0][0])/sum([conttab[0][0],conttab[1][0]]))
        
print (time.strftime("%m/%d/%Y %H:%M:%S"))
