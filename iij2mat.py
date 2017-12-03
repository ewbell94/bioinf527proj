from sys import argv

f=open(argv[1])

AAseq=""
scoremat=[]
for line in f:
    if line.strip()[0]=="#":
        continue
    parts=line.strip().split()
    if len(parts)==0:
        continue
    try:
        scorerow=[int(i) for i in parts]
        scoremat.append(scorerow)
    except ValueError:
        AAseq="".join(parts)
f.close()

for i in range(len(scoremat)):
    for j in range(i+1,len(scoremat)):
        scoremat[i].append(scoremat[j][i])

matname=argv[1].split("/")[-1].split(".")[0]

g=open(matname+".mat","w")
for i in AAseq:
    g.write("\t"+i)
g.write("\n")
for i in range(len(AAseq)):
    g.write(AAseq[i])
    for num in scoremat[i]:
        g.write("\t"+str(num))
    g.write("\n")
g.close()