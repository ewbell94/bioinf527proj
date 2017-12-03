import random

aastring="ARNDCQEGHILKMFPSTWYVBZX*"
mat=[]

def addnum(numcounts,num):
    try:
        val=numcounts[num]
        numcounts[num]=val+1
    except:
        numcounts[num]=1
    return numcounts

def normalize(numcounts):
    tot=0.
    for i in numcounts.keys():
        tot+=numcounts[i]
    for i in numcounts.keys():
        numcounts[i]=numcounts[i]/tot
    return numcounts

def drawnum(numcounts):
    r=random.random()
    thresh=0.
    for i in numcounts.keys():
        thresh+=numcounts[i]
        if r < thresh:
            return i

offcounts={}
matchcounts={}
f=open("blosum62.mat")
f.readline()
lines=f.readlines()
f.close()
for i in range(22):
    row=[int(n) for n in lines[i].strip().split()[1:]]
    for j in range(i+1):
        if i==j:
            matchcounts=addnum(matchcounts,row[j])
        else:
            offcounts=addnum(offcounts,row[j])

matchcounts=normalize(matchcounts)
offcounts=normalize(offcounts)

for i in range(22):
    row=[]
    for j in range(22):
        if i==j:
            row.append(drawnum(matchcounts))
            mat.append(row)
            break
        else:
            row.append(drawnum(offcounts))

for i in range(len(mat)):
    for j in range(i+1,len(mat)):
        mat[i].append(mat[j][i])

mat.append([int(i) for i in lines[22].strip().split()[1:]])
mat.append([int(i) for i in lines[23].strip().split()[1:]])

for i in range(22):
    mat[i].append(mat[22][i])
    mat[i].append(mat[23][i])

for i in mat:
    print(i)

f=open("randblos.mat","w")
for i in aastring:
    f.write("\t"+i)
f.write("\n")

for i in range(len(mat)):
    f.write(aastring[i])
    for j in mat[i]:
        f.write("\t"+str(j))
    f.write("\n")


