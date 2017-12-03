#Written by Eric Bell


#The Box object makes traceback easier by having a pointer to a previous box and how to get there
class Box(object):
    def __init__(self,val=1000,move="N",prev=None,gapcount=1):
        self.val=val #numeric value of the square
        self.move=move #move it took to get to this square
        self.prev=prev #pointer to the previous box
        self.gapcount=gapcount #when gaps happen, how big the jump is
    
    def tostring(self):
        print(self.val,self.move)
        
#printAlign does traceback and prints out the resulting alignment     
def printAlign(alignment,seq1,seq2,mutmat,gapstart,gapextend,argv):
    currbox=alignment[len(alignment)-1][len(alignment[0])-1]
    seq1i=len(seq1)-1
    seq2i=len(seq2)-1
    topstring=""
    midstring=""
    botstring=""
    aligncount=0
    idcount=0
    score=currbox.val
    predscore=0
    #traceback by repeatedly taking the "prev" square
    while currbox.prev:
        #currbox.tostring()
        #if match/mismatch move was taken
        if currbox.move == "M":
            aligncount+=1
            topstring+=seq1[seq1i]
            botstring+=seq2[seq2i]
            predscore+=mutmat[seq1[seq1i]+seq2[seq2i]]
            if seq1[seq1i]==seq2[seq2i]:
                midstring+=":"
                idcount+=1
            else:
                midstring+=" "
            seq1i-=1
            seq2i-=1
        #if horizontal gap was taken
        elif currbox.move == "H":
            for i in range(currbox.gapcount):
                topstring+=seq1[seq1i]
                if botstring[-1:]!="-":
                    predscore-=gapstart
                else:
                    predscore-=gapextend
                botstring+="-"
                midstring+=" "
                seq1i-=1
        #if vertical gap was taken
        elif currbox.move == "V":
            for i in range(currbox.gapcount):
                if topstring[-1:]!="-":
                    predscore-=gapstart
                else:
                    predscore-=gapextend
                topstring+="-"
                botstring+=seq2[seq2i]
                midstring+=" "
                seq2i-=1
        #this should never happen
        else:
            print("Something went horribly wrong")
        currbox=currbox.prev
    #print("Length of Sequence 1: "+str(len(seq1))+" -> "+argv[1])
    #print("Length of Sequence 2: "+str(len(seq2))+" -> "+argv[2])
    #print("Alignment score: "+str(score)+" ("+str(predscore)+" predicted)")
    #print("Aligned Length: "+str(aligncount))
    #print("Identical Length: "+str(idcount))
    #print("Sequence Identity: "+str(float(idcount)/float(len(seq2)))+" ("+str(idcount)+"/"+str(len(seq2))+")")
    #print("")
    #these strings need to be reversed because they were constructed 3'->5'
    #print(topstring[::-1])
    #print(midstring[::-1])
    #print(botstring[::-1])
    #numstring=""
    #for i in range(1,len(botstring)+1):
    #    numstring+=str(i%10)
    #print(numstring)
    return float(idcount)/float(len(seq2))

#makeAlign performs sequence alignment between seq1 and seq2 using the Gotoh/NW algorithm
def makeAlign(seq1,seq2,gapstart,gapextend,mutmat):
    alignment=[]
    H=[]
    V=[]
    #set up S, H, and V and initialize
    for i in range(len(seq2)+1):
        row=[]
        hrow=[]
        vrow=[]
        for j in range(len(seq1)+1):
            row.append(Box())
            hrow.append(0)
            vrow.append(0)
        alignment.append(row)
        H.append(hrow)
        V.append(vrow)
    Vjump=[0]*len(alignment[0])
    for i in range(len(alignment)):
        Hjump=0
        for j in range(len(alignment[0])):
            if i == 0: #initialize top row
                if j == 0:
                    alignment[i][j].val=0
                else:
                    alignment[i][j].prev=alignment[i][j-1]
                    alignment[i][j].move="H"
                    if alignment[i][j].prev.move=="H":
                        alignment[i][j].val=alignment[i][j].prev.val-gapextend
                    else:
                        alignment[i][j].val=alignment[i][j].prev.val-gapstart
                    H[i][j]=alignment[i][j].val #No matches have been made so S_0j=H_0j
                    V[i][j]=float("-inf") #V[0][j] is not defined, this value ensures it doesn't get used 
            elif j == 0: #initialize first column
                alignment[i][j].prev=alignment[i-1][j]
                alignment[i][j].move="V"
                if alignment[i][j].prev.move=="V":
                    alignment[i][j].val=alignment[i][j].prev.val-gapextend
                else:
                    alignment[i][j].val=alignment[i][j].prev.val-gapstart
                H[i][j]=float("-inf") #H[i][0] is not defined, this value ensures it doesn't get used
                V[i][j]=alignment[i][j].val #No matches have been made so S_i0=V_i0
            else: #for the rest of the cells
                H[i][j]=max([alignment[i][j-1].val-gapstart,H[i][j-1]-gapextend])
                if H[i][j]==alignment[i][j-1].val-gapstart:
                    Hjump=j-1
                V[i][j]=max([alignment[i-1][j].val-gapstart,V[i-1][j]-gapextend])
                if V[i][j]==alignment[i-1][j].val-gapstart:
                    Vjump[j]=i-1
                hval=H[i][j]
                vval=V[i][j]
                mval=alignment[i-1][j-1].val+mutmat[seq1[j-1]+seq2[i-1]]
                if mval > hval and mval > vval: #if blosum match/mismatch is the best move
                    alignment[i][j].prev=alignment[i-1][j-1]
                    alignment[i][j].move="M"
                    alignment[i][j].val=mval
                elif vval > hval: #if vertical is the better gap, it will take it
                    alignment[i][j].val=vval
                    alignment[i][j].prev=alignment[Vjump[j]][j]
                    alignment[i][j].move="V"
                    alignment[i][j].gapcount=i-Vjump[j]
                else: #if all values are equal, a horizontal gap will be taken
                    alignment[i][j].val=hval
                    alignment[i][j].prev=alignment[i][Hjump]
                    alignment[i][j].move="H"
                    alignment[i][j].gapcount=j-Hjump
    #this prints out matrices for debugging                
    '''S=[]
    for i in range(len(alignment)):
        row=[]
        for j in range(len(alignment[0])):
            row.append(alignment[i][j].val)
        S.append(row)
    for i in S:
        print(i)
    print("")
    for i in H:
        print(i)
    print("")
    for i in V:
        print(i)'''
    return alignment

#parseMutMat parses the mutation matrix file into a dictionary
#Note that redundancy is preserved for ease of implementation despite using more memory
def parseMutMat(matname):
    f=open(matname)
    mutmat={}
    toprow=f.readline().strip().split()
    for line in f:
        parts=line.strip().split()
        for i in range(1,len(parts)):
            mutmat[toprow[i-1]+parts[0]]=int(parts[i])
    f.close()
    return mutmat

#This is the function that gets called when the program runs            
def main(argv):
    if len(argv) < 3 or len(argv) > 5:
        print("\nIncorrect number of arguments\n")
        print("Program Arguments:\n"+
              "1. Sequence 1 (either string or FASTA filename)\n"+
              "2. Sequence 2 (either string or FASTA filename)\n"+
              "3. Filename of mutation matrix (optional, defaults to blosum62.mat)\n"+
              "4. Gap opening and extension penalty (optional, defaults to 11,1)\n"
        )
        exit(1)
    if "." in argv[1]: #if filename
        f=open(argv[1])
        f.readline()
        seq1=""
        for line in f:
            seq1+=line.strip()
    else: #otherwise interpret as sequence
        seq1=argv[1]
    if "." in argv[2]: #if filename
        f=open(argv[2])
        f.readline()
        seq2=""
        for line in f:
            seq2+=line.strip()
    else: #otherwise interpret as sequence
        seq2=argv[2]
    if len(argv) >= 4: #specify mutation matrix
        mutmat=parseMutMat(argv[3])
        if len(argv) == 5: #specify gap penalties (positive number)
            parts=argv[4].split(",")
            gapstart=int(parts[0])
            gapextend=int(parts[1])
        else:
            gapstart=11
            gapextend=1
    else: #otherwise use defaults
        mutmat=parseMutMat("blosum62.mat")
        gapstart=11
        gapextend=1
    alignment=makeAlign(seq1,seq2,gapstart,gapextend,mutmat) #do the alignment
    return printAlign(alignment,seq1,seq2,mutmat,gapstart,gapextend,argv) #traceback and print out the alignment

#from sys import argv
#main(argv)
