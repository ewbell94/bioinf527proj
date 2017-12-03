randblos<-read.csv("randblos.csv")
blosum62<-read.csv("blosum62.csv")
new99<-read.csv("new99.csv")
mats<-c(randblos,blosum62,new99)
labels<-c("randblos","blosum62","new99")
accmat<-as.data.frame(matrix(ncol=3,nrow=3))
precmat<-as.data.frame(matrix(ncol=3,nrow=3))
recmat<-as.data.frame(matrix(ncol=3,nrow=3))
meanmat<-as.data.frame(matrix(ncol=3,nrow=3))
rownames(meanmat)<-labels
rownames(accmat)<-labels
rownames(precmat)<-labels
rownames(recmat)<-labels
colnames(meanmat)<-colnames(blosum62)
colnames(accmat)<-labels
colnames(precmat)<-labels
colnames(recmat)<-labels

for(i in 1:3){
      meanmat[i,1]<-mean(mats[3*(i-1)+1]$accuracy)
      meanmat[i,2]<-mean(mats[3*(i-1)+2]$precision)
      meanmat[i,3]<-mean(mats[3*(i-1)+3]$recall)
}
meanmat

for(i in 2:3){
      for(j in 1:(i-1)){
      	    accmat[i,j]<-t.test(mats[3*(i-1)+1]$accuracy,mats[3*(j-1)+1]$accuracy)$p.value
	    precmat[i,j]<-t.test(mats[3*(i-1)+2]$precision,mats[3*(j-1)+2]$precision)$p.value
	    recmat[i,j]<-t.test(mats[3*(i-1)+3]$recall,mats[3*(j-1)+3]$recall)$p.value
      }
}

accmat
precmat
recmat