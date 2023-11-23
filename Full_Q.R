fq <- function(X){
  r<-dim(X)[1]
  c<-dim(X)[2]
  X.n<-matrix(rep(0,r*c),r,c)
  Y<-X.n
  for(i in 1:c){
    Y[,i]<-sort(X[,i])
  }
  head(Y)
  m<-apply(Y,1,mean)
  X.1<-as.vector(0)
  for(j in 1:c){
    X.1[order(X[,j])]<-m
    X.n[,j]<-as.matrix(X.1,r,1)
  }
  head(X.n)
  #         [,1]     [,2]      [,3]
  #[1,] 3.614902 1.595378 0.8513045
  colnames(X.n)<-colnames(X)
  rownames(X.n)<-rownames(X)
  return(X.n)
}