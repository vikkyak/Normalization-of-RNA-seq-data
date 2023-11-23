uq <- function(X){
  
  #excluding zero counts in each sample
  UQ<-function(y){
    quantile(y, 0.75)
  }
  X<-X+0.1
  upperQ<-apply(X,2,UQ)
  f.uq<-upperQ/mean(upperQ)
  upq.res<-scale(X,center=FALSE,scale=f.uq) 
  return(upq.res)
}
