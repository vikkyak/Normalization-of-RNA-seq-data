med <- function(X){
  
  MED<-function(y){
    median(y[y>0])
  }
  X<-X+0.1
  med<-apply(X,2, MED)
  f.med<-med/mean(med)
  med.res<-scale(X,center=FALSE, scale=f.med)
  return(med.res)
}
