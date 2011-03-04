require(biclust)
# example data set according to the simulation study in Lee et al. 2010
# generate artifical data set and a correspondig biclust object
u <- c(10,9,8,7,6,5,4,3,rep(2,17),rep(0,75))
v <- c(10,-10,8,-8,5,-5,rep(3,5),rep(-3,5),rep(0,34))
u <- u/sqrt(sum(u^2)) 
v <- v/sqrt(sum(v^2))
d <- 50
X <- (d*u%*%t(v)) + matrix(rnorm(100*50),100,50)
params <- info <- list()
RowxNumber <- matrix(rep(FALSE,100),ncol=1)
NumberxCol <- matrix(rep(FALSE,50),nrow=1)
RowxNumber[u!=0,1] <- TRUE 
NumberxCol[1,v!=0] <- TRUE
Number <- 1
ressim <- BiclustResult(params,RowxNumber,NumberxCol,Number,info)

#perform s4vd biclustering 
system.time( ress4vd <- biclust(X,method=BCs4vd,pcerv=0.5,pceru=0.5,ss.thr=c(0.6,0.65),steps=500,pointwise=F,nbiclust=1,savepath=TRUE))
#perform s4vd biclustering with fast pointwise stability selection
system.time( ress4vdpw <- biclust(X,method=BCs4vd,pcerv=0.5,pceru=0.5,ss.thr=c(0.6,0.65),steps=500,pointwise=T,nbiclust=1))
#perform ssvd biclustering
system.time(resssvd <- biclust(X,BCssvd,K=1))
#agreement of the results with the simulated bicluster
jaccardind(ressim,ress4vd)
jaccardind(ressim,ress4vdpw)
jaccardind(ressim,resssvd)

#heatmap plot
BCheatmap(X,ress4vd)
#stability paths
stabpath(ress4vd,1)
#selection probabilitys for the pointwise stability selection
stabpath(ress4vdpw,1)
#parallel coordinates
parallelCoordinates(X,ress4vd,1,plotBoth=TRUE, compare=TRUE,order = TRUE)

#lung cancer data set   Bhattacharjee et al. 2001
data(lung)
set.seed(12)
res2 <- s4vd(lung,pcerv=.5,pceru=0.01,ss.thr=c(0.6,0.65)
,start.iter=3,size=0.632,cols.nc=T,steps=100,pointwise=T
,merr=0.0001,iter=100,nbiclust=10,col.overlap=F)
BCheatmap(lung,res2)