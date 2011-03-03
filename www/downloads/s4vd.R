#The function performs biclustering of the data matrix by sparse singular value decomposition with nested stability selection.
#arguments:
#
#	x 				The matrix to be clustered.
#   steps			Number of subsamples used to perform the stability selection
#   pcerv			Per comparsion wise error rate to control the number of falsely selected right singular vector coefficients (columns/samples).
#	pceru			Per comparsion wise error rate to control the number of falsely selected left singular vector coefficients (rows/genes).
#	ss.thr			Range of the cutoff threshold (relative selection frequency) for the stability selection.
#	size			Size of the subsamples used to perform the stability selection.  
#	gamm			Weight parameter for the adaptive LASSO, nonnegative constant (default = 0, LASSO).
#	iter			Maximal number of iterations to fit a single bicluster.
#	nbiclust		Maximal number of biclusters. 
#	merr			Threshold to decide convergence. 
#	cols.nc			Allow for negative correlation of columns (samples) over rows (genes).
#	rows.nc			Allow for negative correlation of rows (genes) over columns (samples).
#	row.overlap		Allow rows to overlap between biclusters. 
#	col.overlap		Allow columns to overlap between biclusters. 
#	row.min			Minimal number of rows.
#	col.min			Minimal number of columns.
#	pointwise		If TRUE performs a fast pointwise stability selection instead of calculating the complete stability path.  
#	start.iter		Number of starting iterations in which the algorithm is not allowed to converge. 
#	savepath		Saves the stability path in order plot the path with the stabpathplot function.
	
#Note that pointwise needs to be TRUE to save the path. For extreme high dimensional data sets (e.g. the lung cancer example) the resulting
#biclust object may exceed the available memory.



s4vd <- function(
		X,
		steps = 100,
		pcerv = 0.1,
		pceru = 0.1,
		ss.thr = c(0.6,0.65),
		size = 0.5,
		gamm = 0,
		iter = 20,
		nbiclust = 10,
		merr = 10^(-4),
		cols.nc=TRUE,
		rows.nc=TRUE,
		row.overlap=TRUE,
		col.overlap=TRUE,
		row.min=1,
		col.min=1,
		pointwise=TRUE,
		start.iter=3,
		savepath=FALSE
){
	MYCALL<-match.call()
	startX <- X
	p.ini <- nrow(X)
	n.ini <- ncol(X)
	rowsin <- rep(TRUE,p.ini)	
	colsin <- rep(TRUE,n.ini)
	stop <- FALSE
	start <- TRUE
	info <- Rows <- Cols <- vc <- uc <- list()
	for(k in 1:nbiclust){
		gc()
		cat("Bicluster",k)
		rows <- rep(FALSE,nrow(startX))
		cols <- rep(FALSE,ncol(startX))
		if(is.null(nrow(X))|is.null(ncol(X))){
			number <- k-1
			stop <- TRUE
			break
		}
		if(nrow(X)==0|ncol(X)==0){
			number <- k-1
			stop <- TRUE
			break
		}
		SVD <- svd(X,nu=1,nv=1)
		v0 <- SVD$v
		u0 <- SVD$u
		d0 <- SVD$d
		vc <- uc <- list()
		if((length(u0)*size)<=2|(length(v0)*size)<=2){
			cat("submatrix to small for resampling","\n")
			number <- k-1
			stop <- TRUE
			break
		}
		if(pointwise){
			for(i in 1:iter){
				if(i > start.iter) start <- FALSE
				uc <- updateu.pw(X,v0,pceru,p.ini,ss.thr,steps,size,gamm,rows.nc,uc$l)
				u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2)) 
				u1[is.na(u1)] <- 0
				vc <- updatev.pw(X,u1,pcerv,n.ini,ss.thr,steps,size,gamm,cols.nc,vc$l)
				v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2)) 
				v1[is.na(v1)] <- 0
				if(uc[[3]] & i > start.iter){
					cat("rows not stable")
					stop <- TRUE
					break}
				if(vc[[3]] & i > start.iter){
					cat("columns not stable")
					stop <- TRUE
					break}
				ud <- sqrt(sum((u0-u1)^2))
				vd <- sqrt(sum((v0-v1)^2))
				#cat("iter: ",i," rows: ",sum(uc$sp>=uc$thr)," cols: ",sum(vc$sp>=uc$thr)
				#		," merr: ",min(c(ud,vd)),"row.thr:",uc[[5]],"col.thr",vc[[5]],"\n")
				cat(".")
				v0 <- v1
				u0 <- u1
				if(min(c(vd,ud)) < merr & i > start.iter)break
			}
		}else{
			for(i in 1:iter){
				if(i > start.iter) start <- FALSE
				uc <- updateu(X,v0,pceru,p.ini,ss.thr,steps,size,gamm,rows.nc,savepath)
				u1 <- uc[[1]]/sqrt(sum(uc[[1]]^2)) 
				u1[is.na(u1)] <- 0
				vc <- updatev(X,u1,pcerv,n.ini,ss.thr,steps,size,gamm,cols.nc,savepath)
				v1 <- vc[[1]]/sqrt(sum(vc[[1]]^2)) 
				v1[is.na(v1)] <- 0
				if(uc[[3]] & i > start.iter){
					cat("rows not stable")
					stop <- TRUE
					break}
				if(vc[[3]] & i > start.iter){
					cat("columns not stable")
					stop <- TRUE
					break}
				ud <- sqrt(sum((u0-u1)^2))
				vd <- sqrt(sum((v0-v1)^2))
				#cat("iter: ",i," rows: ",sum(uc$sp>=uc$thr)," cols: ",sum(vc$sp>=uc$thr)
				#		," merr: ",min(c(ud,vd)),"row.thr:",uc[[5]],"col.thr",vc[[5]],"\n")
				cat(".")
				v0 <- v1
				u0 <- u1
				if(min(c(vd,ud)) < merr & i > start.iter)break
			}
		}
		stableu <- uc$sp >= uc$thr
		stablev <- vc$sp >= vc$thr
		d0 <- as.numeric(t(u0)%*%X%*%v0)
		u0[!stableu] <- 0
		v0[!stablev] <- 0
		rows[rowsin] <- u0!=0
		cols[colsin] <- v0!=0
		Rows[[k]] <- rows
		Cols[[k]] <- cols
		if(stop){
			number <- k-1
			break
		}
		if(i==iter){
			number <- k-1
			stop <- TRUE
			cat("Fail to converge! Increase the number of iterations !","\n")
			gc()
			break
		}
		if(!row.overlap){
			rowsin[rows] <- FALSE
			X <- startX[rowsin,colsin]
			info[[k]] <- list(vc,uc,layer=list(u0,v0,d0))
		} 
		if(!col.overlap){
			colsin[cols] <- FALSE
			X <- startX[rowsin,colsin]
			info[[k]] <- list(vc,uc,layer=list(u0,v0,d0))
		} 
		if(row.overlap&col.overlap){
			temp <- svd(X[rows,cols]) 
			#X <- X - (d0*u0%*%t(v0))
			X[rows,cols] <- X[rows,cols] - (temp$d[1]*temp$u[,1]%*%t(temp$v[,1]))
			info[[k]] <- list(vc,uc,layer=list(u0,v0,d0))
		}
		cat("\n")
	}
	if(!stop) number <- k
	params <- list(steps = steps, pcerv=pcerv, pceru=pceru, iter=iter, ss.thr=ss.thr, size=size, gamm=gamm, row.overlap=row.overlap, col.overlap=col.overlap,
			rows.nc=rows.nc, cols.nc=cols.nc, nbiclust=nbiclust, merr=merr, row.min=row.min, col.min=col.min, pointwise=pointwise, start.iter=start.iter, savepath=savepath, Call=MYCALL)  
	RowxNumber=t(matrix(unlist(Rows),byrow=T,ncol=length(Rows[[1]])))
	NumberxCol=matrix(unlist(Cols),byrow=T,ncol=length(Cols[[1]]))
	RowxNumber <- matrix(RowxNumber[,1:number],ncol=number)
	NumberxCol <- matrix(NumberxCol[1:number,],nrow=number)
	Number <- number
	return(BiclustResult(params,RowxNumber,NumberxCol,Number,info))
}


#update v
updatev <- function(X,u0,pcer,n.ini,ss.thr,steps,size,gamm,cols.nc=FALSE,savepath=FALSE){
	n.ini  <-	n <- ncol(X)
	err <- pcer*n.ini
	ols <- t(X)%*%u0	
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	if(savepath) selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(cols.nc){
		for(l in 1:length(lambdas)){
			temp <- adaLasso.nc(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}else{
		for(l in 1:length(lambdas)){
			temp <- adaLasso(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}
	thr <- thrall[ls]
	if(thr > ss.thr[2]){
		while(pcer <= 0.5){
			pcer <- pcer + 0.01	
			thrall <- ((qs^2/((pcer*n.ini)*n.ini))+1)/2
			thr <- thrall[ls]
			if(thr < ss.thr[2])break
		}
	} 
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambdas[ls]/(abs(ols)^gamm)  
	vc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	if(savepath){
		return(list(vc=vc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,selprobpath=selprobpath))
	}else{
		return(list(vc=vc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls))
	}
}

#update u
updateu <- function(X,v0,pcer,p.ini,ss.thr,steps,size,gamm,rows.nc=FALSE,savepath=FALSE,start=FALSE){
	p.ini <- p <- nrow(X)
	err <- pcer*p.ini
	ols <- X%*%v0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	if(savepath) selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(rows.nc){
		for(l in 1:length(lambdas)){
			temp <- adaLasso.nc(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}else{
		for(l in 1:length(lambdas)){
			temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			if(savepath)selprobpath[,l] <- sp
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l]>=ss.thr[1]){
				ls <- l
				break
			}
		}
	}
	thr <- thrall[ls]
	if(thr > ss.thr[2]){
		while(pcer <= 0.5){
			pcer <- pcer + 0.01	
			thrall <- ((qs^2/((pcer*p.ini)*p.ini))+1)/2
			thr <- thrall[ls]
			if(thr < ss.thr[2])break
		}
	} 
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambdas[ls]/(abs(ols)^gamm)  
	uc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	if(savepath){
		return(list(uc=uc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,selprobpath=selprobpath))
	}
	else{
		return(list(uc=uc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls))
	}	
}

#update v pointwise

updatev.pw <- function(X,u0,pcer,n.ini,ss.thr,steps,size,gamm,cols.nc=FALSE,l=NULL,start=FALSE){
	n.ini <- n <- ncol(X)
	err <- pcer*n.ini
	ols <- t(X)%*%u0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	#search for a lambda
	l.min <- 1
	l.max <- length(lambdas)
	if(cols.nc){
		for(g in 1:(length(lambdas))){
			temp <- adaLasso.nc(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){          
				l.min <- l
				if(l == length(lambdas))break   
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					if(thrall[l+1]>1) ls <-l
					temp <- adaLasso.nc(t(X),u0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*n.ini))+1)/2
					break
				}
				l <- min(length(lambdas),l.max,l + ceiling(length(lambdas)/(g+1))) 
				while(thrall[l]!=0){  
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				l.max <- l
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l.min,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0 ){ 
					l <- l+1
					if(l == length(l))break
				} 
			}
		}
	}else{
		for(g in 1:(length(lambdas))){
			temp <- adaLasso(t(X),u0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*n.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){          
				l.min <- l
				if(l == length(lambdas))break   
				if(thrall[l+1]> ss.thr[2]){
					ls <- l +1 
					if(thrall[l+1]>1) ls <-l
					temp <- adaLasso(t(X),u0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*n.ini))+1)/2
					break
				}
				l <- min(length(lambdas),l.max,l + ceiling(length(lambdas)/(g+1))) 
				while(thrall[l]!=0 ){  
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				l.max <- l
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l.min,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0 ){ 
					l <- l+1
					if(l == length(lambdas))break
				} 
			}
		}
	}	
	#thr <- thrall[ls]
	#if(thr > ss.thr[2]){
	#	while(pcer <= 0.5){
	#		pcer <- pcer + 0.01	
	#		thrall <- ((qs^2/((pcer*n.ini)*n.ini))+1)/2
	#		thr <- thrall[ls]
	#		if(thr < ss.thr[2])	break
	#	}
	#}
	thr <- ((qs[ls]^2/((pcer*n.ini)*n.ini))+1)/2
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	vc <- rep(0,n)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	vc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	return(list(vc=vc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,delta=delta))
}


#update u pointwise
updateu.pw <- function(X,v0,pcer,p.ini,ss.thr,steps,size,gamm,rows.nc=FALSE,l=NULL,start=FALSE){
	p.ini <- p <- nrow(X)
	err <- pcer*p.ini
	ols <- X%*%v0
	stop <- FALSE
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	qs <- numeric(length(lambdas)) 
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	#search for a lambda
	l.min <- 1
	l.max <- length(lambdas)
	if(rows.nc){
		for(g in 1:(length(lambdas))){
			temp <- adaLasso.nc(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				l.min <- l
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					if(thrall[l+1]>1) ls <-l
					temp <- adaLasso.nc(X,v0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*p.ini))+1)/2
					break
				}
				l <- min(length(lambdas),l.max,l + ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0 ){  # if thr for current lambda available decrease lambda 
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				l.max <- l
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l.min,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0){ 
					l <- l+1
					if(l == length(l))break
				} 
			}
		}
	}else{
		for(g in 1:(length(lambdas))){
			temp <- adaLasso(X,v0,lambdas[l],steps,size,gamm)
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			if(thrall[l] >= ss.thr[1] & thrall[l] <= ss.thr[2]){
				ls <- l
				break
			} 
			if(thrall[l] < ss.thr[1]){
				l.min <- l
				if(l == length(lambdas))break
				if(thrall[l+1]> ss.thr[2]){
					ls <- l+1
					if(thrall[l+1]>1) ls <-l
					temp <- adaLasso(X,v0,lambdas[ls],steps,size,gamm)
					t <- temp!=0
					qs[ls] <- mean(colSums(t))
					sp <- rowMeans(t)
					thrall[ls] <- ((qs[ls]^2/(err*p.ini))+1)/2
					break
				}
				l <- min(length(lambdas),l.max,l + ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0 ){  
					l <- l-1
					if(l == 0)break
				} 
			}
			if(thrall[l] > ss.thr[2]){ 
				l.max <- l
				if(l == 0)break
				if(thrall[l-1]<ss.thr[1]&thrall[l-1]!=0){
					ls <- l
					break
				}
				l <- max(1,l.min,l - ceiling(length(lambdas)/(g+1)))
				while(thrall[l]!=0){ 
					l <- l+1
					if(l == length(lambdas))break
				} 
			}
		}
	}
	#thr <- thrall[l]
	#if(thr > ss.thr[2]){
	#	while(pcer <= 0.5){
	#		pcer <- pcer + 0.01	
	#		thrall <- ((qs^2/((pcer*p.ini)*p.ini))+1)/2
	#		thr <- thrall[ls]
	#		if(thr < ss.thr[2])break
	#	}
	#}
	thr <- ((qs[ls]^2/((pcer*p.ini)*p.ini))+1)/2
	stable <- which(sp>=thr)
	if(length(stable)==0)stop <- TRUE
	uc <- rep(0,p)
	delta <- lambdas[l]/(abs(ols)^gamm)  
	uc <- (sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta))
	return(list(uc=uc,sp=sp,stop=stop,qs=qs,thr=thr,l=ls,delta=delta))
}

#adaptive Lasso 
adaLasso.nc <- function(X,b,lambda,steps,size,gamm=0){
	subsets <- sapply(1:steps,function(x){sample(1:length(b),length(b)*size)})
	res <- sapply(1:steps,adaLassoSteps.nc,subsets,X,b,lambda,gamm)
	return(res)
}
adaLassoSteps.nc <- function(index,subsets,X,b,lambda,gamm){
	ols <- X[,subsets[,index]]%*%b[subsets[,index]]
	delta <- lambda/(abs(ols)^gamm)                        
	ols <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	ols[is.na(ols)] <- 0
	return(ols)
}
adaLasso <- function(X,b,lambda,steps,size,gamm=0){
	subsets <- sapply(1:steps,function(x){sample(1:length(b),length(b)*size)})
	res <- sapply(1:steps,adaLassoSteps,subsets,X,b,lambda,gamm)
	return(res)
}
adaLassoSteps <- function(index,subsets,X,b,lambda,gamm){
	ols <- X[,subsets[,index]]%*%b[subsets[,index]]
	delta <- lambda/(abs(ols)^gamm)                        
	ols <- sign(ols)*(abs(ols)>=delta)*(abs(ols)-delta)
	ols[is.na(ols)] <- 0
	mostof <- sign(sum(sign(ols)))
	if(mostof==0) mostof <- 1
	ols[which(sign(ols) != mostof) ] <- 0
	ols[is.na(ols)] <- 0
	return(ols)
}

########################ssvd algorithm

ssvdBC <- function(X,K=10,threu = 1, threv = 1, gamu = 0, gamv =0 , merr = 10^(-4), niter = 100){
		MYCALL <- match.call()
		res <- list()
		RowxNumber <- matrix(nrow=nrow(X),ncol=K)
		NumberxCol <- matrix(ncol=ncol(X),nrow=K)
		for(k in 1:K){ 
			res[[k]]	<- ssvd(X,threu = 1, threv = 1, gamu = 0, gamv =0,  u0 = svd(X)$u[,k], v0 = svd(X)$v[,k], merr = 10^(-4), niter = 100)
			if(res[[k]]$stop){
				K <- k-1
				break
			}
			RowxNumber[,k] <- res[[k]][[1]]!=0
			NumberxCol[k,] <- res[[k]][[2]]!=0
			#as.numeric(t(u0)%*%X%*%v0)
			d <- as.numeric(t(res[[k]][[1]])%*%X%*%res[[k]][[2]])
			res[[k]][[4]] <- d
			X <- X - (d*res[[k]][[1]]%*%t(res[[k]][[2]]))
		}
		params <- list(K=K,threu=threu,threv=threv,gamu=gamv,merr=merr,niter=niter,Call=MYCALL)
		Number <- K
		info <- list(res=res)
		return(BiclustResult(params,RowxNumber,NumberxCol,Number,info))
}

## Implementation of the sparse SVD
#
#	work with thresh.R 
#
#  	Input variables:
#       X - argument (n x d matrix)
#       threu = type of penalty (thresholding rule) for the left 
#       	    singular vector,
#				  1 = (Adaptive) LASSO (default)
#                 2 = hard thresholding
#                
#       threv = type of penalty (thresholding rule) for the right 
#               singular vector,
#                 1 = (Adaptive) LASSO (default)
#                 2 = hard thresholding
#
#       gamu = weight parameter in Adaptive LASSO for the 
#              left singular vector, nonnegative constant (default = 0, LASSO)
#
#       gamv = weight parameter in Adaptive LASSO for the  
#              right singular vector, nonnegative constant (default = 0, LASSO)
#
#       u0,  v0 = initial values of left/right singular vectors                
#                 (default = the classical SVs)
# 
#       merr = threshold to decide convergence (default = 10^(-4))
#
#       niter = maximum number of iterations (default = 100)
#
#   Output: 
#       u = left sparse singular vector
#       v = right sparse singaulr vector
#       iter = number of iterations to achieve the convergence


ssvd = function(X,threu = 1, threv = 1, gamu = 0, gamv =0,  u0 = svd(X)$u[,1], v0 = svd(X)$v[,1], merr = 10^(-4), niter = 100){
    n = dim(X)[1]
    d = dim(X)[2]
	stop <- FALSE
    ud = 1;
    vd = 1;
    iter = 0;
    SST = sum(X^2);

    while (ud > merr || vd > merr) {
        iter = iter+1;
		cat("iter: ",iter,"\n")
		cat("v: ",length(which(v0!=0)),"\n")
		cat("u: ",length(which(u0!=0)),"\n")
        # Updating v
        z =  t(X)%*% u0;
        winv = abs(z)^gamv;
        sigsq = abs(SST - sum(z^2))/(n*d-d);

        tv = sort(c(0, abs(z*winv)));
        rv = sum(tv>0);
        Bv = rep(1,d+1)*Inf;
    
        for (i in 1:rv){
            lvc  =  tv[d+1-i];
            temp1 = which(winv!=0);
            
            temp2 = thresh(z[temp1], type = threv, delta = lvc/winv[temp1]);
            vc = rep(0,d)
            vc[temp1] = temp2;
            Bv[i] = sum((X - u0%*%t(vc))^2)/sigsq + i*log(n*d)
        }
    
        Iv = min(which(Bv==min(Bv)))
        temp = sort(c(0, abs(z*winv)));        lv = temp[d+1-Iv]
    
        temp2 = thresh(z[temp1],type = threv, delta = lv/winv[temp1]);
        v1 = rep(0,d)
        v1[temp1] = temp2;
        v1 = v1/sqrt(sum(v1^2)) #v_new
	
		cat("v1", length(which(v1!=0)) ,"\n" )
		#str(X)
		#str(v1)
        # Updating u
        z = X%*%v1;
        winu = abs(z)^gamu;
        sigsq = abs(SST - sum(z^2))/(n*d-n)

        tu = sort(c(0,abs(z*winu)));
        ru = sum(tu>0);
        Bu = rep(1,n+1)*Inf;

        for (i in 1:ru){
            luc  =  tu[n+1-i];
            temp1 = which(winu!=0);
            temp2 = thresh(z[temp1], type = threu, delta = luc/winu[temp1]);
            uc = rep(0,n)
            uc[temp1] = temp2;
            Bu[i] = sum((X - uc%*%t(v1))^2)/sigsq + i*log(n*d)
        }
        
        Iu = min(which(Bu==min(Bu)));
        temp = sort(c(0, abs(z*winu)))
        lu = temp[n+1-Iu];
	  temp2 = thresh(z[temp1],delta = lu/winu[temp1]);
        u1 = rep(0,n);
        u1[temp1] = temp2;
        u1 = u1/sqrt(sum(u1^2));

    
        ud = sqrt(sum((u0-u1)^2));
        vd = sqrt(sum((v0-v1)^2));

        if (iter > niter){
        print("Fail to converge! Increase the niter!")
		stop <- TRUE
	  break
        }
        
	u0 = u1;
	v0 = v1;
    }

return(list(u = u1, v = v1, iter = iter,stop=stop))
}


# heatmap
BCheatmap <- function(
		X,res,
		cexR=.75,
		cexC=.75,
		axisR=FALSE,
		axisC=TRUE,
		heatcols = diverge_hcl(25, h = c(260, 0), c = 80, l = c(30, 100), power = 1.5,gamma = 2.4, fixup = TRUE),
		clustercols= rainbow_hcl(res@Number, c = 100, l = 50),
		allrows=F,
		allcolumns=T
)
{
	number <- res@Number
	layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE),widths=c(8,2),heights=c(1,1))
	if(number==1){
		rowmat <- res@RowxNumber
		colmat <- t(res@NumberxCol)
		roworder <- c(which(res@RowxNumber[,1]),which(!res@RowxNumber[,1]))
		colorder <- c(which(res@NumberxCol[1,]),which(!res@NumberxCol[1,]))
		X <- X[roworder,colorder]
		par(omi=c(.2,0,0,.5)) 
		image(t(X)[,nrow(X):1],col=heatcols,x=c(1:ncol(X)),y=c(1:nrow(X)),axes=F,ylab="",xlab="")
		if(axisC)axis(1,1:ncol(X),labels=colnames(X),las=2,line = -0.5, tick = 0,cex.axis = cexC)
		if(axisR)axis(4,1:nrow(X),labels=rownames(X),las=2,line = -0.5, tick = 0,cex.axis = cexR)
		rin1 <- which(roworder %in% which(rowmat[,1]))
		cin1 <- which(colorder %in% which(colmat[,1]))
		nr <- length(roworder)
		nc <- length(colorder)
		if(allrows) nr <- nrow(X)
		if(allcolumns) nc <- ncol(X)
		xl <- 0.5
		yb <- nr-length(rin1)+0.5
		xr <- length(cin1)+0.5
		yt <- nr+0.5
		rect(xleft=xl,ybottom=yb,xright=xr,ytop=yt,density=0,angle=25,lwd=2.5,col=clustercols[1])
	}else{
	rowmat <- res@RowxNumber
	overlap <- rowSums(rowmat)
	roworder <- which(overlap==number) # in all
	if(number>2){
		for(i in 1:(number-2)){
			innext <- intersect(which(rowmat[,i]),which(rowmat[,i+1]))
			nooverlap <- which(rowmat[,i]&rowSums(rowmat)==1)
			for(l in 1:(number-1-i)){ 
				temp   <- intersect(which(rowmat[,i]),which(rowmat[,i+1+l])) 
				temp   <- temp[!intersect(which(rowmat[,i]),which(rowmat[,i+1+l])) %in% innext]
				roworder <- unique(c(roworder,temp))
			}
			roworder <- unique(c(roworder,nooverlap))
			roworder <- unique(c(roworder,innext))
		}
	}
	innext <- intersect(which(rowmat[,number-1]),which(rowmat[,number]))
	nooverlap <- which(rowmat[,number-1]&rowSums(rowmat)==1)
	roworder <- unique(c(roworder,nooverlap))
	roworder <- unique(c(roworder,innext))
	nooverlap <- which(rowmat[,number]&rowSums(rowmat)==1)
	roworder <- unique(c(roworder,nooverlap))
	if(allrows) roworder <- c(roworder,which(!1:nrow(rowmat)%in%roworder)) 
	colmat <- t(res@NumberxCol)
	overlap <- rowSums(colmat)
	colorder <- which(overlap==number) # in all
	if(number>2){
		for(i in 1:(number-2)){
			innext <- intersect(which(colmat[,i]),which(colmat[,i+1]))
			nooverlap <- which(colmat[,i]&rowSums(colmat)==1)
			for(l in 1:(number-1-i)){ 
				temp   <- intersect(which(colmat[,i]),which(colmat[,i+1+l])) 
				temp   <- temp[!intersect(which(colmat[,i]),which(colmat[,i+1+l])) %in% innext]
				colorder <- unique(c(colorder,temp))
			}
			colorder <- unique(c(colorder,nooverlap))
			colorder <- unique(c(colorder,innext))
		}
	}	
	innext <- intersect(which(colmat[,number-1]),which(colmat[,number]))
	nooverlap <- which(colmat[,number-1]&rowSums(colmat)==1)
	colorder <- unique(c(colorder,nooverlap))
	colorder <- unique(c(colorder,innext))
	nooverlap <- which(colmat[,number]&rowSums(colmat)==1)
	colorder <- unique(c(colorder,nooverlap))
	if(allcolumns) colorder <- c(colorder,which(!1:nrow(colmat)%in%colorder)) 
	X <- X[roworder,colorder]
	par(mar=c(4, 1, 1, 3))# c(bottom, left, top, right)
	image(t(X)[,nrow(X):1],col=heatcols,x=c(1:ncol(X)),y=c(1:nrow(X)),axes=F,ylab="",xlab="")
	if(axisC)axis(1,1:ncol(X),labels=colnames(X),las=2,line = -0.5, tick = 0,cex.axis = cexC)
	if(axisR)axis(4,1:nrow(X),labels=rownames(X),las=2,line = -0.5, tick = 0,cex.axis = cexR)
	rin1 <- which(roworder %in% which(rowmat[,1]))
	cin1 <- which(colorder %in% which(colmat[,1]))
	nr <- length(roworder)
	nc <- length(colorder)
	if(allrows) nr <- nrow(X)
	if(allcolumns) nc <- ncol(X)
	xl <- 0.5
	yb <- nr-length(rin1)+0.5
	xr <- length(cin1)+0.5
	yt <- nr+0.5
	rect(xleft=xl,ybottom=yb,xright=xr,ytop=yt,density=0,angle=25,lwd=2.5,col=clustercols[1])
	for(i in 2:number){
			rin <- which(roworder %in% which(rowmat[,i]))
			rstart <- numeric()
			rstop <- numeric()
			e <- 1
			rstart[e] <- rin[1]
			for(j in 2:length(rin)){
				if(rin[j-1]-rin[j]!=-1){
					rstop[e] <- rin[j-1]
					e <- e+1
					rstart[e] <- rin[j]
				}
			}
			rstop[e] <- rin[j]
			
			cin <- which(colorder %in% which(colmat[,i]))
			cstart <- numeric()
			cstop <- numeric()
			e <- 1
			cstart[e] <- cin[1]
			for(j in 2:length(cin)){
				if(cin[j-1]-cin[j]!=-1){
					cstop[e] <- cin[j-1]
					e <- e+1
					cstart[e] <- cin[j]
				}
			}
			cstop[e] <- cin[j]
			
			for(j in 1:length(rstart)){
				for(k in 1:length(cstart)){
					xl <- cstart[k] - 0.5
					yb <- nr - rstop[j] + .5
					xr <- cstop[k]+0.5
					yt <- nr - rstart[j] + 1.5
					rect(xleft=xl,ybottom=yb,xright=xr,ytop=yt,density=0,angle=45*i,lwd=2.5,col=clustercols[i])
				}
			} 
		}
	}
	min.raw <- min(X)
	max.raw <- max(X)
	z <- seq(min.raw, max.raw, length=length(heatcols))
	image(z=t(matrix(z, ncol=1)),col=heatcols, 
			xaxt="n", yaxt="n")
	axis(4,at=seq(0,1,by=.5),labels=c(round(min.raw,digits=1),0,round(max.raw,digits=1)))
}

# stabpathplot

stabpath <- function(res,number){
	#if(!class(res@Parameters$Method)=="BCs4vd"){
	#	stop("object is not of class BCs4vd")
	#}
	if(res@Parameters$savepath&!res@Parameters$pointwise){
		vc <- res@info[[number]][[1]]
		uc <- res@info[[number]][[2]]
		par(mfrow=c(1,2),omi=c(0.25, 0.25, 0.5, 0.25)) #c(bottom, left, top, right)
		n <- length(vc$qs)-1
		p <- length(uc$qs)-1
		pcerv <- res@Parameters$pcerv
		lv <- res@info[[number]][[1]]$l 
		lu <- res@info[[number]][[2]]$l
		thrv <- ((vc[[4]]^2/((n*pcerv)*n))+1)/2
		redv <- which(vc[[7]][,lv]>thrv[lv])
		pceru <- res@Parameters$pceru
		thru <- ((uc[[4]]^2/((p*pceru)*p))+1)/2
		redu <- which(uc[[7]][,lu]>thru[lu])
		colsv <- rep("black",n)
		colsu <- rep("black",p)
		colsv[redv] <- "red"
		colsu[redu] <- "red"
		matplot(t(vc[[7]]),type="l",col=colsv,lty=1,ylab="selection probability",xlab=expression(paste(lambda[v])),main="stability path columns",ylim=c(0,1))
		abline(v=lv,col="darkred",lwd=2)
		abline(h=thrv[lv],col="darkred")
		lines(thrv,col="darkred",lwd=3,lty=2)
		legend(-(n*0.15), 1.05, c(paste("PCER ",pcerv)),text.col = c("darkred"),bty="n")
		matplot(t(uc[[7]]),type="l",col=colsu,lty=1,ylab="selection probability",xlab=expression(paste(lambda[u])),main="stability path rows",ylim=c(0,1))
		abline(v=lu,col="darkred",lwd=2)
		abline(h=thru[lu],col="darkred")
		lines(thru,col="darkred",lwd=3,lty=2)
		legend(-(p*0.15), 1.05, c(paste("PCER ",pceru)),text.col = c("darkred"),bty="n")
		title(paste("Stability Paths Bicluster: ",number),outer=T)
	}else{
		vc <- res@info[[number]][[1]]
		uc <- res@info[[number]][[2]]
		par(mfrow=c(1,2),omi=c(0.25, 0.25, 0.5, 0.25)) #c(bottom, left, top, right)
		n <- length(vc$qs)-1
		p <- length(uc$qs)-1
		pcerv <- res@Parameters$pcerv
		pceru <- res@Parameters$pceru
		lv <- res@info[[number]][[1]]$l 
		lu <- res@info[[number]][[2]]$l
		redv <- which(vc[[2]]>=res@info[[number]][[1]]$thr)
		redu <- which(uc[[2]]>=res@info[[number]][[2]]$thr)
		colsv <- rep("black",n)
		colsu <- rep("black",p)
		colsv[redv] <- "red"
		colsu[redu] <- "red"
		plot(vc[[2]],col=colsv,cex=.5,ylim=c(0,1),ylab="selection probability",main="stability path columns",xlab="columns")
		abline(h=res@info[[number]][[1]]$thr,col="darkred")
		legend(-(n*0.15), 1.05, c(paste("PCER ",pcerv)),text.col = c("darkred"),bty="n")
		plot(uc[[2]],col=colsu,cex=.5,ylim=c(0,1),ylab="selection probability",main="stability path rows",xlab="rows")
		abline(h=res@info[[number]][[2]]$thr,col="darkred")
		legend(-(p*0.15), 1.05, c(paste("PCER ",pceru)),text.col = c("darkred"),bty="n")
	}
}
