#' Calculate the MSE values of the DLPCA method
#'
#' @param V is the right singular matrix
#' @param X is the original data matrix
#' @param n is the sample size
#' @param p is the number of variables
#' @param m is the number of eigenvalues
#' @param K is the number of nodes
#' @param L is the number of subgroups
#' @return MSEpca
#' @export
#'
#' @examples 
#' data(Application)
#' X=Application
#' n=nrow(Application);p=ncol(Application)
#' m=5;L=4;K=4
#' DLPCA_result=DLPCA(X=X,n=n,p=p,m=m,K=K,L=L)                  
#' V=DLPCA_result$V
#' MSEpca_result=MSEpca(V=V,X=X,n=n,p=p,m=m,K=K,L=L) 
#' MSE_PCA=MSEpca_result$MSEpca
DLPCA=function(X=X,n=n,p=p,m=m,K=K,L=L){
  nk=ceiling(n/K);nl=ceiling(nk/L);s0=min(nl,p)
  MeanXk=matrix(0,L,p)
  MMSER=MMSES=MMSEX=matrix(rep(0,1*K),ncol=K)  
 lll=0
 time=system.time(while(lll==0){ 
mr=c(sample(1:n,n,replace=F),sample(1:n,K*nk-n,replace=F))
  S1=matrix(0,p,p)
  for (i in 1: K){
    Rr=matrix(0,0,p)
    Rk=matrix(0,L*s0,p)
    mri=mr[(i-1)*nk+(1:(nk))]
    Xik=X[mri,]    
    meanXik=(1/nk)*(matrix(1,1,nk)%*%Xik) 
XCik=Xik-matrix(1,nk,1)%*%meanXik 
Vikhat=t(svd(XCik)$v)  
Vikhatm=Vikhat[,1:m]  
SikhatR=Vikhatm%*%t(Vikhatm) 
    Sik=(1/(nk-1))*t(XCik)%*%(XCik)
    mrr=c(sample(1:nk,nk,replace=F),sample(1:nk,L*nl-nk,replace=F))
    for (l in 1:L){
      mrk=mrr[(l-1)*nl+(1:nl)]
      Xkl=XCik[mrk,] 
      meanXkl=(1/nl)*(matrix(1,1,nl)%*%Xkl)
      Rkl=qr.R(qr(Xkl)) 
      Rk[(s0*(l-1)+(1:s0)),]=Rkl
      MeanXk[l,]=meanXkl
    }
    RR=Rk
    MeanXl=sqrt(nl)*(MeanXk-matrix(1,L,1)%*%meanXik)
    h=L
    c=0
    while (h>1){
      c=1+c
      s1=min(p,nl*(2^(c-1)))
      s2=min(p,nl*(2^(c)))
      if (h/2==ceiling(h/2)){
        h=h/2
        R11=matrix(0,h*s2,p)
        for (j in 1:h){
          R1=rbind(Rk[(((j-1)*s1)+(1:s1)),],Rk[(((h+j-1)*s1)+(1:s1)),])
          R11[(s2*(j-1)+(1:s2)),]=qr.R(qr(R1))
        }
        Rk=R11
      }else{
        h=floor(h/2)
        R11=matrix(0,h*s2,p)
        for (j in 1:h){
          R1=rbind(Rk[(((j-1)*s1)+(1:s1)),],Rk[(((h+j-1)*s1)+(1:s1)),])
          R11[(s2*(j-1)+(1:s2)),]=qr.R(qr(R1))
        }
        Rk=R11
      }
    }
    Rk
    h
    if (nrow(Rr)>0){
      Rk=qr.R(qr(rbind(Rk,Rr)))
    }
    R=qr.R(qr(rbind(MeanXl,Rk)))
    V=svd(R)$v
    Vm=V[,1:m]                
    XCikhat=XCik%*%Vm%*%t(Vm) 
Viktilde=t(svd(XCikhat)$v)  
Viktildem=Viktilde[,1:m]  
SiktildeR=Viktildem%*%t(Viktildem) 
MMSEikR=(1/((nk-1)))*sum(diag(t(SikhatR-SiktildeR)%*%(SikhatR-SiktildeR)))
MMSER[i]=MMSEikR
    MMSEXik=(1/(p*(nk-1)))*sum(diag(t(XCik-XCikhat)%*%(XCik-XCikhat)))
    MMSEX[i]=MMSEXik
    Sikhat=(1/(nk-1))*t(XCikhat)%*%XCikhat
    MMSESik=(1/((nk-1)))*sum(diag(t(Sik-Sikhat)%*%(Sik-Sikhat)))
    MMSES[i]=MMSESik
    S1=S1+Sikhat
  }
  lll=1
  })
  Vsig=svd(X)$v    
  Vmsig=Vsig[,1:m]   
  Xm=X%*%Vmsig%*%t(Vmsig)
  sigm=(1/(n-1))*t(Xm)%*%(Xm)
  Smean=S1/K
  MSES=min(MMSES);MSEX=min(MMSEX) ; MSER=min(MMSER)          
  wMSES=which.min(MMSES);wMSEX=which.min(MMSEX); wMSER=which.min(MMSER)
  return(list(time=time,V=V,Vm=Vm,Smean=Smean,MMSER=MMSER,MMSES= MMSES,MMSEX=MMSEX,MSES=MSES,MSEX=MSEX,MSER=MSER,wMSES=wMSES,wMSEX=wMSEX,wMSER=wMSER, sigm= sigm))
}
