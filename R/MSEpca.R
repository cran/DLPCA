#' Caculate the MSE values of the DLPCA method
#'
#' @param V is the right singular matrix
#' @param X is the orignal data set
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

MSEpca <-
function(V=V,X=X,n=n,p=p,m=m,K=K,L=L){
  B=V%*%t(V)%*%t(X)
  U=NULL;Q=NULL;MSEE=NULL
  for(i in 1:n){
    for(j in 1:p){
      E=t(V[,(m+1):p])%*%B[,i]
      U[j]=(V[j,(m+1):p]%*%E)^2    
    }
    Q[i]=sum(U)                                
  }
  L=rep(1,1,n); 
  nk=ceiling(n/K)
  for (k in 1:K){
    SSE=sum(Q[(((k-1)*nk)+(1:nk))])
    MSE=(1/(p*nk))*SSE
    MSEE[k]=MSE
  }
  MSEpca=min(MSEE)
  return(list(MSEpca=MSEpca))
}
