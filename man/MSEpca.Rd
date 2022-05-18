\name{MSEpca}
\alias{MSEpca}
\title{MSE on PCA}
\description{Caculate the MSE value on PCA}
\usage{
MSEpca(V = V, X = X, n = n, p = p, m = m, K = K, L = L)
}
\arguments{
  \item{V}{is the right singular matrix}
  \item{X}{is the orignal data set}
  \item{n}{is the sample size}
  \item{p}{is the number of variables}
  \item{m}{is the number of eigenvalues}
  \item{K}{is the number of nodes}
  \item{L}{is the number of subgroups}
}

\value{
 \item{MSEpca}{the MSE value on PCA}
}

\examples{
data(Application)
X=Application
n=nrow(Application);p=ncol(Application)
m=5;L=4;K=4
DLPCA_result=DLPCA(X=X,n=n,p=p,m=m,K=K,L=L)                  
V=DLPCA_result$V
MSEpca_result=MSEpca(V=V,X=X,n=n,p=p,m=m,K=K,L=L) 
MSE_PCA=MSEpca_result$MSEpca
}
