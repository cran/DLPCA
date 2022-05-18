\name{DLPCA}
\alias{DLPCA}
\title{Distributed local PCA}
\description{Caculate the estimator on the DLPCA method}
\usage{
DLPCA(X = X, n = n, p = p, m = m, K = K, L = L)
}
\arguments{
  \item{X}{is the orignal data set}
  \item{n}{is the sample size}
  \item{p}{is the number of variables}
  \item{m}{is the number of eigenvalues}
  \item{K}{is the number of nodes}
  \item{L}{is the number of subgroups}
}

\value{
  \item{time}{is the time costs}
  \item{V}{is the right singular matrix}
  \item{Vm}{is the m-right singular matrix}
  \item{Smean}{is the mean covariance matrix}
  \item{MMSES}{is the mean MSE values on the covariance matrix}
  \item{MMSEX}{is the mean MSE values on X}
  \item{MSES}{is the min MSE values on the covariance matrix}
  \item{MSEX}{is the min MSE values  on X}
  \item{wMSES}{is the location of the min MSE values on the covariance matrix}
  \item{wMSEX}{is the location of the min MSE values on X}
}
\examples{
data(Application)
X=Application
n=nrow(Application);p=ncol(Application)
m=5;L=4;K=4
DLPCA_result=DLPCA(X=X,n=n,p=p,m=m,K=K,L=L)        
}

