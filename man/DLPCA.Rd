\name{DLPCA}
\alias{DLPCA}
\title{Distributed local PCA}
\description{Calculate the estimator on the DLPCA method}
\usage{
DLPCA(X = X, n = n, p = p, m = m, K = K, L = L)
}
\arguments{
  \item{X}{is the original data matrix}
  \item{n}{is the sample size}
  \item{p}{is the number of variables}
  \item{m}{is the number of eigenvalues}
  \item{K}{is the number of nodes}
  \item{L}{is the number of subgroups}
}

\value{
  \item{time}{is the time cost}
  \item{V}{is the right singular matrix}
  \item{Vm}{is the m-right singular matrix}
  \item{Smean}{is the mean covariance matrix}
  \item{MMSER}{is the mean MSE values of the robust covariance matrix sub-estimators}
  \item{MMSES}{is the mean MSE values of the covariance matrix sub-estimators}
  \item{MMSEX}{is the mean MSE values of the sub-estimators of the matrix X }
  \item{MSER}{is the min MSE values of the robust covariance matrix sub-estimators}
  \item{MSES}{is the min MSE values of the covariance matrix sub-estimators}
  \item{MSEX}{is the min MSE values  of the sub-estimators of the  matrix  X}
  \item{wMSER}{is the location of the min MSE values of the robust covariance matrix sub-estimators}
  \item{wMSES}{is the location of the min MSE values of the covariance matrix sub-estimators}
  \item{wMSEX}{is the location of the min MSE values of the sub-estimators of the matrix X}
  \item{sigm}{is the estimator of the covariance matrix of the matrix X}
}
\examples{
data(Application)
X=Application
n=nrow(Application);p=ncol(Application)
m=5;L=4;K=4
DLPCA_result=DLPCA(X=X,n=n,p=p,m=m,K=K,L=L)        
}

