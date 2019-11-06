#' @import mclust clue
#' @importFrom clusteval sim_normal
#' @title Random Projection Ensemble Clustering Algorithm
#'
#' @description This function allows to run the RPEClu algorithm.
#'
#' @param x A numeric high-dimensional matrix where rows correspond to observations and columns correspond to variables.
#' @param true.cl A vector containing the true cluster membership. If supplied, the Adjusted Rand Index (ARI) of the predicted clustering is also returned. By default is set to NULL.
#' @param g The number of clusters.
#' @param d The dimension of the projected space. If is \code{NULL} (default option), then \eqn{d = \lceil c* log(g) \rceil.}{d = ceil(c* log(g)).}
#' @param c The constant which governs the dimension of the projected space if \eqn{d} is not provided. The default is set to 10.
#' @param B The number of generated projections; the default is 1000.
#' @param B.star The number of \emph{base} models to retain in the final ensemble; the default is 100.
#' @param modelNames A vector of character strings indicating the models to be fitted in the EM phase of the GMM. The help file for \link[mclust]{Mclust} describes the available options.
#' @param seed A single value indicating the initializing seed for random generations.
#' @param diagonal A logical value indicating whether the conditional covariate matrix has a restricted form, i.e. it is diagonal. The default is FALSE.
#' @param ensmethod A character string specifying the method for computing the clustering consensus. See the \link[clue]{cl_consensus} help file for available options. The default is \code{DWH}.
#' @param verb A logical controlling if the progress number of projections is displayed during the fitting procedure. By default is \code{FALSE}.
#' @return The output components are as follows:
#' \item{\code{ensemble}}{A list including:
#' \describe{
#'  \item{\code{label.vec}}{The vector of labels predicted by the ensemble of size \code{B.star}.}
#'  \item{\code{ari}}{The corresponding ARI (if \code{true.cl} is not \code{NULL}).}}}
#' \item{individual}{A list including:
#'  \describe{\item{\code{label.vec}}{The vectors of labels predicted by each \emph{base} model.}
#'   \item{\code{ari}}{The corresponding ARI (if \code{true.cl} is not \code{NULL}).}
#'   \item{\code{bic}}{The BIC associated to each \emph{base} model.}
#'   \item{\code{bic.GMM}}{The BIC associated to the Gaussian mixture fitted on each projected data.}
#'   \item{\code{bic.reg}}{The BIC for the linear regression of the \eqn{(p-d)} last columns of \eqn{Y*} on the first \eqn{d} ones.}}}
#' @references Anderlucci, Fortunato, Montanari (2019) <arXiv:1909.10832>
#' @examples
#' \donttest{
#' data(Meat)
#' out.clu <- RPGMMClu(Meat$x, Meat$y, g=5, B=1000, B.star=100, verb=TRUE)}
#'
#' data <- sim_normal(n = rep(100, 2), p = 100, rho = rep(0.1, 2), delta = 0.5, sigma2 = 1, seed = 106)
#' out.clu <- RPGMMClu(data$x, data$y, g=2, B=10, B.star=5, verb=TRUE)
#'
#' @export

RPGMMClu<-function(x, true.cl=NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, modelNames = NULL, diagonal = FALSE, ensmethod="DWH", seed = 101, verb = FALSE){

  p<-ncol(x)
  if(is.null(d)) d<-ceiling(c*log(g)) else d = d
  n<-nrow(x)

  set.seed(seed)
  RPbase<-generateRP(p=ncol(x),d=d,B=B)
  index=matrix(1:(d*B),d,B)

  Bic<-Ari<-bic.c<-bic.nc<-loglik.c<-det.s.ybar.y<-NULL
  cl.m<-matrix(NA, n, B)

  for (b in 1:B){
    A<-as.matrix(RPbase[,index[,b]])
    A.bar <- qr.Q(qr(A),complete=TRUE)[,(d+1):p]
    y<-x%*%A
    y.bar<-x%*%A.bar
    y.star<-cbind(y,y.bar)
    out<-Mclust(y,g, verbose=FALSE, modelNames = modelNames)
    loglik.c[b]<-out$loglik
    cl.m[,b]<-out$cl
    bic.c[b]<-out$bic # BIC for the d "relevant variables", i.e. the d RPs

    #regression model of the (p-d) not relevant variables on the d relevant ones
    X<-cbind(matrix(1,n),y)
    Y<-y.bar
    B<-solve(t(X)%*%X)%*%t(X)%*%Y
    s.ybar.y<-1/(n-1)*t(Y)%*%(diag(n)-X%*%solve(t(X)%*%X, tol=1e-30)%*%t(X))%*%Y
    if(diagonal) s.ybar.y<-diag(diag(s.ybar.y))
    m<-ncol(y.bar)
    if(diagonal) det.o<-prod(diag(s.ybar.y)) else det.o<-det(s.ybar.y)
    if((det.o)<1e-323) det.s.ybar.y[b]<-1e-323 else if(det.o==Inf) det.s.ybar.y[b]<-1e+308 else det.s.ybar.y[b]<-det.o

    loglik.nc<-(-(m*n)/2)*log((2*pi))-(n/2)*log(det.s.ybar.y[b])-0.5*sum(diag(((Y-X%*%B)%*%solve(s.ybar.y,tol=1e-30)%*%t(Y-X%*%B))))

    # BIC for the (p-d) "not relevant variables"
    if (diagonal) k<-m*(d+1)+m else k<-m*(d+1)+(m*(m+1))/2  # free parameters for the regression model
    bic.nc[b]<-2*loglik.nc-k*log(n)


    Bic[b]<-bic.c[b]+bic.nc[b]
    if(!is.null(true.cl)) Ari[b]=adjustedRandIndex(out$classification,true.cl)

    if(verb) print(b)
  }

  # Ensemble ####

  cl.ens.1<-data.frame(cl.m[,order(Bic,decreasing=TRUE)[1:B.star]])
  cl.ens.1.2<-lapply(cl.ens.1,function(x) as.cl_membership(x))
  cl.consensus<-apply(cl_consensus(cl.ens.1.2,method=ensmethod)$.Data,1,which.max)
  if(!is.null(true.cl)) ari<-adjustedRandIndex(cl.consensus,true.cl)


  names(ari)<-paste0("B.star=",B.star)
  ensemble<-list(ari = ari, label.vec = cl.consensus)
  individual<-list(label.vec = cl.m, ari = Ari, bic = Bic, bic.GMM=bic.c, bic.reg=bic.nc)
  return(list(ensemble = ensemble, individual = individual))
}
