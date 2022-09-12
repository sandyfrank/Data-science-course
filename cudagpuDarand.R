

#' Do Differential Analysis with Random Reference Genes
#'
#' Implement the `cudagpuDArand` procedure for transcriptomic data with GPU. This is GPU version of `DArand` function.
#' The procedure is based on random and repeated selection of subsets of reference genes as described in the paper cited below. Observed counts data are normalized with counts from the subset and a differential analysis is used to detect differentially expressed genes. Thought repetitions, the number times a gene is detected is recorded and the final selection is determined from p-values computed under Binomial distribution and adjusted with the Holm's correction.
#'
#' @param X a two-dimensional array (or data.frame) containing the expression table of n individuals in rows and m gene expressions in columns.
#' @param nA integer, number of individuals of the first category, should be smaller than n
#' @param k integer, number of random genes selected (default `k = ceiling(log2(m))`) as reference genes.
#' @param alpha numeric, global test level (default 0.05)
#' @param eta numeric, inner test level (default 0.05)
#' @param beta numeric, inner type II error (default 0.1)
#' @param r integer, number of random 'reference' set selected (the default 1000)
#' @param with.info logical, if `TRUE` results are displayed (the default `FALSE`)
#' @param clog numeric, constant (default 1) controlling the gaussian approximation of the test statistic (in Negative Binomial and Poisson case) .
#' @param use.multi.core logical, if `TRUE` (the default) parallel computing with `mclapply` is used.
#' @param step integer, only used when use.Iter is TRUE to get information on the number of iterations (default 0). Not for use.
#' @param scales numeric, only used for simulation of oracle purpose (default `NULL`). Not for use.
#' @param use.Iter logical, applies iterative procedure  (default FALSE)
#' @param set.seed numeric, set random seed (as is in [`set.seed`][set.seed] function for random number generation ), here default is `NULL`.
#'
#' @details The expression table should be organized in such a way that individuals are represented in rows and genes in columns of the `X` array. Furthermore, in the current version, the procedure provides a differential analysis comparing exactly two experimental conditions. Hence, lines from 1 to `nA` should correspond to the first condition and the remaining lines to the second condition.
#'
#' @details In the inner part of the procedure, called further *randomization*, scaling factors are estimated using a normalization subset of `k` genes randomly selected from  all m genes.  These `k` genes are used as reference genes.  The normalized data are compared between the experimental conditions within an approximately gaussian test for Poisson  or negative-binomial counts as proposed in the methodology cited below. For this inner test the type I (`eta`) and the type II (`beta`) errors should be specified, otherwise the default values will be used.  Since true reference genes (*housekeeping genes*) are unknown, the inner part is repeated `r` times.
#'
#' @details Through all `r` randomization, for each gene, the number of detections (*i.e.* the number of randomizations when a given gene is identified as differentially expressed) is collected. For these detection counts, the corresponding p-values are computed under the Binomial distribution. The finale detection uses the p-values and,  owing to Holm's correction, controls FWER at specified level `alpha`.
#'
#' @details  The maximal number of discoveries is limited to Delta -  the parameter that is a function of `eta`, `beta` and the probability of selecting a subset containing at least one differentially expressed gene leading to a wrong normalization (see [`select_prob`][select_prob] ) . If `use.Iter` is TRUE (the default), the maximal number of discoveries is limited (per iteration) to Delta. The procedure is iterated as long as the number of discoveries is equal to the value of Delta computed in the iteration. Starting from step=1, at each iteration the one-type error is  halved `alpha=alpha/2` to ensure the overall test level respects the initial alpha.
#'
#' @details `clog` is a constant that controls gaussian approximation of the test statistic for the count data arising from Negative Binomial or Poisson distribution. The constant should be ajusted to keep the probability `1-5*n^(-clog)` high while shift term `1+sqrt(clog*n)` low.
#'
#' @return position vector of the gene expressions found as differentially expressed.
#' @export
#'
#' @author Sandy Frank 
#' @references \emph{Differential analysis in Transcriptomic: The strengh of randomly picking 'reference' genes}, D. Desaulle, C. Hoffman, B. Hainque and Y. Rozenholc.
#' <https://arxiv.org/abs/2103.09872>
#'
#' @examples
#'
#' L = build_example(m=500,m1=25,nA=6,fold=20,mu0=100,use.scales=FALSE,nb.size=Inf)
#' DArand(L$X,L$nA,alpha=0.05)
#'
#'


library(dplyr)
library(ggplot2)
library(Rtsne)
library(cluster)
library(dbscan)
library(plotly)
library(microbenchmark)
library(tidyr)

library(tidyverse)
library(broom)




dyn.load("gpu_norm10.so")

cudagpuDArand = function(X,nA,k=NULL,alpha=0.05,
                eta=0.05,beta=0.1,r=r,with.info=FALSE,
                clog=1,step=0,scales=NULL, 
                use.Iter=TRUE,set.seed=NULL) {
  
  if(!is.null(set.seed)) set.seed(set.seed)
  if (!(is.data.frame(X) | is.array(X)) ) stop('X must be an array or a data.frame object')
  
  if (nA>=dim(X)[1]) stop('nA is larger than the number of individuals n.')
  if (is.null(k)) k = ceiling(log2(dim(X)[2]))
  
  if (k/dim(X)[2]>0.2) warning(paste0('Number of reference genes: ',k,'  is large compared to the total number of genes: ',dim(X)[2],'.'))
  
  if (with.info & step>=0 & use.Iter) print(paste('Step',step<-step+1,': use',k,'random references'))
  
  if (use.Iter) alpha=alpha/2 # start with alpha/2 to ensure global procedure with level alpha
  
  n = dim(X)[1]; nB = n-nA # population A with size nA
  m = dim(X)[2] # number of genes
  
  # proba for at least one DE gene to be selected
  pi0d = 1-select_prob(m,k) # \pi_d^0
  pi1d = c(0,pi0d[1:(m-k)])
  
  log_correct = (1+sqrt(clog*log(n)))^2
  
  q2_eta = qnorm(1-eta/2)^2
  
  threshold = log_correct * q2_eta
  
  results = matrix(rep(0,m*r), nrow = m, ncol = r) # O if not tested neither detected
  
  # genes = build_example(m,m1=m1,n1=nA,fold=fold)$X
  
  select = replicate(r,{S=sample(m,k,replace=F);sort(S,decreasing=T)})
  
  genes = X ##

  
  norm_D<-function(Genes, Select,sizenA){
    
    N=nrow(Genes)
    M=ncol(Genes)
    
    K=nrow(Select)
    R=ncol(Select)
    CC = matrix(rep(0, M*R), nrow=M, ncol=R) 
    
    res = .C("gpu_norm10", 
             Genes=as.double(Genes),
             Select = as.double(Select),
             C = as.double(CC), 
             nA=as.integer(sizenA), N=N, M=M, K=K, R=R)
    
    res
    
  }
  
  
  
  res = norm_D(genes, select,nA)
  
  Cgpu=matrix(res$C, nrow=m, ncol=r)
  
  tmp = results
  
  tmp = (Cgpu > threshold)
  
  Rj = rowSums(tmp)
  
  # Holm's detections and number of tests
  tri.holm = order(Rj,decreasing=TRUE)
  # non DE gene detection rate when d varies
  error.rate = (1-pi0d)*eta + pi0d # Binom(rj,error.rate)
  error.rate = pmin(error.rate*(m-k)/m, 1) # Binom(r,error.rate)
  
  #if (use.Delta) Delta = sum((1-beta)*(1-pi1d) > eta*(1-pi0d)+pi0d) else Delta = Inf #DD we always need to use Delta is required by the method (validity condition)
  Delta = sum((1-beta)*(1-pi1d) > eta*(1-pi0d)+pi0d)
  
  # next we compute sorted p-values with the good non DE gene detection rate
  pvals.holm = rep(1, l=m)
  for (d in 1:min(m,Delta)) {
    jnd = tri.holm[d]
    pvals.holm[d] = stats::pbinom(Rj[jnd],r,error.rate[d],lower.tail=F)
  }
  
  # now we find the detections after randomization:
  d.hat = 0
  while (pvals.holm[d.hat+1]<alpha/(m-d.hat) & d.hat<min(Delta,m)) d.hat = d.hat+1
  if (d.hat>0) detect.holm = tri.holm[1:d.hat] else detect.holm=integer(0)
  
  # display information
  if (with.info)
    print(paste0(ifelse(use.Iter,'With','Without'),' iterative procedure (delta=',Delta,'): ',d.hat,' founds'))
  
  
  # do iterations when Delta discoveries
  if (use.Iter & (d.hat==Delta)) {
    stay = setdiff(1:m,detect.holm)
    found = cudagpuDArand(X=X[,stay,drop=F],nA=nA,k=k,alpha=alpha,eta=eta,beta=beta,r=r,with.info=with.info,clog=clog,step=step)
    detect.holm = c(detect.holm,stay[found])
  }
  
  # sort detections
  if (d.hat>0) detect.holm=sort(detect.holm)
  
  if (with.info & use.Iter & step<=1) print(paste(length(detect.holm),'detected genes'))
  
  # return detections
  detect.holm
}


library(DArand)


m = 20000; m1 = 200 ;  fold =  0 ; k = NULL
r = 1000 ; nA = nB = 60 ; eta = 0.05 ; 
q2_eta = qnorm(1-eta/2)^2 ; c.log = 2.0

X = build_example(m,m1=m1,n1=nA,fold=fold)$X




### Trying and timing

dyn.load("gpu_norm10.so")

library(DArand)

m = 20;
m1 = 2
fold = 0
k = NULL
r = 10
nA = nB = 60
eta = 0.05
q2_eta = qnorm(1-eta/2)^2
c.log = 2.0


X = build_example(m,m1=m1,n1=nA,fold=fold)$X

n = c(1,2,6)
ff= lapply(n, function(n) {a =microbenchmark(n + n, check='identical')

})

ffsummary <- do.call(rbind,(ff))


n = c(3,6,2)
res = microbenchmark(
  ff= lapply(n, function(n) {a =microbenchmark(n + n, check='identical')
})
)

n = c(4,7,6)
res = lapply(n, function(n){a=microbenchmark(n+n ,unit = NULL, times = 2L)

})
print(unlist(res), unit= 's')



tidy.microbenchmark <- function(x, unit, ...){
  summary(x, unit = unit)
}

res_tidy = tidy(res) %>% 
  mutate(expr = as.character(expr)) %>% 
  separate(expr, c("func","n"), remove = FALSE)

res_tidy

ggplot(res_tidy, aes(x = n, y = mean, group = func, col = func)) +
  geom_line() +
  geom_point() +
  labs(y = "Runtime", x = "n")









microbenchmark(cudagpuDArand(X,nA,k,alpha=0.05,eta=0.05,beta=0.1,r,with.info=FALSE,clog=1,step=0,scales=NULL, use.Iter=TRUE,set.seed=NULL), times = 2L)



microbenchmark(DArand(X,nA,k,alpha=0.05,eta=0.05,beta=0.1,r,with.info=FALSE, clog=1,step=0,scales=NULL, use.Iter=TRUE,set.seed=NULL),times = 2L)




cudagpuDArand(X,nA,k,alpha=0.05,eta=0.05,beta=0.1,r,with.info=FALSE, clog=1,step=0,scales=NULL, use.Iter=TRUE,set.seed=NULL)



DArand(X,nA,k,alpha=0.05,eta=0.05,beta=0.1,r,with.info=FALSE,clog=1, step=0,scales=NULL, use.Iter=TRUE,set.seed=NULL)







DArand1 = function(X, n1, k = NULL, alpha = 0.05, eta = 0.05, beta = 0.1, 
          r = 1000, with.info = FALSE, clog = 1, use.multi.core = TRUE, 
          step = 0, scales = NULL, use.Iter = TRUE, set.seed = NULL) 
{
  if (!is.null(set.seed)) 
    set.seed(set.seed)
  if (!(is.data.frame(X) | is.array(X))) 
    stop("X must be an array or a data.frame object")
  if (n1 >= dim(X)[1]) 
    stop("n1 is larger than the number of individuals n.")
  if (is.null(k)) 
    k = ceiling(log2(dim(X)[2]))
  if (k/dim(X)[2] > 0.2) 
    warning(paste0("Number of reference genes: ", k, 
                   "  is large compared to the total number of genes: ", 
                   dim(X)[2], "."))
  if (with.info & step >= 0 & use.Iter) 
    print(paste("Step", step <- step + 1, ": use", 
                k, "random references"))
  if (use.Iter) 
    alpha = alpha/2
  n = dim(X)[1]
  n2 = n - n1
  m = dim(X)[2]
  pi0d = 1 - select_prob(m, k)
  pi1d = c(0, pi0d[1:(m - k)])
  Rj = rj = rep(0, m)
  do.randomization = function(i) {
    select = sample(seq_len(m), k)
    tested = setdiff(seq_len(m), select)
    if (is.null(scales)) 
      s.hat = n * rowSums(X[, select, drop = FALSE])/sum(X[, 
                                                           select, drop = FALSE])
    else s.hat = scales
    Y = 2 * sqrt(X[, tested, drop = F]/s.hat)
    Y1.bar = colMeans(Y[1:n1, , drop = F])
    Y2.bar = colMeans(Y[-(1:n1), , drop = F])
    S1 = n1 * rowMeans((t(Y[1:n1, , drop = F]) - Y1.bar)^2)/(n1 - 
                                                               1)
    S2 = n2 * rowMeans((t(Y[-(1:n1), , drop = F]) - Y2.bar)^2)/(n2 - 
                                                                  1)
    pvals = 2 * stats::pnorm(abs(Y1.bar - Y2.bar)/sqrt(S1/n1 + 
                                                         S2/n2)/(1 + sqrt(clog * log(n))), lower.tail = FALSE)
    detect = which(pvals < eta)
    if (use.multi.core) {
      detect = tested[detect]
      return(list(select = select, detect = detect))
    }
    if (!use.multi.core) {
      if (length(detect) > 0) 
        Rj[tested[detect]] <<- Rj[tested[detect]] + 1
      return(NULL)
    }
  }
  if (!use.multi.core) 
    tmp = lapply(1:r, do.randomization)
  if (use.multi.core) {
    tmp = parallel::mclapply(1:r, do.randomization)
    lapply(tmp, function(l) {
      if (length(l$detect) > 0) 
        Rj[l$detect] <<- Rj[l$detect] + 1
    })
  }
  tri.holm = order(Rj, decreasing = TRUE)
  error.rate = (1 - pi0d) * eta + pi0d
  error.rate = pmin(error.rate * (m - k)/m, 1)
  Delta = sum((1 - beta) * (1 - pi1d) > eta * (1 - pi0d) + 
                pi0d)
  pvals.holm = rep(1, l = m)
  for (d in 1:min(m, Delta)) {
    jnd = tri.holm[d]
    pvals.holm[d] = stats::pbinom(Rj[jnd], r, error.rate[d], 
                                  lower.tail = F)
  }
  d.hat = 0
  while (pvals.holm[d.hat + 1] < alpha/(m - d.hat) & d.hat < 
         min(Delta, m)) d.hat = d.hat + 1
  if (d.hat > 0) 
    detect.holm = tri.holm[1:d.hat]
  else detect.holm = integer(0)
  if (with.info) 
    print(paste0(ifelse(use.Iter, "With", "Without"), 
                 " iterative procedure (delta=", Delta, "): ", 
                 d.hat, " founds"))
  if (use.Iter & (d.hat == Delta)) {
    stay = setdiff(1:m, detect.holm)
    found = DArand(X = X[, stay, drop = F], n1 = n1, k = k, 
                   alpha = alpha, eta = eta, beta = beta, r = r, with.info = with.info, 
                   clog = clog, step = step)
    detect.holm = c(detect.holm, stay[found])
  }
  if (d.hat > 0) 
    detect.holm = sort(detect.holm)
  if (with.info & use.Iter & step <= 1) 
    print(paste(length(detect.holm), "detected genes"))
  detect.holm
}

