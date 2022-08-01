

#' Do Differential Analysis with Random Reference Genes
#'
#' Implement the `DArand` procedure for transcriptomic data.
#' The procedure is based on random and repeated selection of subsets of reference genes as described in the paper cited below. Observed counts data are normalized with counts from the subset and a differential analysis is used to detect differentially expressed genes. Thought repetitions, the number times a gene is detected is recorded and the final selection is determined from p-values computed under Binomial distribution and adjusted with the Holm's correction.
#'
#' @param X a two-dimensional array (or data.frame) containing the expression table of n individuals in rows and m gene expressions in columns.
#' @param G a grouping variable, shoud be a factor
#' @param n1 integer, number of individuals of the first category, should be smaller than n
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
#' @param n_vec A vector for number of individuals for each categories 
#' @details The expression table should be organized in such a way that individuals are represented in rows and genes in columns of the `X` array. Furthermore, in the current version, the procedure provides a differential analysis comparing exactly two experimental conditions. Hence, lines from 1 to `n1` should correspond to the first condition and the remaining lines to the second condition.
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
#' L = build_example(m=500,m1=25,n1=6,fold=20,mu0=100,use.scales=FALSE,nb.size=Inf)
#' DArand(L$X,L$n1,alpha=0.05)
#'
#'
#'



LMrand = function(X, crtl, n_vec, k=NULL,alpha=0.05,eta=0.05,beta=0.1,r,with.info=FALSE,clog=1,use.multi.core=TRUE,step=0,scales=NULL, use.Iter=TRUE,set.seed=NULL) {
  
  if(!is.null(set.seed)) set.seed(set.seed)
  if (!(is.data.frame(X) | is.array(X)) ) stop('X must be an array or a data.frame object')
  
  #if (n1>=dim(X)[1]) stop('n1 is larger than the number of individuals n.')
  if (is.null(k)) k = ceiling(log2(dim(X)[2]))
  
  if (k/dim(X)[2]>0.2) warning(paste0('Number of reference genes: ',k,'  is large compared to the total number of genes: ',dim(X)[2],'.'))
  
  if (with.info & step>=0 & use.Iter) print(paste('Step',step<-step+1,': use',k,'random references'))
  
  if (use.Iter) alpha=alpha/2 # start with alpha/2 to ensure global procedure with level alpha
  
  n = dim(X)[1]
  m = dim(X)[2] # number of genes
  
  
  G = unlist(sapply(1:length(n_vec), function(l) factor(rep(LETTERS[l],n_vec[l]))))
  G = relevel(factor(G), ref = crtl) ##DD ref = Ctrl
  
  
  ##DD counts for modalities of G variable w
  
  # proba for at least one DE gene to be selected
  pi0d = 1-select_prob(m,k) # \pi_d^0
  pi1d = c(0,pi0d[1:(m-k)])
  
  # Holm's detections and number of tests
  
  Rj = lapply(1:(length(n_vec)-1), function(i) rep(0,m)) ##DD selon le nb de modalités -1 (sans ref)
  
  ########## randomization and detection
  do.randomization = function(i) {
    
    # draw normalization subset of "housekeeping genes"
    select = sample(seq_len(m), k)
    tested = setdiff(seq_len(m), select)
    
    # estimate standarized scales from selected genes oracle
    if (is.null(scales)) s.hat = n*rowSums(X[,select,drop=FALSE])/sum(X[,select,drop=FALSE]) else s.hat = scales
    
    # normalize non selected genes
    Y = 2*sqrt(X[,tested,drop=F]/s.hat)
    
    # linear regression model
    
    lm = sapply(1:dim(Y)[2], function(j) summary(lm(Y[,j]~G)))
    
    pval_Gi = function(k){P = sapply(1:dim(Y)[2], function(j) lm[,j]$coefficients[k+1,4])}
    
    pval_G = lapply( 1:(length(n_vec)-1), function(k)  pval_Gi(k))
      
    
    # les indices des detectés
    
    detect = lapply(1:(length(n_vec)-1), function(i) which(pval_G[[i]]<eta))
    
    
    

    ## les indices de tested detectés 

    if (use.multi.core) {
 
      detect = lapply(1:(length(n_vec)-1), function(i) tested[detect[[i]]])
      return(list(select=select,detect=detect))
    }
    
    if (!use.multi.core) {
      # add 1 to detected genes
      lapply(1:(length(n_vec)-1), function(i) {if (length(detect[[i]])>0){Rj[[i]][tested[detect[[i]]]] <<- Rj[[i]][tested[detect[[i]]]]+1}else {return(NULL)}})
    }
  }
  
  
  
  # run random normalization
  if (!use.multi.core) tmp = lapply(1:r, do.randomization)
  

  
  if (use.multi.core) {
    tmp = parallel::mclapply(1:r,do.randomization) # list of length r returning the two vectors select and detect
    lapply(1:(length(n_vec)-1), function(i){lapply(tmp,function(l) { if (length(l$detect[[i]])>0){ Rj[[i]][l$detect[[i]]] <<- Rj[[i]][l$detect[[i]]]+1}})})
  }
  
  tri.holm = lapply(1:(length(n_vec)-1), function(i) order(Rj[[i]],decreasing=TRUE))
  
  # non DE gene detection rate when d varies
  error.rate = (1-pi0d)*eta + pi0d # Binom(rj,error.rate)
  error.rate = pmin(error.rate*(m-k)/m, 1) # Binom(r,error.rate)
  
  #if (use.Delta) Delta = sum((1-beta)*(1-pi1d) > eta*(1-pi0d)+pi0d) else Delta = Inf #DD we always need to use Delta is required by the method (validity condition)
  Delta = sum((1-beta)*(1-pi1d) > eta*(1-pi0d)+pi0d)
  
  # next we compute sorted p-values with the good non DE gene detection rate

  pvals.holm = lapply(1:(length(n_vec)-1), function(l) rep(1, l=m))
  
  
  
  lapply(1:(length(n_vec)-1), function(l) {for (d in 1:min(m,Delta)) { pvals.holm[[l]][d] <<- stats::pbinom(Rj[[l]][tri.holm[[l]][d]],r,error.rate[d],lower.tail=F) }})
  
 
  
  # now we find the detections after randomization:
  d.hat = rep(0, (length(n_vec)-1) )
  
  lapply(1:(length(n_vec)-1), function(i) {while (pvals.holm[[i]][d.hat[i]+1]< alpha/(m-d.hat[i]) & d.hat[i]< min(Delta,m)) d.hat[i] <<- d.hat[i]+1})
  
  detect.holm = lapply(1:(length(n_vec)-1), function(i) {if (d.hat[i]>0) dh = tri.holm[[i]][1:d.hat[i]] else dh = integer(0)})
  
  
  # display information
  if (with.info)
    print(paste0(ifelse(use.Iter,'With','Without'),' iterative procedure (delta=',Delta,'): ',d.hat,' founds'))
  
  
  # do iterations when Delta discoveries in each group
  

  
   
   lapply(1:(length(n_vec)-1), function(i){if (use.Iter & (d.hat[[i]]==Delta)) {
     stay  = setdiff(1:m,detect.holm[[i]])
     found = unlist(LMrand(X, crtl, n_vec, k=NULL,alpha=0.05,eta=0.05,beta=0.1,r,with.info=FALSE,clog=1,use.multi.core=TRUE,step=0,scales=NULL, use.Iter=TRUE,set.seed=NULL))
     detect.holm[[i]] <<- c(detect.holm[[i]],stay[found])
   }})
  
  
  # sort detections in each group
  
  
  lapply(1:(length(n_vec)-1), function(i){ if (d.hat[[i]]>0) detect.holm[[i]] <<- sort(detect.holm[[i]])})
  
  lapply(1:(length(n_vec)-1), function(i){  if (with.info & use.Iter & step<=1) print(paste(length(detect.holm[[i]]),'detected genes in', G[[i]]))})

  detect.holm
  
}








################################ Trying with exemples############################################################


crtl = "D" # modalité de référence de la variable G
L = build_example(m=1000,m1=50,n1=18,fold=500,mu0=100,use.scales=FALSE,nb.size=Inf)

X = L$X
k=5
scales =NULL


n_vec = c(6,6,6,6,6,6)

#' @param crtl reference category for the variable G


r = 100

crtl = "D"
k = 5
use.multi.core=TRUE
n = dim(X)[1]
m = dim(X)[2]
select = sample(seq_len(m), k)
tested = setdiff(seq_len(m), select)


r=20
LMrand(X, crtl = LETTERS[4] , n_vec, k=NULL,alpha=0.05,eta=0.05,beta=0.1,r=r,with.info=FALSE,clog=1,use.multi.core=TRUE,step=0,scales=NULL, use.Iter=TRUE,set.seed=NULL)
  
alpha=0.05
eta=0.05
beta=0.1

with.info=FALSE
clog=1
use.multi.core=TRUE
step=0
scales=NULL
use.Iter=TRUE
set.seed=NULL


