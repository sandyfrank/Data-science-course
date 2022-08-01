library(DArand)

m = 50; m1 = 10
fold = 5000
k = 3
r = 25
nA = nB = 60

eta = 0.05
q2_eta = qnorm(1-eta/2)^2
c.log = 2.0
log_correct = (1.0+sqrt(c.log*log(nA+nB)))^2
threshold = log_correct * q2_eta

genes = build_example(m,m1=m1,n1=nA,fold=fold)$X

select = replicate(r,{S=sample(m,k,replace=F);sort(S,decreasing=T)})

dyn.load("gpu_norm10.so")

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

system.time(norm_D(genes, select,nA))


res = norm_D(genes, select,nA)

Cgpu=matrix(res$C, nrow=m, ncol=r)

n=nA+nB;
IndA = 1:nA
IndB = nA+(1:nB)
C0 = NULL


system.time(for(i in 1:r){
  ind = select[,i]## vec
  ref = genes[,ind]## mat
  s_hat = rowSums(ref)/sum(ref)*nrow(genes)
  #s_hat = rowSums(ref)#vec
  
  tested = 1:m
  for (j in ind) {
    if (j < m) tested[j:(m-1)] = tested[(j+1):m]
  }
  tested = tested[1:(m-k)] # NOW it is working properly
  Gi = genes[,tested]<<<
  
  for (lig in 1:n) Gi[lig,] = 2*(Gi[lig,]/s_hat[lig])^(1/2)
  
  GiA = Gi[IndA,]
  GiB = Gi[IndB,]
  
  GAbar = colMeans(GiA)## meanA
  GBbar = colMeans(GiB) ## meanB
  
  for (lig in 1:nA) GiA[lig,] = GiA[lig,]-t(GAbar) # now GiA is centered
  for (lig in 1:nB) GiB[lig,] = GiB[lig,]-t(GBbar) # now GiB is centered
  
  S2A = colSums(GiA^2)/(nA-1) ## varA
  S2B = colSums(GiB^2)/(nB-1) ## varB
  Var = S2A/nA + S2B/nB ## var
  
  #stat2_to_one = (GAbar - GBbar)^2 / Var / threshold
  stat2_to_one = (GAbar - GBbar)^2 / Var ## teststat
  
  
  output=rep(0,m)
  #output[tested] = GAbar ## max(abs(C0-Cgpu)) [1] 2.842171e-14
  #output[tested] = GBbar ## max(abs(C0-Cgpu)) [1] 1.421085e-14
  #output[tested] = S2A ## max(abs(C0-Cgpu)) [1] 8.778756e-12
  #output[tested] = S2B ## max(abs(C0-Cgpu)) [1] 1.615597e-12
  output[tested] = stat2_to_one ##  max(abs(C0-Cgpu)) [1] 4.672256e-08
  
  ##C0=cbind(C0,s_hat)
  C0=cbind(C0,output)
  
})

max(abs(C0-Cgpu))
