// works but results of test statistic seems not correct 
extern "C"  void gpu_normBlocks(double *Genes, double *Select, double *C, int *nA, int *N, int *M, int *K, int *R); 

__global__ void
normalization(const double *Genes, const double *Select, double *C, int nA, int N, int M, int K, int R)
{
 int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i < R)
  {
   int indk0=0;
   double tmp = 0;
   for(int k=0; k<K; k++){
  	    indk0 = int(Select[k+i*K]-1);
  		for(int n=0; n<N; n++){
  			tmp+=Genes[n+indk0*N];
  		}
   }
   
   int indt = 0; //index of tested gene 
   int indk = 0; //index of normalizing gene


   double shat = 0; 
   double u = 0; double sum1=0; double sum2=0;
   double meanA=0; double varA = 0;
   double meanB=0; double varB = 0;



int val = 0;

   for(int x=1; x<(M+1); x++){ 
   		val = 0; // val =0 if normalizing and val = 1 if to be tested
   		for(int k = 0; k<K; k++){ 
   			if(Select[k+i*K] == x){val=val+1;break;}
   		}
   		if(val==0){//if gene 'x' to be tested
   			indt = int(x-1); //index of t-th genes tested
   		
   			sum1=0; sum2=0;
   			for(int n=0; n<nA; n++){ // samples or individuals

   			 shat = 0;
   				for(int k=0; k<K; k++){ //normalization for the n-th subject
   					indk = int(Select[k+i*K]-1);
  					shat+=Genes[n+indk*N];
   				}
   				
   				shat = shat/tmp*N;
   				u = 2*pow(Genes[n+indt*N]/shat,0.5); 
   				sum1 += u; 
   				sum2 += u*u;
   			}
   			meanA = sum1/nA;
   			varA = (sum2 -sum1*sum1/nA)/(nA-1);
   		
   			sum1=0; sum2=0;
   			for(int n=nA; n<N; n++){ // samples or individuals
   			shat = 0;
   				for(int k=0; k<K; k++){ //normalization for the n-th subject
   					indk = int(Select[k+i*K]-1);
  					shat+=Genes[n+indk*N];
   				}
   				shat = shat/tmp*N;
   				u = 2*pow(Genes[n+indt*N]/shat,0.5); 
   				sum1 += u; 
   				sum2 += u*u;
   			}
   			meanB = sum1/(N-nA);
   			varB = (sum2 -sum1*sum1/(N-nA))/(N-nA-1);
   	
   		
   		   C[indt+i*M] = (meanA-meanB)*(meanA-meanB)/(varA/nA+varB/(N-nA));
   		
   		}//end if.
   }//end for(int x...
   
   
   
 
   
   
  }
}


void gpu_normBlocks(double *Genes, double *Select, double *C, int *nA, int *N, int *M, int *K, int *R) 
{
  // Device Memory
  
  double *d_Genes, *d_Select, *d_C, *tmp, *shat, *u, *sum1, *sum2, *meanA, *meanB, *varA, *varB;
  
  // Define the execution configuration
  double THREADS = 1;
  
  double n_blocksx = ceil(*R/THREADS); 
  dim3 threadPerBlock(THREADS);
  dim3 numBlocks(n_blocksx);
    
  // Allocate output array
  cudaMalloc((void**)&d_Genes, *N * *M * sizeof(double));
  cudaMalloc((void**)&d_Select, *K * *R * sizeof(double));
  cudaMalloc((void**)&d_C, *M * *R * sizeof(double));
  cudaMalloc((void**)&tmp, sizeof(double));
  cudaMalloc((void**)&shat, sizeof(double));
  cudaMalloc((void**)&u, sizeof(double)); cudaMalloc((void**)&sum1, sizeof(double)); cudaMalloc((void**)&sum2, sizeof(double));
  cudaMalloc((void**)&meanA, sizeof(double)); cudaMalloc((void**)&meanB, sizeof(double));
  cudaMalloc((void**)&varA, sizeof(double)); cudaMalloc((void**)&varB, sizeof(double));
  
  // copy data to device
  cudaMemcpy(d_Genes, Genes, *N * *M * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Select, Select, *K * *R * sizeof(double), cudaMemcpyHostToDevice);
  
  // GPU vector normalization
  normalization<<<numBlocks,threadPerBlock>>>(d_Genes, d_Select, d_C, *nA, *N, *M, *K, *R);
  
  // Copy output
  cudaMemcpy(C, d_C, *M * *R * sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(d_Genes);
  cudaFree(d_Select);
  cudaFree(d_C);
  cudaFree(tmp);  cudaFree(shat);
  cudaFree(u); cudaFree(sum1); cudaFree(sum2);
  cudaFree(meanA);cudaFree(meanB);cudaFree(varA);cudaFree(varB);
}

//Bild shared object 
//https://forums.developer.nvidia.com/t/shared-library-creation/4776/8
//
// nvcc --ptxas-options=-v --compiler-options '-fPIC' -o gpu_normBlocks.so --shared gpu_normBlocks.cu
