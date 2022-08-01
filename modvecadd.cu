

// This enable the CUDA code to be call from R ( Wrapper function in R creation)

extern "C" void gvectorAdd(double *A, double *B, double *C, int *n);



// This is kernel : executed on the device

__global__ void
vectorAdd(const double *A, const double *B, double *C, int numElements)
{ double A2 = 0  ; double B2 = 0 ; 
  int i = blockDim.x * blockIdx.x + threadIdx.x;
 
  if(i < numElements)
  { 
    A2   = (A[i]*A[i])/numElements ;
    B2   = (B[i]*B[i])/numElements ;

    C[i] = A2 + B2 ;
  }
}




// main code configuration needed to launch the kernel

void gvectorAdd(double *A, double *B, double *C, int *n) 
{
  // Device Memory

  double *d_A, *d_B, *d_C;


  // Define the execution configuration

  double THREADS = 1024;
  
  double n_blocksx = ceil(*n/THREADS); 

  dim3 threadPerBlock(THREADS);
  dim3 numBlocks(n_blocksx);

  // Allocate memory on the device


  cudaMalloc((void**)&d_A, *n * sizeof(double));
  cudaMalloc((void**)&d_B, *n * sizeof(double));
  cudaMalloc((void**)&d_C, *n * sizeof(double));

  // copy data from host to device

  cudaMemcpy(d_A, A, *n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, B, *n * sizeof(double), cudaMemcpyHostToDevice);

  // Launching the kernel

  vectorAdd<<<numBlocks,threadPerBlock>>>(d_A, d_B, d_C, *n);
  
  // Copy output from device back to the host
  
  cudaMemcpy(C, d_C, *n * sizeof(double), cudaMemcpyDeviceToHost);

  // Free device memory 
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);
}

// Compiling the whole using nvcc + creating the shared object 

// nvcc --ptxas-options=-v --compiler-options '-fPIC' -o modvecadd.so --shared modvecadd.cu





