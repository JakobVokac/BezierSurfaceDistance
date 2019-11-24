#include <stdio.h>

__global__ void kernel(){

  if (threadIdx.x == 0) printf("I am thread 0 in block %d\n", blockIdx.x);
}

int main(){

  kernel<<<2,1024>>>();
  cudaDeviceSynchronize();
}
