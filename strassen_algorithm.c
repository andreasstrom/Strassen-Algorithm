/*
Main file for matrix multiplication, utilizing the Strassen algorithm
Author: Andreas Ström
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

#include "strassen_help.h"

static double get_wall_seconds(void){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

int main(int argc, char *argv[]){

    if(argc<6){
      printf("Input: mat_size[int>0] print[1/0] verify[1/0] min_val[int] max_val[int]\n");
      printf("Tresh hold for hybrid algorithm and thread count specified in strassen_help.h\n");
      exit(-1);
    }
    
    //initialize
    int size = atoi(argv[1]);
    const int print=atoi(argv[2]), verify=atoi(argv[3]), min=atoi(argv[4]), max=atoi(argv[5]);

    if(min>=max || size<0){printf("mat_size < 0  or  min_val >= max_val\n"); exit(-1);}
    
    if(size & (size-1) || size==1){
      size = pow(2,1+(int)(log(size)/log(2)));
      printf("\nOverwrite mat_size into a power of 2, new mat_size: %dx%d\n", size, size);
    }

    //allocate and fill matrices
    const int n = size;
    const int len = n*n;
    int* A = calloc(len, sizeof(int)); fill_matrix(A, len, min, max);
    int* B = calloc(len, sizeof(int)); fill_matrix(B, len, min, max);
    int* C = calloc(len, sizeof(int));
    
    //start and time the Strassen algorithm
    double t = get_wall_seconds();
    strass_init(&((p_strass_init_t){.size = n, .dst = C, .a = A, .b=B, .threads=THREADS}));
    t = get_wall_seconds() - t;
    
    //print matrices (if specified)
    if(print){
      printf("\nMatrix A:\n");
      print_matrix(A, len);
      printf("\nMatrix B:\n");
      print_matrix(B,len);
      printf("\nMatrix C:\n");
      print_matrix(C,len);
    }
    printf("\nStrassen timing: %10.10fs, Matrix sizes: %dx%d \n", t, size, size);

    //verify results (if specified)
    if(verify){
      printf("\nVerifying results:\n");
      int* D = calloc(len, sizeof(int));
      t = get_wall_seconds();
      if (verify_result(A, B, C, D, n)){
        t = get_wall_seconds() - t;
        printf("results OK, Naive Timing: %10.10fs\n\n", t);
      }else{printf("what MISMATCH in results\n");}
      free(D);
    }
    free(A); free(B); free(C);

    return 0;
}