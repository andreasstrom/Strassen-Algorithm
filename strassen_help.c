/*
File containing functions used by strassen Algorithm
Author: Andreas Ström
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

#include "strassen_help.h"

static int* transpose(int* dst, int* a, int size){
    //func to transpose matrices stored as 1D arrays
    for (int i = 0; i < size; i++){
       for (int j = 0; j < size; j++){
          dst[i*size+j] = a[j*size+i];
       }
    }
    return dst;
}

static int* naive(int* dst, int* a, int* b, int n){
    //func to compute matrix multiplication for "small enough" matrices
    //faster than cache blocking for "small enough" matrices
    int tmp;
    int* tmpMat = calloc(n*n, sizeof(int));
    tmpMat = transpose(tmpMat, b, n);
    for(int i=0;i<n; i++){
        for(int j=0; j<n; j++){
            tmp = 0;
            for(int k=0; k<n; k++){
                tmp += a[i*n+k] * tmpMat[j*n+k];
            }
            
            dst[i*n+j]=tmp;
        }
    }
    free(tmpMat);
    return dst;
}

static int* add(int* dst, int* a, int* b, int size){
    //func for element wise matrix addition
    for(int i=0; i<size; i++){
        dst[i] = a[i] + b[i];
    }
    return dst;
}

static int* sub(int* dst, int* a, int* b, int size){
    //func for element wise matrix subtraction
    for(int i=0; i<size; i++){
        dst[i] = a[i] - b[i];
    }
    return dst;
}

 void* strass_init(void* arguments){
    //main func for strassen algorithm
    p_strass_init_t * args = arguments;
    if(args->size<=TRESH){return naive(args->dst, args->a, args->b, args->size);}
 
    int n = args->size * 0.5;
    int len = n*n;

    //allocate space
    int* a11 = calloc(len, sizeof(int)); int* a12 = calloc(len, sizeof(int));
    int* a21 = calloc(len, sizeof(int)); int* a22 = calloc(len, sizeof(int));
    int* b11 = calloc(len, sizeof(int)); int* b12 = calloc(len, sizeof(int));
    int* b21 = calloc(len, sizeof(int)); int* b22 = calloc(len, sizeof(int));
    int* c11 = calloc(len, sizeof(int)); int* c12 = calloc(len, sizeof(int));
    int* c21 = calloc(len, sizeof(int)); int* c22 = calloc(len, sizeof(int));
    int* m1 = calloc(len, sizeof(int)); int* m2 = calloc(len, sizeof(int)); 
    int* m3 = calloc(len, sizeof(int));
    int* m4 = calloc(len, sizeof(int)); int* m5 = calloc(len, sizeof(int));
    int* m6 = calloc(len, sizeof(int)); int* m7 = calloc(len, sizeof(int));
    int* tmp1 = calloc(len, sizeof(int)); int* tmp2 = calloc(len, sizeof(int)); 

    //split matrices
    int len2 = args->size*args->size * 0.5;
    for(int i=0; i<n; i++){   
        for(int j=0;j<n;j++){
            a11[i*n+j] = args->a[i*args->size+j]; a12[i*n+j] = args->a[i*args->size+j+n];
            a21[i*n+j] = args->a[i*args->size+len2+j]; a22[i*n+j] = args->a[i*args->size+len2+j+n];

            b11[i*n+j] = args->b[i*args->size+j]; b12[i*n+j] = args->b[i*args->size+j+n];
            b21[i*n+j] = args->b[i*args->size+len2+j]; b22[i*n+j] = args->b[i*args->size+len2+j+n];
        }
    }

    //allocate matrices used for parallelization
    int* tmp3 = calloc(len, sizeof(int)); int* tmp4 = calloc(len, sizeof(int));
    int* tmp5 = calloc(len, sizeof(int)); int* tmp6 = calloc(len, sizeof(int));
    int* tmp7 = calloc(len, sizeof(int)); int* tmp8 = calloc(len, sizeof(int));
    int* tmp9 = calloc(len, sizeof(int)); int* tmp10 = calloc(len, sizeof(int));

    //create the threads and m-matrices
    int residual = args->threads%7;
    pthread_t threads[7];
    if(args->threads==7){
        pthread_create(&threads[0], NULL, &strass, &((p_strass_t){.dst = m1,
            .a = add(tmp1, a11, a22, len), .b=add(tmp2, b11, b22, len), .size = n}));
        pthread_create(&threads[1], NULL, &strass, &((p_strass_t){.dst = m2,
            .a = add(tmp3, a21, a22, len), .b=b11, .size = n}));
        pthread_create(&threads[2], NULL, &strass, &((p_strass_t){.dst = m3,
            .a = a11, .b= sub(tmp4, b12, b22, len), .size = n}));
        pthread_create(&threads[3], NULL, &strass, &((p_strass_t){.dst = m4,
            .a = a22, .b= sub(tmp5, b21, b11, len), .size = n}));
        pthread_create(&threads[4], NULL, &strass, &((p_strass_t){.dst = m5,
            .a = add(tmp6, a11, a12, len), .b=b22, .size = n}));
        pthread_create(&threads[5], NULL, &strass, &((p_strass_t){.dst = m6,
            .a = sub(tmp7, a21, a11, len), .b=add(tmp8, b11, b12, len), .size = n}));
        pthread_create(&threads[6], NULL, &strass, &((p_strass_t){.dst = m7,
            .a = sub(tmp9, a12, a22, len), .b=add(tmp10, b21, b22, len), .size = n}));

        for(int i=0; i<7; i++)
                pthread_join(threads[i],NULL);

    }else if(args->threads>7){
        pthread_create(&threads[0], NULL, &strass_init, &((p_strass_init_t){.threads = 7, 
            .dst = m1, .a = add(tmp1, a11, a22, len), .b=add(tmp2, b11, b22, len), .size = n}));
        pthread_create(&threads[1], NULL, &strass_init, &((p_strass_init_t){.threads = 7, 
            .dst = m2, .a = add(tmp3, a21, a22, len), .b=b11, .size = n}));
        pthread_create(&threads[2], NULL, &strass_init, &((p_strass_init_t){.threads = 7, 
            .dst = m3, .a = a11, .b= sub(tmp4, b12, b22, len), .size = n}));
        pthread_create(&threads[3], NULL, &strass_init, &((p_strass_init_t){.threads = 7, 
            .dst = m4, .a = a22, .b= sub(tmp5, b21, b11, len), .size = n}));
        pthread_create(&threads[4], NULL, &strass_init, &((p_strass_init_t){.threads = 7, 
            .dst = m5, .a = add(tmp6, a11, a12, len), .b=b22, .size = n}));
        pthread_create(&threads[5], NULL, &strass_init, &((p_strass_init_t){.threads = 7, 
            .dst = m6, .a = sub(tmp7, a21, a11, len), .b=add(tmp8, b11, b12, len), .size = n}));
        pthread_create(&threads[6], NULL, &strass_init, &((p_strass_init_t){.threads = 7, 
            .dst = m7, .a = sub(tmp9, a12, a22, len), .b=add(tmp10, b21, b22, len), .size = n}));

        for(int i=0; i<7; i++)
                pthread_join(threads[i],NULL);
                
    }else{
        strass(&((p_strass_t){.dst = m1, .a = add(tmp1, a11, a22, len),
            .b=add(tmp2, b11, b22, len), .size = n}));
        strass(&((p_strass_t){.dst = m2, .a = add(tmp1, a21, a22, len), .b=b11, .size = n}));
        strass(&((p_strass_t){.dst = m3, .a = a11, .b= sub(tmp1, b12, b22, len), .size = n}));
        strass(&((p_strass_t){.dst = m4, .a = a22, .b= sub(tmp1, b21, b11, len), .size = n}));
        strass(&((p_strass_t){.dst = m5, .a = add(tmp1, a11, a12, len), .b=b22, .size = n}));
        strass(&((p_strass_t){.dst = m6, .a = sub(tmp1, a21, a11, len),
            .b=add(tmp2, b11, b12, len), .size = n}));
        strass(&((p_strass_t){.dst = m7, .a = sub(tmp1, a12, a22, len),
            .b=add(tmp2, b21, b22, len), .size = n}));
    }
    
    free(tmp3); free(tmp4);
    free(tmp5); free(tmp6);
    free(tmp7); free(tmp8);
    free(tmp9); free(tmp10);

    //last part of the Strassen algorithm, adding the resulting m-matrices
    add(c11, add(tmp1,m1,m4,len), sub(tmp2,m7,m5,len), len);
    add(c12, m3, m5, len);
    add(c21, m2, m4, len);
    add(c22, add(tmp1,m1,m3,len), sub(tmp2,m6,m2,len), len);

    //combine results
    for(int i=0; i<n; i++){
        for(int j=0;j<n;j++){
            args->dst[i*args->size+j] = c11[i*n+j]; args->dst[i*args->size+j+n] = c12[i*n+j];
            args->dst[i*args->size+len2+j]=c21[i*n+j];args->dst[i*args->size+len2+j+n]=c22[i*n+j];
        }
    }
    
    //free space
    free(a11); free(a12); free(a21); free(a22);
    free(b11); free(b12); free(b21); free(b22);
    free(c11); free(c12); free(c21); free(c22);
    free(m1); free(m2); free(m3);
    free(m4); free(m5);
    free(m6); free(m7);
    free(tmp1); free(tmp2);
}

void* strass(void* arguments){
    //main func for strassen algorithm
    p_strass_t * args = arguments;
    if(args->size<=TRESH){return naive(args->dst, args->a, args->b, args->size);}

    int n = args->size * 0.5;
    int len = n*n; 

    //allocate space
    int* a11 = calloc(len, sizeof(int)); int* a12 = calloc(len, sizeof(int));
    int* a21 = calloc(len, sizeof(int)); int* a22 = calloc(len, sizeof(int));
    int* b11 = calloc(len, sizeof(int)); int* b12 = calloc(len, sizeof(int));
    int* b21 = calloc(len, sizeof(int)); int* b22 = calloc(len, sizeof(int));
    int* c11 = calloc(len, sizeof(int)); int* c12 = calloc(len, sizeof(int));
    int* c21 = calloc(len, sizeof(int)); int* c22 = calloc(len, sizeof(int));
    int* m1 = calloc(len, sizeof(int)); int* m2 = calloc(len, sizeof(int)); 
    int* m3 = calloc(len, sizeof(int));
    int* m4 = calloc(len, sizeof(int)); int* m5 = calloc(len, sizeof(int));
    int* m6 = calloc(len, sizeof(int)); int* m7 = calloc(len, sizeof(int));
    int* tmp1 = calloc(len, sizeof(int)); int* tmp2 = calloc(len, sizeof(int)); 

    //split matrices
    int len2 = args->size*args->size * 0.5;
    for(int i=0; i<n; i++){   
        for(int j=0;j<n;j++){
            a11[i*n+j] = args->a[i*args->size+j]; a12[i*n+j] = args->a[i*args->size+j+n];
            a21[i*n+j] = args->a[i*args->size+len2+j]; a22[i*n+j] = args->a[i*args->size+len2+j+n];

            b11[i*n+j] = args->b[i*args->size+j]; b12[i*n+j] = args->b[i*args->size+j+n];
            b21[i*n+j] = args->b[i*args->size+len2+j]; b22[i*n+j] = args->b[i*args->size+len2+j+n];
        }
    }

    //recursively utilize strassen algorithm
    strass(&((p_strass_t){.dst = m1, .a = add(tmp1, a11, a22, len),
        .b=add(tmp2, b11, b22, len), .size = n}));
    strass(&((p_strass_t){.dst = m2, .a = add(tmp1, a21, a22, len), .b=b11, .size = n}));
    strass(&((p_strass_t){.dst = m3, .a = a11, .b= sub(tmp1, b12, b22, len), .size = n}));
    strass(&((p_strass_t){.dst = m4, .a = a22, .b= sub(tmp1, b21, b11, len), .size = n}));
    strass(&((p_strass_t){.dst = m5, .a = add(tmp1, a11, a12, len), .b=b22, .size = n}));
    strass(&((p_strass_t){.dst = m6, .a = sub(tmp1, a21, a11, len),
        .b=add(tmp2, b11, b12, len), .size = n}));
    strass(&((p_strass_t){.dst = m7, .a = sub(tmp1, a12, a22, len),
        .b=add(tmp2, b21, b22, len), .size = n}));

    //last part of the Strassen algorithm, adding the resulting m-matrices
    add(c11, add(tmp1,m1,m4,len), sub(tmp2,m7,m5,len), len);
    add(c12, m3, m5, len);
    add(c21, m2, m4, len);
    add(c22, add(tmp1,m1,m3,len), sub(tmp2,m6,m2,len), len);

    //combine results
    for(int i=0; i<n; i++){
        for(int j=0;j<n;j++){
            args->dst[i*args->size+j] = c11[i*n+j]; args->dst[i*args->size+j+n] = c12[i*n+j];
            args->dst[i*args->size+len2+j]=c21[i*n+j];args->dst[i*args->size+len2+j+n]=c22[i*n+j];
        }
    }
    //free space
    free(a11); free(a12); free(a21); free(a22);
    free(b11); free(b12); free(b21); free(b22);
    free(c11); free(c12); free(c21); free(c22);
    free(m1); free(m2); free(m3);
    free(m4); free(m5);
    free(m6); free(m7);
    free(tmp1); free(tmp2);
}

void fill_matrix(int* theMatrix, int n, int min, int max){
    //func to fill matrices with pseudo-random numbers
    for(int i = 0 ; i < n ; i++)
        theMatrix[i]= min + rand() % (max-min);
}

void print_matrix(int* theMatrix, int n){
    //func for printing matrices
    int n2 = sqrt(n);
    for(int i = 0 ; i < n2; i++){
        for(int j = 0 ; j < n2 ; j++){
	        printf("% 4d " , theMatrix[i*n2+j]);
        }
        putchar('\n');
    }
}

int verify_result(int*a, int* b, int* c, int* d, int size){
    //func for verifying results from strassen algorithm vs naive algorithm
    naive(d,a,b,size);
    for (int i = 0; i < size; i++) {
        if(c[i]!=d[i]){return 0;}
    }
    return 1;
}