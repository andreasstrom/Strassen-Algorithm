/*
Header file for strassen Algorithm
Author: Andreas Ström
*/

#define TRESH 128 //found for my system by trial and error

//will run on 1, 7 or 49 threads depending on input [<7, ==7 or >7 respectively]
#define THREADS 49

//structs used for parallelization
typedef struct p_strass{
    int* dst;
    int* a;
    int* b;
    int size;
} p_strass_t;

typedef struct p_strass_init{
    int* dst;
    int* a;
    int* b;
    int threads;
    int size;
} p_strass_init_t;

void fill_matrix(int* theMatrix, int n, int min, int max);
void print_matrix(int* theMatrix, int n);
void* strass(void* arguments);
void* strass_init(void* arguments);
int verify_result(int*a, int* b, int* c, int* d, int size);