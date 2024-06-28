#include <cstdio>

#define NX 102400

int main(void)
{
    long vecA[NX];
    long sum, psum, sumex;
    int i;
    
    
    sum = 0.0;
#pragma omp parallel default(none) shared(sum, vecA)
{

    /* Initialization of the vectors */
    #pragma omp for
    for (i = 0; i < NX; i++) {
        vecA[i] = (long) i + 1;
    }

    sum = 0.0;

    /* TODO: Parallelize computation */
    #pragma omp for reduction(+:sum)
    for (i = 0; i < NX; i++) {
        sum += vecA[i];
    }
}
    printf("Sum: %ld\n", sum);

    sumex = (long) NX * (NX + 1) / ((long) 2);
    printf("Arithmetic sum formula (exact): %ld\n", sumex);
    printf("Difference: %li\n", sumex - sum);

    return 0;
}
