#include <cstdio>
#include <omp.h>
int main()
{
    printf("Hello world!\n");
#pragma omp parallel
    {
        int tot = omp_get_num_threads();
        int tid = omp_get_thread_num();

        printf("Total active threads: %i, running on thread %i.\n", tot, tid);
    }
    return 0;
}
