#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>
#include <unistd.h>

int mpi_rank_global = -1;

void set_mpi_rank(int rank)
    __attribute__((no_instrument_function));

void __cyg_profile_func_enter(void *this_fn, void *call_site)
    __attribute__((no_instrument_function));

void __cyg_profile_func_exit(void *this_fn, void *call_site)
    __attribute__((no_instrument_function));

void set_mpi_rank(int rank)
{
    mpi_rank_global = rank;
}

void __cyg_profile_func_enter(void *this_fn, void *call_site)
{
    void *buffer[10];
    int nptrs = backtrace(buffer, 10);
    char **strings = backtrace_symbols(buffer, nptrs);
    if (strings != NULL) {
        // Print only the top of the stack (first relevant function)
        printf("[Rank %d] entering \n\t%s\n\t%s\n", mpi_rank_global, strings[1], strings[2]); 
        fflush(stdout);
        free(strings);
    }
}

void __cyg_profile_func_exit(void *this_fn, void *call_site)
{
    void *buffer[10];
    int nptrs = backtrace(buffer, 10);
    char **strings = backtrace_symbols(buffer, nptrs);
    if (strings != NULL) {
        // Print only the top of the stack (first relevant function)
        printf("[Rank %d] exiting \n\t%s\n\t%s\n", mpi_rank_global, strings[1], strings[2]); 
        fflush(stdout);
        free(strings);
    }
}