#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <liboptarith/primes.h>
#include "clgrp.h"

int main() {
    // Parameters causing the hang
    long D = -503316492500;
    int h_star = 400;
    int init_pow = 1371;
    long ell = 5;

    printf("Reproducing hang with D=%ld, h_star=%d, init_pow=%d, ell=%ld\n", D, h_star, init_pow, ell);

    // Initialize hash tables (mimicking clgrp_ell.c logic)
    // h_max approximation for table size. 
    // Using a safe upper bound logic or just the value from the code.
    // In code: h_upper_bound(-D_max * ell^4). 
    // Here D is already D_sub (-D_max * ell^4 ish).
    // So h_upper_bound(D) (absolute value)
    
    int h_max = h_upper_bound(D); 
    printf("h_max bound: %d\n", h_max);
    
    int table_size = next_prime((((int) sqrt(h_max)) << 1) - 1);
    printf("Table size: %d\n", table_size);

    htab_t R, Q;
    htab_init(&R, table_size, sizeof(form_t), &hash_form_t, &eq_form_t, &del_form_t);
    htab_init(&Q, table_size, sizeof(form_t), &hash_form_t, &eq_form_t, &del_form_t);

    int result[MAX_RANK];
    
    printf("Calling compute_group_bjt...\n");
    int rank = compute_group_bjt(result, D, init_pow, h_star, ell, &R, &Q);
    
    printf("Computation finished. Rank: %d\n", rank);
    printf("Result: %d\n", result[0]);

    htab_clear(&R);
    htab_clear(&Q);

    return 0;
}
