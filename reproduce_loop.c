#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <liboptarith/primes.h>
#include "clgrp.h"

int main() {
    long D = -503316492500;
    int h_star = 400;
    int init_pow = 1371;
    long ell = 5;

    // Setup Group
    s64_qform_group_t group;
    s64_qform_group_init(&group);
    s64_qform_group_set_discriminant_s64(&group, D);
    group.conductor_ell = ell;

    group_pow_t gp;
    group_pow_init(&gp, &group.desc.group);

    // Setup Generator g (using same logic as clgrp.c)
    // We know from logs: g = (220263, -21698, 571802)
    s64_qform_t g;
    g.a = 220263;
    g.b = -21698;
    g.c = 571802;
    
    // Setup g_inv
    s64_qform_t g_inv;
    s64_qform_set(&group, &g_inv, &g);
    s64_qform_inverse(&group, &g_inv); // Should be (220263, 21698, 571802)

    printf("g = (%d, %d, %ld)\n", g.a, g.b, g.c);
    printf("g_inv = (%d, %d, %ld)\n", g_inv.a, g_inv.b, g_inv.c);

    // Run Loop up to 400
    s64_qform_t a;
    s64_qform_set_id(&group, &a); // Start at Identity
    
    printf("Start: a = (%d, %d, %ld)\n", a.a, a.b, a.c);

    for (int i = 1; i <= 405; i++) {
        s64_qform_compose(&group, &a, &a, &g_inv);
        
        if (i == 399 || i == 400 || i == 401) {
             printf("Step %d: a = (%d, %d, %ld)\n", i, a.a, a.b, a.c);
             if (a.a == 1) printf("  -> IDENTITY FOUND by a.a==1\n");
             if (s64_qform_is_id(&group, &a)) printf("  -> IDENTITY FOUND by s64_qform_is_id\n");
        }
    }

    return 0;
}
