#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "clgrp.h"
#include "clgrp_ell.h"
#include "functions.h"

#define MAX_LINE_LENGTH 1024
#define MAX_INVARIANTS 20

int verify_input_files_exist(const char *folder, int a, int m, long files)
{
    char name[512];

    for (long i = 0; i < files; i++)
    {
        sprintf(name, "%s/cl%dmod%d/cl%dmod%d.%ld.gz", folder, a, m, a, m, i);
        if (access(name, F_OK) == -1)
        {
            fprintf(stderr, "Missing input file: %s\n", name);
            return 0;
        }
    }

    return 1;
}

void process_clgrp_file(const int index, const long D_total,
                        const char *folder, const int a, const int m,
                        const long ell, int ** h_factors)
{
    char input_cmd[512], output_name[512], output_dir[512];
    char line[MAX_LINE_LENGTH], data[256];
    FILE *infd, *outfd;

    struct timeval begin, end;
    unsigned long exec_time;

    /* Check if output file already exists */
    sprintf(output_name, "%s/cl%dmod%dl%ld/cl%dmod%dl%ld.%d.gz",
            folder, a, m, ell, a, m, ell, index);
    if (access(output_name, F_OK) != -1)
    {
        printf("Output file %s already exists, skipping.\n", output_name);
        return;
    }

    /* Open input file via gunzip */
    sprintf(input_cmd, "gunzip -c %s/cl%dmod%d/cl%dmod%d.%d.gz",
            folder, a, m, a, m, index);
    infd = popen(input_cmd, "r");
    if (infd == NULL)
    {
        fprintf(stderr, "Unable to open input file for index %d\n", index);
        return;
    }




    const long D_max = (index + 1) * D_total * m;
    // compute an upper bound on the size of the table
    int h_max = h_upper_bound(-D_max * ell * ell * ell * ell);
    int table_size = next_prime((((int) sqrt(h_max)) << 1) - 1);

    if (table_size == -1)
    {
        perror("Not enough primes in liboptarith/primes.h\n");
        fflush(stderr);
        exit(1);
    }

    htab_t R, Q;
    htab_init(&R, table_size, sizeof(form_t), &hash_form_t, &eq_form_t, &del_form_t);
    htab_init(&Q, table_size, sizeof(form_t), &hash_form_t, &eq_form_t, &del_form_t);
    sprintf(output_dir, "%s/cl%dmod%dl%ld", folder, a, m, ell);
    mkdir(output_dir, 0744);
    sprintf(output_name, "%s/cl%dmod%dl%ld/cl%dmod%dl%ld.%d",
            folder, a, m, ell, a, m, ell, index);
    outfd = fopen(output_name, "w");
    if (outfd == NULL)
    {
        fprintf(stderr, "Unable to create output file %s\n", output_name);
        pclose(infd);
        return;
    }

    /* Calculate starting discriminant */
    long D = (long)index * D_total * m + a;
    int dist, h;
    // int input_invariants[MAX_INVARIANTS];
    int result[MAX_INVARIANTS];
    int /*input_rank,*/ output_rank;
    char kron;
    char output_line[MAX_LINE_LENGTH];
    long D_sub;

    int init_pow = 1, h_fac, h_temp, h_fac_total;
    const int * h_cur_factors;

    gettimeofday(&begin, NULL);

    /* Process each line */
    while (fgets(line, MAX_LINE_LENGTH, infd) != NULL)
    {
        /* Parse input line: dist h c1 c2 ... ct */
        char *token = strtok(line, " \t\n");
        if (token == NULL) continue;

        dist = atoi(token);

        token = strtok(NULL, " \t\n");
        if (token == NULL) continue;
        h = atoi(token);

        // /* Parse invariant factors */
        // input_rank = 0;
        // while ((token = strtok(NULL, " \t\n")) != NULL && input_rank < MAX_INVARIANTS)
        // {
        //     input_invariants[input_rank++] = atoi(token);
        // }

        /* Update discriminant */
        D += (long)dist * m;

        /* Compute Kronecker symbol (D/ell) */
        kron = kronecker_symbol(-D, ell);

        D_sub = D;
        if (kron == 0) {
            // ramified
            h *= ell;
            D_sub *= ell * ell;
        } else if (kron == -1) {
            // inert
            h *= (ell + 1) * ell;
            D_sub *= ell * ell * ell * ell;
        } else if (kron == 1) {
            // split
            h *= (ell - 1) * ell;
            D_sub *= ell * ell * ell * ell;
        }
        /* Compute class structure of order of index ell^2 */
				init_pow = 1;
				h_cur_factors = h_factors[h];
				h_fac_total = h_cur_factors[0];
				h_temp = h;

        // fprintf(stderr, "DEBUG: h=%d, D=%ld, D_sub=%ld, kron=%d, h_fac_total=%d\n", h, D, D_sub, kron, h_fac_total);
        // fflush(stderr);
				for (int i = 1; i <= h_fac_total; i++)
				{
					h_fac = h_cur_factors[i];
            // fprintf(stderr, "\th_fac=%d\n", h_fac);
            // fflush(stderr);
					h_temp /= h_fac;
					if (h_temp % h_fac != 0)
					{
						init_pow *= h_fac;
					}
				}

				// h /= init_pow;
        // fprintf(stderr, "DEBUG: init_pow=%d, h=%d\n", init_pow, h);
        // fflush(stderr);
        output_rank = compute_group_bjt(result, -D_sub, init_pow, h, &R, &Q);
        // h = result[0] * init_pow;
        result[1] *= init_pow;
        


        /* Format output line: dist kron c1 c2 ... ct */
        sprintf(output_line, "%d\t%d\t", dist, (int)kron);
        for (int r = 1; r < output_rank; r++)
        {
            // fprintf(stderr, "\t%d\n", result[r]);
            // fflush(stderr);
            sprintf(data, "%d ", result[r]);
            strcat(output_line, data);
        }
            // fprintf(stderr, "\t%d\n", result[output_rank]);
            // fflush(stderr);
        sprintf(data, "%d\n", result[output_rank]);
        strcat(output_line, data);

        fputs(output_line, outfd);
    }

    pclose(infd);
    fclose(outfd);

    /* Compress output file */
    sprintf(input_cmd, "gzip %s", output_name);
    system(input_cmd);

    htab_clear(&R);
    htab_clear(&Q);
    gettimeofday(&end, NULL);
    exec_time = (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
    printf("index=%d, ell=%ld, took %.3f\n", index, ell, exec_time / 1e6);
    fflush(stdout);
}
