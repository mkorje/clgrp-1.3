#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "clgrp_ell.h"
#include "sieve.h"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    if (argc != 7)
    {
        fprintf(stderr, "Format: mpirun -np [#procs] ./clgrp_ell [D_max] [files] [a] [m] [ell] [folder]\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "  D_max  - maximum |discriminant|\n");
        fprintf(stderr, "  files  - number of input files (must divide D_max)\n");
        fprintf(stderr, "  a      - congruence class (|D| = a mod m)\n");
        fprintf(stderr, "  m      - modulus (8 or 16)\n");
        fprintf(stderr, "  ell    - prime for Kronecker symbol and order computation\n");
        fprintf(stderr, "  folder - base folder containing cl[a]mod[m]/ directories\n");
        MPI_Finalize();
        exit(1);
    }

    long D_max = atol(argv[1]);
    const long files = atol(argv[2]);
    const int a = atoi(argv[3]);
    const int m = atoi(argv[4]);
    const long ell = atol(argv[5]);
    const char *folder = argv[6];
    const long D_total = D_max / (files * m);

    int * primes;
    int ** h_factors;

    int i, myrank, idx = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        /* Master process: verify input files and distribute work */

        printf("clgrp_ell: D_max=%ld, files=%ld, a=%d, m=%d, ell=%ld, folder=%s\n",
               D_max, files, a, m, ell, folder);
        printf("D_total=%ld\n", D_total);
        fflush(stdout);

        /* Verify all input files exist before starting */
        if (!verify_input_files_exist(folder, a, m, files))
        {
            fprintf(stderr, "Error: Not all input files exist. Aborting.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        printf("All %ld input files verified.\n", files);
        fflush(stdout);

        int num_procs;
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
        num_procs--;

        if (num_procs == 0)
        {
            fprintf(stderr, "Error: Need at least 2 MPI processes (1 master + 1 worker).\n");
            MPI_Finalize();
            exit(1);
        }

        /* Send initial work to all workers */
        for (i = 0; i < num_procs; i++)
        {
            MPI_Send(&i, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);
        }

        /* Distribute remaining work as workers complete */
        for (i = num_procs; i < files; i++)
        {
            MPI_Recv(&idx, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&i, 1, MPI_INT, idx, 0, MPI_COMM_WORLD);
        }

        /* Send termination signal to all workers */
        for (i = 0, idx = -1; i < num_procs; i++)
        {
            MPI_Send(&idx, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);
        }

        printf("All files processed.\n");
    }
    else
    {
        int h_max = 0;
        // Ramare's bound
        h_max = (1/M_PI) * sqrt(D_max) * (0.5 * log(D_max) + 2.5 - log(6)) + 1;
        h_max *= ell * (ell + 1);
        // fprintf(stderr, "DEBUG: h_max=%d\n", h_max);
        // fflush(stderr);

        long D_root = sqrt(D_max * ell * ell * ell * ell);
        long prime_bound = (D_root > h_max) ? D_root : h_max;

        primes = (int *) malloc(((int) (1.25506 * prime_bound / log(prime_bound))) * sizeof(int));
        prime_sieve(prime_bound, primes);
        primes = (int *) realloc(primes, (2 + primes[0]) * sizeof(int));

        long temp = 1;

        // compute maximal number of prime factors of a class number
        for (i = 1; temp < h_max; i++)
        {
            temp *= primes[i];
        }

        const int h_max_factors = i;
        h_factors = (int **) malloc(h_max * sizeof(int *));

        for (i = 0; i < h_max; i++)
        {
            h_factors[i] = (int *) malloc(h_max_factors * sizeof(int));
        }

        regular_sieve(h_max, h_max, h_factors, primes, 0);

        /* Worker process: receive file indices and process them */

        MPI_Recv(&idx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        while (idx != -1)
        {
            process_clgrp_file(idx, D_total, folder, a, m, ell, h_factors);
            MPI_Send(&myrank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(&idx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    MPI_Finalize();

    return 0;
}
