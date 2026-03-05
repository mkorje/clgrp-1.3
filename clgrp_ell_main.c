#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>

#ifdef WITH_PARI
#include <pari/pari.h>
#endif

#include "clgrp_ell.h"

static const int congruences[4][2] = {{3, 8}, {7, 8}, {4, 16}, {8, 16}};
#define NUM_CONGRUENCES 4

static const int special_indices[] = {0, 255, 511};
#define NUM_SPECIAL ((int)(sizeof(special_indices) / sizeof(special_indices[0])))

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    if (argc != 5)
    {
        fprintf(stderr, "Format: mpirun -np [#procs] ./clgrp_ell [D_max] [files] [ell] [folder]\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "  D_max  - maximum |discriminant|\n");
        fprintf(stderr, "  files  - number of input files (must divide D_max)\n");
        fprintf(stderr, "  ell    - prime for Kronecker symbol and order computation\n");
        fprintf(stderr, "  folder - base folder containing cl[a]mod[m]/ directories\n");
        MPI_Finalize();
        exit(1);
    }

    long D_max = atol(argv[1]);
    const long files = atol(argv[2]);
    const long ell = atol(argv[3]);
    const char *folder = argv[4];
    const long total_work = NUM_CONGRUENCES * NUM_SPECIAL;

    int myrank, idx = 0;
    int work_item[3]; /* {file_index, a, m} */

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        /* Master process: verify input files and distribute work */

        printf("clgrp_ell: D_max=%ld, files=%ld, ell=%ld, folder=%s\n",
               D_max, files, ell, folder);
        fflush(stdout);

        /* Verify input files exist for special indices before starting */
        for (int ci = 0; ci < NUM_CONGRUENCES; ci++)
        {
            int a = congruences[ci][0], m = congruences[ci][1];
            for (int si = 0; si < NUM_SPECIAL; si++)
            {
                char name[512];
                sprintf(name, "%s/cl%dmod%d/cl%dmod%d.%d.gz", folder, a, m, a, m, special_indices[si]);
                if (access(name, F_OK) == -1)
                {
                    fprintf(stderr, "Missing input file: %s\n", name);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
        }

        printf("All %ld input files verified.\n", total_work);
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

        int active_workers = 0;

        /* Send initial work to all workers */
        for (long j = 0; j < num_procs; j++)
        {
            if (j < total_work)
            {
                int ci = j / NUM_SPECIAL;
                work_item[0] = special_indices[j % NUM_SPECIAL];
                work_item[1] = congruences[ci][0];
                work_item[2] = congruences[ci][1];
                MPI_Send(work_item, 3, MPI_INT, j + 1, 0, MPI_COMM_WORLD);
                active_workers++;
            }
            else
            {
                int term[3] = {-1, -1, -1};
                MPI_Send(term, 3, MPI_INT, j + 1, 0, MPI_COMM_WORLD);
            }
        }

        /* Distribute remaining work as workers complete */
        for (long j = num_procs; j < total_work; j++)
        {
            MPI_Recv(&idx, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int ci = j / files;
            work_item[0] = j % files;
            work_item[1] = congruences[ci][0];
            work_item[2] = congruences[ci][1];
            MPI_Send(work_item, 3, MPI_INT, idx, 0, MPI_COMM_WORLD);
        }

        /* Wait for all active workers to finish, then terminate them all
         * at once so they reach MPI_Finalize at roughly the same time
         * (avoids PMIx collective timeout from staggered Finalize calls) */
        int *finished = (int *)malloc(active_workers * sizeof(int));
        for (int i = 0; i < active_workers; i++)
        {
            MPI_Recv(&finished[i], 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (int i = 0; i < active_workers; i++)
        {
            int term[3] = {-1, -1, -1};
            MPI_Send(term, 3, MPI_INT, finished[i], 0, MPI_COMM_WORLD);
        }
        free(finished);

        printf("All files processed.\n");
    }
    else
    {
        		#ifdef WITH_PARI
		pari_init(1000000, 0);
		#endif

        // Ramare's bound
        int h_max = (1/M_PI) * sqrt(D_max) * (0.5 * log(D_max) + 2.5 - log(6)) + 1;
        h_max *= ell * (ell + 1);

        // Smallest prime factor sieve
        int *spf = (int *) malloc(h_max * sizeof(int));
        for (int j = 0; j < h_max; j++) spf[j] = j;
        for (int j = 2; (long)j * j < h_max; j++)
        {
            if (spf[j] == j)
            {
                for (int k = j * j; k < h_max; k += j)
                {
                    if (spf[k] == k) spf[k] = j;
                }
            }
        }

        /* Worker process: receive work items and process them */

        MPI_Recv(work_item, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        while (work_item[0] != -1)
        {
            int file_idx = work_item[0];
            int a = work_item[1];
            int m = work_item[2];
            long D_total = D_max / (files * m);
            process_clgrp_file(file_idx, D_total, folder, a, m, ell, spf);
            MPI_Send(&myrank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(work_item, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    MPI_Finalize();

    return 0;
}
