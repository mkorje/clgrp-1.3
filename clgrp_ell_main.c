#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#ifdef WITH_PARI
#include <pari/pari.h>
#endif

#include "clgrp_ell.h"

static const int congruences[4][2] = {{3, 8}, {7, 8}, {4, 16}, {8, 16}};
#define NUM_CONGRUENCES 4

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
    const long total_work = NUM_CONGRUENCES * files;

    int myrank, idx = 0;
    int work_item[3]; /* {file_index, a, m} */

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        /* Master process: verify input files and distribute work */

        printf("clgrp_ell: D_max=%ld, files=%ld, ell=%ld, folder=%s\n",
               D_max, files, ell, folder);
        fflush(stdout);

        /* Verify all input files exist before starting */
        for (int ci = 0; ci < NUM_CONGRUENCES; ci++)
        {
            int a = congruences[ci][0], m = congruences[ci][1];
            if (!verify_input_files_exist(folder, a, m, files))
            {
                fprintf(stderr, "Error: Not all input files exist for a=%d, m=%d. Aborting.\n", a, m);
                MPI_Abort(MPI_COMM_WORLD, 1);
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
                int ci = j / files;
                work_item[0] = j % files;
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

        /* Wait for all active workers to finish and terminate them */
        while (active_workers > 0)
        {
            MPI_Recv(&idx, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int term[3] = {-1, -1, -1};
            MPI_Send(term, 3, MPI_INT, idx, 0, MPI_COMM_WORLD);
            active_workers--;
        }

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
