/*=============================================================================

    This file is part of CLGRP.

    CLGRP is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    CLGRP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CLGRP; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/

/*
 * Count ell-rank statistics from clgrp_ell output files.
 *
 * C/MPI port of the ell-count Rust crate.
 *
 * Usage: mpirun -np [#procs] ./ell_count [folder] [ell] [d_max] [files]
 *
 *   folder  - base folder containing cl[a]mod[m]/ and cl[a]mod[m]l[ell]/
 *             subdirectories
 *   ell     - prime ell
 *   d_max   - maximum |discriminant|
 *   files   - number of files per congruence class
 *
 * Rank 0 is the coordinator; ranks 1..n-1 are workers.
 * Workers process one (congruence class, file index) pair at a time and send
 * back a flat array of counts.  Rank 0 accumulates cumulative statistics and
 * prints a CSV to stdout.
 */

#include <mpi.h>

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE 4096
#define MAX_EXP 20
#define NUM_CLASSES 4

/*
 * Flat counts array layout (51 uint64_t values):
 *
 *  [0]        total
 *  [1]        inert
 *  [2]        inert_cyclic
 *  [3]        inert_rank_inc
 *  [4]        inert_cyc_rank_inc
 *  [5..25]    inert_cyc_no_inc_by_exp[0..MAX_EXP]
 *  [26]       ramified
 *  [27]       ramified_cyclic
 *  [28]       ramified_rank_inc
 *  [29]       ramified_cyc_rank_inc
 *  [30..50]   ramified_cyc_no_inc_by_exp[0..MAX_EXP]
 */
#define COUNTS_N (5 + (MAX_EXP + 1) + 4 + (MAX_EXP + 1)) /* 51 */

#define C_TOTAL          0
#define C_INERT          1
#define C_INERT_CYC      2
#define C_INERT_RANK_INC 3
#define C_INERT_CYC_RI   4
#define C_INERT_BY_EXP   5                   /* [5..25] */
#define C_RAMIFIED       (5 + MAX_EXP + 1)  /* 26 */
#define C_RAM_CYC        (C_RAMIFIED + 1)   /* 27 */
#define C_RAM_RANK_INC   (C_RAMIFIED + 2)   /* 28 */
#define C_RAM_CYC_RI     (C_RAMIFIED + 3)   /* 29 */
#define C_RAM_BY_EXP     (C_RAMIFIED + 4)   /* [30..50] */

/* Response buffer sent from worker to coordinator:
 *   resp[0]       = work_idx (cast to uint64_t)
 *   resp[1..51]   = counts[COUNTS_N]
 */
#define RESP_N (1 + COUNTS_N) /* 52 */

/* The four congruence classes processed, matching the Rust crate. */
static const int CLASSES[NUM_CLASSES][2] = {
	{3, 8}, {7, 8}, {4, 16}, {8, 16}
};

/* ell-adic valuation: largest v such that ell^v divides n. */
static uint32_t ell_val(uint64_t n, uint64_t ell)
{
	if (n == 0)
		return 0;
	uint32_t v = 0;
	while (n % ell == 0) {
		n /= ell;
		v++;
	}
	return v;
}

/*
 * Process one (congruence class, file index) pair and write 51 aggregate
 * counts into `counts`.
 *
 * Reads two gzipped files in lockstep:
 *   fundamental: {folder}/cl{a}mod{m}/cl{a}mod{m}.{file_idx}.gz
 *                line format: dist h c1 c2 ... ct
 *   ell:         {folder}/cl{a}mod{m}l{ell}/cl{a}mod{m}l{ell}.{file_idx}.gz
 *                line format: dist kron c1 c2 ... ct
 */
static void process_file(
	const char *folder, int a, int m, long ell, long file_idx, long files,
	long d_max, uint64_t counts[COUNTS_N])
{
	memset(counts, 0, COUNTS_N * sizeof(uint64_t));

	long d_total   = d_max / (files * (long)m);
	long starting_d = file_idx * d_total * (long)m + (long)a;

	char fund_cmd[2048], ell_cmd[2048];
	snprintf(fund_cmd, sizeof(fund_cmd),
		"gunzip -c '%s/cl%dmod%d/cl%dmod%d.%ld.gz'",
		folder, a, m, a, m, file_idx);
	snprintf(ell_cmd, sizeof(ell_cmd),
		"gunzip -c '%s/cl%dmod%dl%ld/cl%dmod%dl%ld.%ld.gz'",
		folder, a, m, ell, a, m, ell, file_idx);

	FILE *fund_f = popen(fund_cmd, "r");
	if (fund_f == NULL) {
		fprintf(stderr,
			"Error opening fundamental file: a=%d m=%d idx=%ld\n",
			a, m, file_idx);
		return;
	}
	FILE *ell_f = popen(ell_cmd, "r");
	if (ell_f == NULL) {
		fprintf(stderr,
			"Error opening ell file: a=%d m=%d l=%ld idx=%ld\n",
			a, m, ell, file_idx);
		pclose(fund_f);
		return;
	}

	char fund_line[MAX_LINE], ell_line[MAX_LINE];
	long fund_d = starting_d;
	long ell_d  = starting_d;

	while (fgets(fund_line, MAX_LINE, fund_f) != NULL &&
	       fgets(ell_line,  MAX_LINE, ell_f)  != NULL)
	{
		/* --- Parse fundamental entry: dist h c1 c2 ... ct --- */
		char *tok = strtok(fund_line, " \t\n");
		if (tok == NULL)
			continue;
		long fund_dist = atol(tok);

		tok = strtok(NULL, " \t\n"); /* skip class number h */
		if (tok == NULL)
			continue;

		int      fund_ell_rank = 0;
		uint32_t fund_ell_exp  = 0;
		while ((tok = strtok(NULL, " \t\n")) != NULL) {
			uint64_t c = strtoull(tok, NULL, 10);
			if (c > 0) {
				uint32_t v = ell_val(c, (uint64_t)ell);
				if (v > 0)
					fund_ell_rank++;
				if (v > fund_ell_exp)
					fund_ell_exp = v;
			}
		}
		fund_d += fund_dist * (long)m;

		/* --- Parse ell entry: dist kron c1 c2 ... ct --- */
		tok = strtok(ell_line, " \t\n");
		if (tok == NULL)
			continue;
		long ell_dist = atol(tok);

		tok = strtok(NULL, " \t\n");
		if (tok == NULL)
			continue;
		int kron = atoi(tok);

		int ell_ell_rank = 0;
		while ((tok = strtok(NULL, " \t\n")) != NULL) {
			uint64_t c = strtoull(tok, NULL, 10);
			if (c > 0 && c % (uint64_t)ell == 0)
				ell_ell_rank++;
		}
		ell_d += ell_dist * (long)m;

		if (fund_d != ell_d) {
			fprintf(stderr,
				"Discriminant mismatch: fund=%ld ell=%ld"
				" (a=%d m=%d idx=%ld)\n",
				fund_d, ell_d, a, m, file_idx);
			break;
		}

		int cyclic   = (fund_ell_rank <= 1);
		int rank_inc = (ell_ell_rank > fund_ell_rank);

		counts[C_TOTAL]++;

		if (kron == -1) {
			counts[C_INERT]++;
			if (cyclic)            counts[C_INERT_CYC]++;
			if (rank_inc)          counts[C_INERT_RANK_INC]++;
			if (cyclic && rank_inc) counts[C_INERT_CYC_RI]++;
			if (cyclic && (!rank_inc || fund_ell_exp == 0)) {
				uint32_t e = (fund_ell_exp <= MAX_EXP) ? fund_ell_exp : MAX_EXP;
				counts[C_INERT_BY_EXP + e]++;
			}
		} else if (kron == 0) {
			counts[C_RAMIFIED]++;
			if (cyclic)            counts[C_RAM_CYC]++;
			if (rank_inc)          counts[C_RAM_RANK_INC]++;
			if (cyclic && rank_inc) counts[C_RAM_CYC_RI]++;
			if (cyclic && (!rank_inc || fund_ell_exp == 0)) {
				uint32_t e = (fund_ell_exp <= MAX_EXP) ? fund_ell_exp : MAX_EXP;
				counts[C_RAM_BY_EXP + e]++;
			}
		}
		/* kron == 1 (split): counted only in total */
	}

	pclose(fund_f);
	pclose(ell_f);
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	if (argc != 5) {
		fprintf(stderr,
			"Format: mpirun -np [#procs] ./ell_count"
			" [folder] [ell] [d_max] [files]\n");
		MPI_Finalize();
		exit(1);
	}

	const char *folder = argv[1];
	const long ell     = atol(argv[2]);
	const long d_max   = atol(argv[3]);
	const long files   = atol(argv[4]);
	const long step    = d_max / files;
	const int  N       = NUM_CLASSES * (int)files; /* total work items */

	int myrank, num_procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	if (myrank == 0)
	{
		/* Coordinator */
		int W       = num_procs - 1;            /* number of workers */
		int initial = (W < N) ? W : N;          /* initial batch size */

		uint64_t *by_index = (uint64_t *)calloc(files * COUNTS_N, sizeof(uint64_t));
		if (by_index == NULL) {
			fprintf(stderr, "Memory allocation failed.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		fprintf(stderr, "ell_count: l=%ld, D_max=%ld, files=%ld, step=%ld\n",
			ell, d_max, files, step);
		fprintf(stderr, "Processing %d file pairs ...\n", N);
		fflush(stderr);

		/* Send initial batch */
		for (int i = 0; i < initial; i++)
			MPI_Send(&i, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);

		uint64_t  resp[RESP_N];
		MPI_Status status;
		int done = 0;

		/* Main distribution loop: one item left to distribute per iteration */
		for (int i = initial; i < N; i++) {
			MPI_Recv(resp, RESP_N, MPI_UINT64_T, MPI_ANY_SOURCE,
				 0, MPI_COMM_WORLD, &status);
			int src      = status.MPI_SOURCE;
			int file_idx = (int)(resp[0]) % (int)files;

			uint64_t *slot = &by_index[(long)file_idx * COUNTS_N];
			for (int k = 0; k < COUNTS_N; k++)
				slot[k] += resp[1 + k];

			done++;
			fprintf(stderr, "\r  %d/%d", done, N);
			fflush(stderr);

			MPI_Send(&i, 1, MPI_INT, src, 0, MPI_COMM_WORLD);
		}

		/* Drain remaining outstanding responses */
		for (int i = 0; i < initial; i++) {
			MPI_Recv(resp, RESP_N, MPI_UINT64_T, MPI_ANY_SOURCE,
				 0, MPI_COMM_WORLD, &status);
			int file_idx = (int)(resp[0]) % (int)files;

			uint64_t *slot = &by_index[(long)file_idx * COUNTS_N];
			for (int k = 0; k < COUNTS_N; k++)
				slot[k] += resp[1 + k];

			done++;
			fprintf(stderr, "\r  %d/%d", done, N);
			fflush(stderr);
		}
		fprintf(stderr, "\n");
		fflush(stderr);

		/* Terminate all workers */
		int sentinel = -1;
		for (int i = 1; i <= W; i++)
			MPI_Send(&sentinel, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

		/* Compute cumulative sums across file indices */
		uint64_t *cumul = (uint64_t *)calloc(files * COUNTS_N, sizeof(uint64_t));
		if (cumul == NULL) {
			fprintf(stderr, "Memory allocation failed.\n");
			free(by_index);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		{
			uint64_t running[COUNTS_N];
			memset(running, 0, sizeof(running));
			for (long fi = 0; fi < files; fi++) {
				uint64_t *src_slot = &by_index[fi * COUNTS_N];
				for (int k = 0; k < COUNTS_N; k++)
					running[k] += src_slot[k];
				memcpy(&cumul[fi * COUNTS_N], running,
				       COUNTS_N * sizeof(uint64_t));
			}
		}

		/* Find max non-zero exponent columns to display */
		int max_inert_exp = 0, max_ram_exp = 0;
		for (long fi = 0; fi < files; fi++) {
			uint64_t *c = &cumul[fi * COUNTS_N];
			for (int e = MAX_EXP; e > max_inert_exp; e--) {
				if (c[C_INERT_BY_EXP + e] > 0) {
					max_inert_exp = e;
					break;
				}
			}
			for (int e = MAX_EXP; e > max_ram_exp; e--) {
				if (c[C_RAM_BY_EXP + e] > 0) {
					max_ram_exp = e;
					break;
				}
			}
		}

		/* CSV header */
		printf("boundary,total,inert,inert_cyclic,inert_rank_inc"
		       ",inert_cyc_rank_inc");
		for (int e = 0; e <= max_inert_exp; e++)
			printf(",inert_cyc_no_inc_e%d", e);
		printf(",ramified,ramified_cyclic,ramified_rank_inc"
		       ",ramified_cyc_rank_inc");
		for (int e = 0; e <= max_ram_exp; e++)
			printf(",ram_cyc_no_inc_e%d", e);
		printf(",total_cyc_rank_inc\n");

		/* CSV data rows — one row per file index */
		for (long fi = 0; fi < files; fi++) {
			uint64_t *c      = &cumul[fi * COUNTS_N];
			long      boundary = (fi + 1) * step;

			printf("%ld,%" PRIu64 ",%" PRIu64 ",%" PRIu64
			       ",%" PRIu64 ",%" PRIu64,
				boundary,
				c[C_TOTAL], c[C_INERT], c[C_INERT_CYC],
				c[C_INERT_RANK_INC], c[C_INERT_CYC_RI]);
			for (int e = 0; e <= max_inert_exp; e++)
				printf(",%" PRIu64, c[C_INERT_BY_EXP + e]);
			printf(",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64,
				c[C_RAMIFIED], c[C_RAM_CYC],
				c[C_RAM_RANK_INC], c[C_RAM_CYC_RI]);
			for (int e = 0; e <= max_ram_exp; e++)
				printf(",%" PRIu64, c[C_RAM_BY_EXP + e]);
			printf(",%" PRIu64 "\n",
				c[C_INERT_CYC_RI] + c[C_RAM_CYC_RI]);
		}
		fflush(stdout);

		free(by_index);
		free(cumul);
	}
	else
	{
		/* Worker */
		int        work_idx;
		uint64_t   resp[RESP_N];
		MPI_Status status;

		MPI_Recv(&work_idx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

		while (work_idx != -1)
		{
			int class_idx = work_idx / (int)files;
			int file_idx  = work_idx % (int)files;
			int a         = CLASSES[class_idx][0];
			int m         = CLASSES[class_idx][1];

			uint64_t counts[COUNTS_N];
			process_file(folder, a, m, ell, (long)file_idx,
				     files, d_max, counts);

			resp[0] = (uint64_t)work_idx;
			memcpy(&resp[1], counts, COUNTS_N * sizeof(uint64_t));

			MPI_Send(resp, RESP_N, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD);
			MPI_Recv(&work_idx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		}
	}

	MPI_Finalize();
	return 0;
}
