#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include "clgrp.h"
#include "functions.h"
#include "sieve.h"

#define MAX_LINE_LENGTH 1024
#define MAX_INVARIANTS 20

// Helper function to extract the directory path (includes trailing '/')
char *get_directory(const char *filepath) {
  const char *last_slash = strrchr(filepath, '/');

  if (last_slash == NULL) {
    // No slashes found; file is in the current directory
    // Return "./" or an empty string depending on your preference
    return strdup("./");
  }

  // Calculate length up to and including the slash
  // Pointer arithmetic: (end_pointer - start_pointer) + 1 for length
  size_t length = last_slash - filepath + 1;

  // Allocate memory (len + 1 for the null terminator)
  char *dir_path = (char *)malloc(length + 1);
  if (dir_path) {
    strncpy(dir_path, filepath, length);
    dir_path[length] = '\0'; // Always null-terminate manually after strncpy
  }

  return dir_path;
}

int process_file(long ell, int m, int a, int index, FILE *infd, FILE *outfd) {
  int *primes;
  int **h_factors;

  int i, myrank, idx = 0;

  long D_block = 268435456; // 2^28

  // long D_root = sqrt(D_max * ell * ell * ell * ell);

  // primes = (int *)malloc(((int)(1.25506 * D_root / log(D_root))) * sizeof(int));
  // prime_sieve(D_root, primes);
  // primes = (int *)realloc(primes, (2 + primes[0]) * sizeof(int));

  // int h_max = 0;
  // long temp = 1;

  // // Ramare's bound
  // h_max = (1 / M_PI) * sqrt(D_max) * (0.5 * log(D_max) + 2.5 - log(6)) + 1;
  // h_max *= ell * (ell + 1);

  // // compute maximal number of prime factors of a class number
  // for (i = 1; temp < h_max; i++) {
  //   temp *= primes[i];
  // }

  // const int h_max_factors = i;
  // h_factors = (int **)malloc(h_max * sizeof(int *));

  // for (i = 0; i < h_max; i++) {
  //   h_factors[i] = (int *)malloc(h_max_factors * sizeof(int));
  // }

  // regular_sieve(h_max, h_max, h_factors, primes, 0);

  // char input_cmd[512], output_name[512], output_dir[512];
  char line[MAX_LINE_LENGTH], data[256];

  // const long D_max = (index + 1) * D_total * m;
  // // compute an upper bound on the size of the table
  int h_max = h_upper_bound(-D_block * (index + 1) * ell * ell * ell * ell);
  int table_size = next_prime((((int)sqrt(h_max)) << 1) - 1);
  // fprintf(stderr, "D_max: %ld, h_max: %ld, max prime: %ld\n", D_max, h_max,
  // (((long) sqrt(h_max)) << 1) - 1);
  if (table_size == -1) {
    perror("Not enough primes in liboptarith/primes.h\n");
    fflush(stderr);
    exit(1);
  }

  htab_t R, Q;
  htab_init(&R, table_size, sizeof(form_t), &hash_form_t, &eq_form_t,
            &del_form_t);
  htab_init(&Q, table_size, sizeof(form_t), &hash_form_t, &eq_form_t,
            &del_form_t);

  /* Calculate starting discriminant */
  long D = (long)index * D_block + a;
  int dist, h;
  int input_invariants[MAX_INVARIANTS];
  int result[MAX_INVARIANTS];
  int input_rank, output_rank;
  char kron;
  char output_line[MAX_LINE_LENGTH];
  long D_sub;

  int init_pow = 1, h_fac, h_temp, h_fac_total;
  const int *h_cur_factors;

  /* Process each line */
  while (fgets(line, MAX_LINE_LENGTH, infd) != NULL) {
    /* Parse input line: dist h c1 c2 ... ct */
    char *token = strtok(line, " \t\n");
    if (token == NULL)
      continue;

    dist = atoi(token);

    token = strtok(NULL, " \t\n");
    if (token == NULL)
      continue;
    h = atoi(token);

    // /* Parse invariant factors */
    // input_rank = 0;
    // while ((token = strtok(NULL, " \t\n")) != NULL && input_rank <
    // MAX_INVARIANTS)
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

    if (D == 4) {
      h /= 2;
    } else if (D == 3) {
      h /= 3;
    }

    /* Compute class structure of order of index ell^2 */
    init_pow = 1;
    h_cur_factors = h_factors[h];
    h_fac_total = h_cur_factors[0];
    h_temp = h;

    for (int i = 1; i <= h_fac_total; i++) {
      h_fac = h_cur_factors[i];
      h_temp /= h_fac;
      if (h_temp % h_fac != 0) {
        init_pow *= h_fac;
      }
    }

    h /= init_pow;
    output_rank = compute_group_bjt(result, -D_sub, init_pow, h, ell, &R, &Q);

    h = result[0] * init_pow;
    result[1] *= init_pow;

    /* Format output line: dist kron c1 c2 ... ct */
    sprintf(output_line, "%d\t%d\t", dist, (int)kron);
    for (int r = 1; r < output_rank; r++) {
      sprintf(data, "%d ", result[r]);
      strcat(output_line, data);
    }
    sprintf(data, "%d\n", result[output_rank]);
    strcat(output_line, data);

    fputs(output_line, outfd);
    fflush(outfd);
    break;
  }

  htab_clear(&R);
  htab_clear(&Q);

  return 0;
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Format: ./clgrp_ell_new [file] [ell]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  file   - input file from LMFDB\n");
    fprintf(stderr,
            "  ell    - prime for Kronecker symbol and order computation\n");
    exit(1);
  }

  const char *input_path = argv[1];
  const long ell = atol(argv[2]);

  int a, m, index;

  /* Parse the filename */
  const char *filename = strrchr(input_path, '/');
  filename = (filename != NULL) ? filename + 1 : input_path;
  if (sscanf(filename, "cl%dmod%d.%d.gz", &a, &m, &index) != 3) {
    fprintf(stderr, "Invalid input path: %s\n", input_path);
    exit(1);
  }

    /* Get directory */
    char *dir = get_directory(input_path);
    if (dir == NULL) {
      fprintf(stderr, "Memory allocation failed.\n");
      exit(1);
    }

    /* Open input file via gunzip */
    char input_cmd[512];
    sprintf(input_cmd, "gunzip -c %s", input_path);
    FILE *infd = popen(input_cmd, "r");
    if (infd == NULL) {
      fprintf(stderr, "Unable to open input file: %s\n", input_path);
      free(dir);
      exit(1);
    }

    /* Open output file */
    char output_path[512];
    snprintf(output_path, sizeof(output_path), "%scl%dmod%dell%ld.%d", dir, a,
             m, ell, index);
    FILE *outfd = fopen(output_path, "w");
    if (outfd == NULL) {
      fprintf(stderr, "Unable to create output file: %s\n", output_path);
      free(dir);
      pclose(infd);
      exit(1);
    }

    /* Print information before starting */
    printf("clgrp_ell_new: ℓ=%ld,    |d|≡%d (mod %d),    %d⋅2²⁸≤|d|<%d⋅2²⁸\n",
           ell, a, m, index, index + 1);
    fflush(stdout);

    /* Timing setup */
    struct timeval begin, end;
    unsigned long exec_time;
    gettimeofday(&begin, NULL);
    
    /* Run */
    process_file(ell, m, a, index, infd, outfd);

    /* Timing results */
    gettimeofday(&end, NULL);
    exec_time =
        (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
    printf("took %.3f\n", exec_time / 1e6);
    fflush(stdout);

    /* Cleanup */
    free(dir);
    pclose(infd);
    fclose(outfd);

    /* Compress output file */
    sprintf(input_cmd, "gzip %s", output_path);
    system(input_cmd);

    return 0;
}
