#ifndef CLGRP_ELL_H
#define CLGRP_ELL_H

#include "functions.h"

/*
 * Verify that all input files exist for the given parameters.
 * Returns 1 if all files exist, 0 otherwise.
 */
int verify_input_files_exist(const char *folder, int a, int m, long files);

/*
 * Process a single clgrp input file and produce output with Kronecker symbols
 * and class structure of the order of index ell^2.
 *
 * Parameters:
 *   index   - file index (0 to files-1)
 *   D_total - discriminants per file divided by m
 *   folder  - base folder for input/output
 *   a       - congruence class (|D| = a mod m)
 *   m       - modulus
 *   ell     - prime for Kronecker symbol and order computation
 */
void process_clgrp_file(const int index, const long D_total,
                        const char *folder, const int a, const int m,
                        const long ell, int ** h_factors);

#endif /* CLGRP_ELL_H */
