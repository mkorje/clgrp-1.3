//! Analyze ℓ-adic growth in class group tabulations
//!
//! Reads fundamental discriminant files (cl[a]mod[m].[index].gz) and
//! index-ℓ² files (cl[a]mod[m]l[ell].[index].gz) to find discriminants
//! where a ℤ/ℓ^N factor grows to ℤ/ℓ^(N+1).

use clap::Parser;
use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "ell_growth")]
#[command(about = "Count discriminants with ℓ-adic growth in class groups")]
struct Args {
    /// Base folder containing cl[a]mod[m]/ and cl[a]mod[m]l[ell]/ directories
    folder: PathBuf,

    /// Prime ℓ for growth analysis
    #[arg(short, long)]
    ell: u64,

    /// Maximum |discriminant|
    #[arg(short = 'D', long)]
    d_max: i64,

    /// Number of files
    #[arg(short, long)]
    files: i64,

    /// Verbose output (print each matching discriminant)
    #[arg(short, long, default_value_t = false)]
    verbose: bool,

    /// Growth detection mode:
    /// - "strict": fund has ℓ^N, loses one, and gains one ℓ^(N+1) in ell
    /// - "any": fund has ℓ^N and ell has ℓ^(N+1) (regardless of other factors)
    /// - "net": total ℓ^(N+1) count increases (fund -> ell)
    #[arg(long, default_value = "strict")]
    mode: String,
}

/// All congruence classes for fundamental discriminants
const CONGRUENCE_CLASSES: [(i32, i32); 4] = [
    (8, 16), // D ≡ 8 mod 16
    (4, 16), // D ≡ 4 mod 16 (equivalently, D ≡ 12 mod 16)
    (3, 8),  // D ≡ 3 mod 8
    (7, 8),  // D ≡ 7 mod 8
];

/// Compute ℓ-adic valuation: largest k such that ℓ^k | n
fn ell_valuation(mut n: u64, ell: u64) -> u32 {
    if n == 0 {
        return 0;
    }
    let mut k = 0;
    while n % ell == 0 {
        n /= ell;
        k += 1;
    }
    k
}

/// Extract the ℓ-adic profile from invariant factors
/// Returns a sorted Vec of valuations (descending), e.g., [2, 1] for ℤ/ℓ²ℤ × ℤ/ℓℤ
fn ell_profile(invariants: &[u64], ell: u64) -> Vec<u32> {
    let mut profile: Vec<u32> = invariants
        .iter()
        .map(|&c| ell_valuation(c, ell))
        .filter(|&v| v > 0)
        .collect();
    profile.sort_by(|a, b| b.cmp(a)); // Descending order
    profile
}

/// Count occurrences of exponent n in the profile
fn count_factor(profile: &[u32], n: u32) -> usize {
    profile.iter().filter(|&&v| v == n).count()
}

/// Parse a line from the fundamental file: "dist h c1 c2 ... ct"
fn parse_fundamental_line(line: &str) -> Option<(i64, u64, Vec<u64>)> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 2 {
        return None;
    }
    let dist: i64 = parts[0].parse().ok()?;
    let _h: u64 = parts[1].parse().ok()?;
    let invariants: Vec<u64> = parts[2..].iter().filter_map(|s| s.parse().ok()).collect();
    Some((dist, _h, invariants))
}

/// Parse a line from the ell file: "dist kron c1 c2 ... ct"
fn parse_ell_line(line: &str) -> Option<(i64, i8, Vec<u64>)> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 2 {
        return None;
    }
    let dist: i64 = parts[0].parse().ok()?;
    let kron: i8 = parts[1].parse().ok()?;
    let invariants: Vec<u64> = parts[2..].iter().filter_map(|s| s.parse().ok()).collect();
    Some((dist, kron, invariants))
}

/// Results from processing one file pair
/// Now tracks counts for ALL values of N
#[derive(Debug, Default, Clone)]
struct FileResults {
    /// Total discriminants processed
    total: u64,
    /// For each N: (with_factor, with_growth)
    /// Key is N, value is (count with ℓ^N factor, count with growth to ℓ^(N+1))
    by_n: HashMap<u32, (u64, u64)>,
    /// Breakdown by (N, Kronecker symbol): ((N, kron), (with_factor, with_growth))
    by_n_kron: HashMap<(u32, i8), (u64, u64)>,
}

impl FileResults {
    fn merge(&mut self, other: FileResults) {
        self.total += other.total;
        for (n, (factor, growth)) in other.by_n {
            let entry = self.by_n.entry(n).or_insert((0, 0));
            entry.0 += factor;
            entry.1 += growth;
        }
        for ((n, kron), (factor, growth)) in other.by_n_kron {
            let entry = self.by_n_kron.entry((n, kron)).or_insert((0, 0));
            entry.0 += factor;
            entry.1 += growth;
        }
    }
}

/// Process a single file pair
fn process_file_pair(
    folder: &PathBuf,
    a: i32,
    m: i32,
    ell: u64,
    index: i64,
    d_total: i64,
    verbose: bool,
    mode: &str,
) -> Result<FileResults, Box<dyn std::error::Error + Send + Sync>> {
    let fund_path = folder
        .join(format!("cl{}mod{}", a, m))
        .join(format!("cl{}mod{}.{}.gz", a, m, index));

    let ell_path = folder
        .join(format!("cl{}mod{}l{}", a, m, ell))
        .join(format!("cl{}mod{}l{}.{}.gz", a, m, ell, index));

    // Open both files
    let fund_file = File::open(&fund_path)?;
    let ell_file = File::open(&ell_path)?;

    let fund_reader = BufReader::new(GzDecoder::new(fund_file));
    let ell_reader = BufReader::new(GzDecoder::new(ell_file));

    let mut results = FileResults::default();

    // Starting discriminant for this file
    let mut d_fund: i64 = index * d_total * (m as i64) + (a as i64);

    // Process lines in parallel (but must be synchronized since discriminant tracking is sequential)
    for (fund_line, ell_line) in fund_reader.lines().zip(ell_reader.lines()) {
        let fund_line = fund_line?;
        let ell_line = ell_line?;

        let Some((dist_fund, _h, fund_invariants)) = parse_fundamental_line(&fund_line) else {
            continue;
        };
        let Some((dist_ell, kron, ell_invariants)) = parse_ell_line(&ell_line) else {
            continue;
        };

        // Sanity check: distances should match
        if dist_fund != dist_ell {
            eprintln!(
                "Warning: distance mismatch at D={}: fund={}, ell={}",
                d_fund, dist_fund, dist_ell
            );
        }

        // Update discriminant
        d_fund += dist_fund * (m as i64);
        results.total += 1;

        // Compute ℓ-profiles
        let fund_profile = ell_profile(&fund_invariants, ell);
        let ell_prof = ell_profile(&ell_invariants, ell);

        // Find the maximum N we need to check (max valuation in either profile)
        let max_n = fund_profile
            .iter()
            .chain(ell_prof.iter())
            .copied()
            .max()
            .unwrap_or(0);

        // Check each N from 1 to max_n
        for target_n in 1..=max_n {
            let fund_n_count = count_factor(&fund_profile, target_n);
            let fund_n1_count = count_factor(&fund_profile, target_n + 1);
            let ell_n_count = count_factor(&ell_prof, target_n);
            let ell_n1_count = count_factor(&ell_prof, target_n + 1);

            // Check if fundamental has factor of order ℓ^N
            if fund_n_count == 0 {
                continue;
            }

            // Record that this discriminant has an ℓ^N factor
            let n_entry = results.by_n.entry(target_n).or_insert((0, 0));
            n_entry.0 += 1;
            let nk_entry = results.by_n_kron.entry((target_n, kron)).or_insert((0, 0));
            nk_entry.0 += 1;

            // Detect growth based on mode
            let growth = match mode {
                "strict" => {
                    // Strict: one ℓ^N factor disappears AND one ℓ^(N+1) factor appears
                    ell_n1_count > fund_n1_count && ell_n_count < fund_n_count
                }
                "any" => {
                    // Any: fund has ℓ^N and ell has ℓ^(N+1)
                    fund_n_count > 0 && ell_n1_count > 0
                }
                "net" => {
                    // Net: total ℓ^(N+1) count increases
                    ell_n1_count > fund_n1_count
                }
                _ => {
                    // Default to strict
                    ell_n1_count > fund_n1_count && ell_n_count < fund_n_count
                }
            };

            if growth {
                n_entry.1 += 1;
                nk_entry.1 += 1;

                if verbose {
                    println!(
                        "D={}: N={}, kron={}, fund_profile={:?}, ell_profile={:?}",
                        d_fund, target_n, kron, fund_profile, ell_prof
                    );
                }
            }
        }
    }

    Ok(results)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    println!("ℓ-adic growth analysis");
    println!("======================");
    println!("folder: {:?}", args.folder);
    println!("ℓ={}", args.ell);
    println!("D_max={}, files={}", args.d_max, args.files);
    println!(
        "Target: all N where ℤ/{}^N ℤ → ℤ/{}^(N+1)ℤ growth",
        args.ell, args.ell
    );
    println!("Detection mode: {}", args.mode);
    println!();

    // Aggregate results across all congruence classes
    let mut grand_total = FileResults::default();
    let mut by_class: Vec<((i32, i32), FileResults)> = Vec::new();

    for (a, m) in CONGRUENCE_CLASSES {
        let d_total = args.d_max / (args.files * m as i64);

        println!("Processing {} mod {} ...", a, m);

        // Process all file pairs in parallel for this congruence class
        let mode = args.mode.clone();
        let results: Vec<_> = (0..args.files)
            .into_par_iter()
            .map(|index| {
                match process_file_pair(
                    &args.folder,
                    a,
                    m,
                    args.ell,
                    index,
                    d_total,
                    args.verbose,
                    &mode,
                ) {
                    Ok(r) => r,
                    Err(e) => {
                        eprintln!("Error processing file {} for {}mod{}: {}", index, a, m, e);
                        FileResults::default()
                    }
                }
            })
            .collect();

        // Merge results for this congruence class
        let mut class_results = FileResults::default();
        for r in results {
            class_results.merge(r);
        }

        grand_total.merge(class_results.clone());
        by_class.push(((a, m), class_results));
    }

    // Output per-class summaries
    println!();
    println!("Results by congruence class");
    println!("===========================");
    for ((a, m), results) in &by_class {
        println!();
        println!("{} mod {}:", a, m);
        println!("  Total discriminants: {}", results.total);

        // Get all N values and sort them
        let mut n_values: Vec<_> = results.by_n.keys().copied().collect();
        n_values.sort();

        for n in n_values {
            let (with_factor, with_growth) = results.by_n.get(&n).copied().unwrap_or((0, 0));
            println!(
                "  N={}: with ℓ^{} factor: {}, with growth to ℓ^{}: {} ({:.2}%)",
                n,
                n,
                with_factor,
                n + 1,
                with_growth,
                if with_factor > 0 {
                    100.0 * with_growth as f64 / with_factor as f64
                } else {
                    0.0
                }
            );

            // Kronecker breakdown for this N
            for kron in [-1i8, 0, 1] {
                if let Some(&(kf, kg)) = results.by_n_kron.get(&(n, kron)) {
                    let kron_name = match kron {
                        -1 => "inert",
                        0 => "ramified",
                        1 => "split",
                        _ => "unknown",
                    };
                    println!(
                        "      kron={:2} ({}): factor={}, growth={} ({:.2}%)",
                        kron,
                        kron_name,
                        kf,
                        kg,
                        if kf > 0 {
                            100.0 * kg as f64 / kf as f64
                        } else {
                            0.0
                        }
                    );
                }
            }
        }
    }

    // Output grand totals
    println!();
    println!("Grand Total (all congruence classes)");
    println!("====================================");
    println!("Total discriminants: {}", grand_total.total);
    println!();

    // Get all N values and sort them
    let mut n_values: Vec<_> = grand_total.by_n.keys().copied().collect();
    n_values.sort();

    println!("Summary table:");
    println!(
        "{:>4} {:>12} {:>12} {:>10}",
        "N", "with_factor", "with_growth", "rate"
    );
    println!("{}", "-".repeat(42));
    for n in &n_values {
        let (with_factor, with_growth) = grand_total.by_n.get(n).copied().unwrap_or((0, 0));
        println!(
            "{:>4} {:>12} {:>12} {:>9.4}%",
            n,
            with_factor,
            with_growth,
            if with_factor > 0 {
                100.0 * with_growth as f64 / with_factor as f64
            } else {
                0.0
            }
        );
    }

    println!();
    println!("Detailed breakdown by N and Kronecker symbol:");
    for n in n_values {
        let (with_factor, with_growth) = grand_total.by_n.get(&n).copied().unwrap_or((0, 0));
        println!();
        println!("N={}: ℤ/{}^{}ℤ → ℤ/{}^{}ℤ", n, args.ell, n, args.ell, n + 1);
        println!(
            "  Total: with_factor={}, with_growth={} ({:.4}%)",
            with_factor,
            with_growth,
            if with_factor > 0 {
                100.0 * with_growth as f64 / with_factor as f64
            } else {
                0.0
            }
        );

        for kron in [-1i8, 0, 1] {
            if let Some(&(kf, kg)) = grand_total.by_n_kron.get(&(n, kron)) {
                let kron_name = match kron {
                    -1 => "inert",
                    0 => "ramified",
                    1 => "split",
                    _ => "unknown",
                };
                println!(
                    "  kron={:2} ({}): factor={}, growth={} ({:.2}%)",
                    kron,
                    kron_name,
                    kf,
                    kg,
                    if kf > 0 {
                        100.0 * kg as f64 / kf as f64
                    } else {
                        0.0
                    }
                );
            }
        }
    }

    Ok(())
}
