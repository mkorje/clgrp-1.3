mod local;

use clap::Parser;
use local::{EllStream, FundamentalStream, LocalSpec};
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Parser)]
#[command(name = "ell-count")]
#[command(about = "Count ℓ-rank statistics from clgrp_ell output files")]
struct Args {
    /// Base folder containing cl[a]mod[m]/ and cl[a]mod[m]l[ell]/ directories
    folder: PathBuf,

    /// Prime ℓ
    #[arg(short, long)]
    ell: u64,

    /// Maximum |discriminant|
    #[arg(short = 'D', long)]
    d_max: i64,

    /// Number of files per congruence class
    #[arg(short, long)]
    files: i64,
}

const CONGRUENCE_CLASSES: [(i32, i32); 4] = [(3, 8), (7, 8), (4, 16), (8, 16)];
const MAX_EXP: usize = 20;

#[derive(Clone)]
struct Counts {
    total: u64,
    inert: u64,
    inert_cyclic: u64,
    inert_rank_inc: u64,
    inert_cyc_rank_inc: u64,
    inert_cyc_no_inc_by_exp: [u64; MAX_EXP + 1],
    ramified: u64,
    ramified_cyclic: u64,
    ramified_rank_inc: u64,
    ramified_cyc_rank_inc: u64,
    ramified_cyc_no_inc_by_exp: [u64; MAX_EXP + 1],
}

impl Default for Counts {
    fn default() -> Self {
        Self {
            total: 0,
            inert: 0,
            inert_cyclic: 0,
            inert_rank_inc: 0,
            inert_cyc_rank_inc: 0,
            inert_cyc_no_inc_by_exp: [0; MAX_EXP + 1],
            ramified: 0,
            ramified_cyclic: 0,
            ramified_rank_inc: 0,
            ramified_cyc_rank_inc: 0,
            ramified_cyc_no_inc_by_exp: [0; MAX_EXP + 1],
        }
    }
}

impl Counts {
    fn add_assign(&mut self, other: &Counts) {
        self.total += other.total;
        self.inert += other.inert;
        self.inert_cyclic += other.inert_cyclic;
        self.inert_rank_inc += other.inert_rank_inc;
        self.inert_cyc_rank_inc += other.inert_cyc_rank_inc;
        self.ramified += other.ramified;
        self.ramified_cyclic += other.ramified_cyclic;
        self.ramified_rank_inc += other.ramified_rank_inc;
        self.ramified_cyc_rank_inc += other.ramified_cyc_rank_inc;
        for i in 0..=MAX_EXP {
            self.inert_cyc_no_inc_by_exp[i] += other.inert_cyc_no_inc_by_exp[i];
            self.ramified_cyc_no_inc_by_exp[i] += other.ramified_cyc_no_inc_by_exp[i];
        }
    }
}

/// Process a single file pair (one fundamental + one ell file), returning
/// aggregate counts for all discriminants in that file.
fn process_file(
    folder: &Path,
    a: i32,
    m: i32,
    ell: u64,
    file_idx: i64,
    files: i64,
    d_max: i64,
) -> Counts {
    let spec = LocalSpec {
        folder: folder.to_path_buf(),
        a,
        m,
        files,
        d_max,
    };

    let fund = match FundamentalStream::open_file(spec.clone(), ell, file_idx) {
        Ok(s) => s,
        Err(e) => {
            eprintln!(
                "Error opening {}: {}",
                spec.fundamental_path(file_idx).display(),
                e
            );
            return Counts::default();
        }
    };
    let ell_s = match EllStream::open_file(spec.clone(), ell, file_idx) {
        Ok(s) => s,
        Err(e) => {
            eprintln!(
                "Error opening {}: {}",
                spec.ell_path(ell, file_idx).display(),
                e
            );
            return Counts::default();
        }
    };

    let mut counts = Counts::default();

    for (fund_result, ell_result) in fund.zip(ell_s) {
        let Ok(fund_entry) = fund_result else {
            eprintln!(
                "Error reading fundamental file {}mod{} index {}",
                a, m, file_idx
            );
            break;
        };
        let Ok(ell_entry) = ell_result else {
            eprintln!("Error reading ell file {}mod{} index {}", a, m, file_idx);
            break;
        };

        assert_eq!(
            fund_entry.abs_discriminant, ell_entry.abs_discriminant,
            "discriminant mismatch at {}mod{} index {}",
            a, m, file_idx
        );

        let cyclic = fund_entry.ell_rank <= 1;
        let rank_inc = ell_entry.ell_rank > fund_entry.ell_rank;

        counts.total += 1;
        match ell_entry.kronecker {
            -1 => {
                counts.inert += 1;
                if cyclic {
                    counts.inert_cyclic += 1;
                }
                if rank_inc {
                    counts.inert_rank_inc += 1;
                }
                if cyclic && rank_inc {
                    counts.inert_cyc_rank_inc += 1;
                }
                if cyclic && (!rank_inc || fund_entry.ell_exponent == 0) {
                    let exp = fund_entry.ell_exponent as usize;
                    counts.inert_cyc_no_inc_by_exp[exp.min(MAX_EXP)] += 1;
                }
            }
            0 => {
                counts.ramified += 1;
                if cyclic {
                    counts.ramified_cyclic += 1;
                }
                if rank_inc {
                    counts.ramified_rank_inc += 1;
                }
                if cyclic && rank_inc {
                    counts.ramified_cyc_rank_inc += 1;
                }
                if cyclic && (!rank_inc || fund_entry.ell_exponent == 0) {
                    let exp = fund_entry.ell_exponent as usize;
                    counts.ramified_cyc_no_inc_by_exp[exp.min(MAX_EXP)] += 1;
                }
            }
            _ => {}
        }
    }

    counts
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let files = args.files;
    let step = args.d_max / files;

    eprintln!(
        "ell_count: ℓ={}, D_max={}, files={}, step={}",
        args.ell, args.d_max, files, step
    );

    // Build work items: (a, m, file_idx) for all classes × files
    let work: Vec<(i32, i32, i64)> = CONGRUENCE_CLASSES
        .iter()
        .flat_map(|&(a, m)| (0..files).map(move |idx| (a, m, idx)))
        .collect();

    let total_work = work.len();
    eprintln!("Processing {} file pairs in parallel ...", total_work);

    // Process all file pairs in parallel → Vec<(file_idx, Counts)>
    let done = AtomicUsize::new(0);
    let per_file: Vec<(i64, Counts)> = work
        .par_iter()
        .map(|&(a, m, idx)| {
            let counts = process_file(&args.folder, a, m, args.ell, idx, files, args.d_max);
            let n = done.fetch_add(1, Ordering::Relaxed) + 1;
            if n % 16 == 0 || n == total_work {
                eprint!("\r  {n}/{total_work} ({:.0}%)", 100.0 * n as f64 / total_work as f64);
            }
            (idx, counts)
        })
        .collect();
    eprintln!();

    eprintln!("Aggregating ...");

    // Sum across congruence classes for each file index
    let mut by_index = vec![Counts::default(); files as usize];
    for (idx, counts) in per_file {
        by_index[idx as usize].add_assign(&counts);
    }

    // Cumulative sums → checkpoint[i] = counts for all |D| < (i+1)*step
    let mut cumulative = Vec::with_capacity(files as usize);
    let mut running = Counts::default();
    for counts in &by_index {
        running.add_assign(counts);
        cumulative.push(running.clone());
    }

    // Determine max exponent columns to display
    let max_inert_exp = (0..=MAX_EXP)
        .rev()
        .find(|&e| cumulative.iter().any(|c| c.inert_cyc_no_inc_by_exp[e] > 0))
        .unwrap_or(0);
    let max_ram_exp = (0..=MAX_EXP)
        .rev()
        .find(|&e| {
            cumulative
                .iter()
                .any(|c| c.ramified_cyc_no_inc_by_exp[e] > 0)
        })
        .unwrap_or(0);

    // CSV header
    print!("boundary,total,inert,inert_cyclic,inert_rank_inc,inert_cyc_rank_inc");
    for e in 0..=max_inert_exp {
        print!(",inert_cyc_no_inc_e{}", e);
    }
    print!(",ramified,ramified_cyclic,ramified_rank_inc,ramified_cyc_rank_inc");
    for e in 0..=max_ram_exp {
        print!(",ram_cyc_no_inc_e{}", e);
    }
    println!(",total_cyc_rank_inc");

    // CSV data — one row per file index
    for (i, counts) in cumulative.iter().enumerate() {
        let boundary = (i as i64 + 1) * step;
        print!(
            "{},{},{},{},{},{}",
            boundary,
            counts.total,
            counts.inert,
            counts.inert_cyclic,
            counts.inert_rank_inc,
            counts.inert_cyc_rank_inc,
        );
        for e in 0..=max_inert_exp {
            print!(",{}", counts.inert_cyc_no_inc_by_exp[e]);
        }
        print!(
            ",{},{},{},{}",
            counts.ramified,
            counts.ramified_cyclic,
            counts.ramified_rank_inc,
            counts.ramified_cyc_rank_inc,
        );
        for e in 0..=max_ram_exp {
            print!(",{}", counts.ramified_cyc_no_inc_by_exp[e]);
        }
        println!(
            ",{}",
            counts.inert_cyc_rank_inc + counts.ramified_cyc_rank_inc
        );
    }

    Ok(())
}
