use clgrp::bjt::compute_group_bjt;
use malachite::base::num::arithmetic::traits::KroneckerSymbol;
use malachite::integer::Integer;
use qform::{S64Form, S64Group};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::Instant;

fn prime_divisors(mut n: u64) -> Vec<u64> {
    let mut divs = Vec::new();
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            divs.push(d);
            while n % d == 0 {
                n /= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        divs.push(n);
    }
    divs
}

fn sieve_primes(limit: usize) -> Vec<u64> {
    if limit < 2 {
        return Vec::new();
    }
    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;
    let mut p = 2;
    while p * p <= limit {
        if is_prime[p] {
            let mut m = p * p;
            while m <= limit {
                is_prime[m] = false;
                m += p;
            }
        }
        p += 1;
    }
    is_prime
        .iter()
        .enumerate()
        .filter_map(|(i, &b)| if b { Some(i as u64) } else { None })
        .collect()
}

const Q_TABLE: [u64; 20] = [
    947, 2269, 3929, 6011, 8447, 11093, 14149, 17393, 20921, 24733, 28807, 33151, 37619, 42533,
    47507, 52859, 58321, 64231, 70099, 76463,
];

const C_TABLE: [f64; 5] = [
    9785883.57965490035712718963623046875,
    62938341.824559591710567474365234375,
    201442194.1823149621486663818359375,
    494557498.147336781024932861328125,
    1013054914.3328609466552734375,
];

fn fundamental_discriminant(d: i64) -> (i64, u64) {
    let mut n = d.unsigned_abs();
    let mut f = 1u64;

    let mut p = 3u64;
    while p * p <= n {
        while n % (p * p) == 0 {
            n /= p * p;
            f *= p;
        }
        p += 2;
    }

    while n % 4 == 0 {
        let q = n / 4;
        if q % 4 == 0 {
            n = q;
            f *= 2;
        } else if q % 4 == 3 {
            n = q;
            f *= 2;
            break;
        } else {
            break;
        }
    }

    (-(n as i64), f)
}

fn h_lower_bound(d: i64, primes: &[u64]) -> i32 {
    let d_abs = d.unsigned_abs();
    let (dk, cond) = fundamental_discriminant(d);

    let r = ((d_abs as f64).log10() as usize) / 5;
    let r_c = r.min(C_TABLE.len() - 1);
    let r_q = r.min(Q_TABLE.len() - 1);
    let c = C_TABLE[r_c];
    let q = Q_TABLE[r_q];
    let q2 = q * 2;

    let dk_int = Integer::from(dk);

    let mut e = 0.0_f64;
    let mut pi = 0_usize;
    let mut p = primes[pi];
    while p < q {
        let kron = (&dk_int).kronecker_symbol(Integer::from(p)) as f64;
        e += ((p as f64) / (p as f64 - kron)).ln();
        pi += 1;
        if pi >= primes.len() {
            break;
        }
        p = primes[pi];
    }

    let mut wt = 1.0_f64;
    for i in q..=p {
        wt -= (i as f64) * (i as f64).ln() / c;
    }

    while p < q2 {
        let kron = (&dk_int).kronecker_symbol(Integer::from(p)) as f64;
        e += wt * ((p as f64) / (p as f64 - kron)).ln();
        pi += 1;
        if pi >= primes.len() {
            break;
        }
        p = primes[pi];
        wt -= ((p - 1) as f64 * ((p - 1) as f64).ln() + p as f64 * (p as f64).ln()) / c;
    }

    let mut h =
        (e.exp() * (d_abs as f64).sqrt()) / (std::f64::consts::SQRT_2 * std::f64::consts::PI);

    if cond > 1 {
        for &p in &prime_divisors(cond) {
            let kron = (&dk_int).kronecker_symbol(Integer::from(p)) as f64;
            h *= 1.0 - kron / p as f64;
        }
    }

    let h_star = h as i32;
    if h_star < 5 {
        h_star + 1
    } else {
        h_star
    }
}

/// Compute the order of `f` in the class group of known order `h`.
fn element_order(group: &S64Group, f: &S64Form, h: u64) -> u64 {
    let mut ord = h;
    for p in prime_divisors(ord) {
        while ord % p == 0 {
            let exp = ord / p;
            if exp > u32::MAX as u64 {
                break;
            }
            if group.is_id(&group.pow_u32(f, exp as u32)) {
                ord /= p;
            } else {
                break;
            }
        }
    }
    ord
}

struct DiscResult {
    d: i64,
    order: u64,
    has_order_4: bool,
}

fn process_disc(d: i64, primes: &[u64]) -> Option<DiscResult> {
    let group = S64Group::new(d)?;
    let f = group.prime_form(2)?;

    let h_star = h_lower_bound(d, primes);
    let result = compute_group_bjt(d, 1, h_star, 0).ok()?;

    let has_order_4 = result.invariants.iter().any(|&inv| inv % 4 == 0);
    let order = element_order(&group, &f, result.h as u64);

    Some(DiscResult {
        d,
        order,
        has_order_4,
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: split2-order4 <n_max>");
        std::process::exit(1);
    }
    let n_max: u64 = args[1].parse().unwrap_or_else(|_| {
        eprintln!("Invalid n_max: {}", args[1]);
        std::process::exit(1);
    });

    let max_bound = 4u128
        .checked_mul(1u128.checked_shl(n_max as u32).unwrap_or_else(|| {
            eprintln!("Overflow: 2^{} too large", n_max);
            std::process::exit(1);
        }))
        .and_then(|v| v.checked_sub(1))
        .unwrap_or_else(|| {
            eprintln!("Overflow computing bound for n_max={}", n_max);
            std::process::exit(1);
        });
    if max_bound > i64::MAX as u128 {
        eprintln!("Bound exceeds i64 range");
        std::process::exit(1);
    }
    let max_bound = max_bound as u64;

    let sieve_limit = 2 * Q_TABLE[Q_TABLE.len() - 1] as usize + 100;
    let primes = sieve_primes(sieve_limit);

    // For ell=2, only D ≡ 1 (mod 8) have (D/2) = 1, i.e. |D| ≡ 7 (mod 8).
    let discs: Vec<i64> = (7..=max_bound)
        .step_by(8)
        .map(|abs_d| -(abs_d as i64))
        .collect();

    let total = discs.len();
    eprintln!(
        "Processing {} discriminants (n_max={}, |D| up to {})",
        total, n_max, max_bound
    );
    let started = Instant::now();
    let processed = AtomicUsize::new(0);
    let errors = AtomicUsize::new(0);

    let workers = thread::available_parallelism().map_or(1, |n| n.get());
    let chunk_size = discs.len().div_ceil(workers.max(1));

    let results: Vec<DiscResult> = thread::scope(|scope| {
        let handles: Vec<_> = discs
            .chunks(chunk_size)
            .map(|chunk| {
                let primes = &primes;
                let processed = &processed;
                let errors = &errors;
                scope.spawn(move || {
                    let mut local = Vec::new();
                    for &d in chunk {
                        match process_disc(d, primes) {
                            Some(r) => local.push(r),
                            None => {
                                errors.fetch_add(1, Ordering::Relaxed);
                            }
                        }
                        let done = processed.fetch_add(1, Ordering::Relaxed) + 1;
                        if done == total || done % 100_000 == 0 {
                            let secs = started.elapsed().as_secs_f64().max(0.001);
                            let rate = done as f64 / secs;
                            eprintln!(
                                "[progress] {}/{} ({:.1}%), {:.0} disc/s",
                                done,
                                total,
                                done as f64 * 100.0 / total as f64,
                                rate,
                            );
                        }
                    }
                    local
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    });

    let err_count = errors.load(Ordering::Relaxed);
    if err_count > 0 {
        eprintln!("{} discriminants failed (skipped)", err_count);
    }
    eprintln!(
        "Done in {:.1}s, {} results",
        started.elapsed().as_secs_f64(),
        results.len()
    );

    for n in 1..=n_max {
        let bound = 4u64.saturating_mul(1u64.checked_shl(n as u32).unwrap_or(u64::MAX));
        let bound = bound.saturating_sub(1);
        let mut count = 0u64;
        let mut any_has_order_4 = false;
        for r in &results {
            if r.order == n && r.d.unsigned_abs() <= bound {
                count += 1;
                if r.has_order_4 {
                    any_has_order_4 = true;
                }
            }
        }
        let all_no_order_4 = !any_has_order_4;
        println!("{}\t{}\t({})", n, all_no_order_4, count);
    }
}
