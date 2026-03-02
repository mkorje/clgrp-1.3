use std::collections::{BTreeSet, HashSet};

use clgrp::bjt::compute_group_bjt;
use malachite::base::num::arithmetic::traits::KroneckerSymbol;
use malachite::integer::Integer;
use qform::S64Group;

/// Iterator over imaginary quadratic order discriminants d < -4 with
/// d ≡ 0, 1 (mod 4) where 2 is split, yielded in increasing order of |d|:
/// -7, -8, -11, -12, -15, -16, -19, -20, …
struct Discriminants {
    abs_d: i64,
}

impl Discriminants {
    fn new() -> Self {
        Discriminants { abs_d: 7 }
    }
}

impl Iterator for Discriminants {
    type Item = i64;

    fn next(&mut self) -> Option<i64> {
        loop {
            self.abs_d = if self.abs_d % 4 == 3 {
                self.abs_d + 1
            } else {
                self.abs_d + 3
            };

            if (&Integer::from(-self.abs_d)).kronecker_symbol(Integer::from(2)) == 1 {
                break;
            }
        }
        Some(-self.abs_d)
    }
}

/// A set of (d, r) pairs from the universe {(d, r) : d > 3, r >= d - 1},
/// represented compactly as diagonals plus individual points.
///
/// A "diagonal" with offset c means (d, d + c) holds for ALL d > 3.
/// Individual points record specific (d, r) values not yet covered by a diagonal.
struct PairCoverage {
    /// Offsets c where (d, d + c) is known for all d > 3.
    diagonals: BTreeSet<i64>,
    /// Individual (d, r) points not covered by any diagonal.
    points: HashSet<(u64, u64)>,
}

impl PairCoverage {
    fn new() -> Self {
        PairCoverage {
            diagonals: BTreeSet::new(),
            points: HashSet::new(),
        }
    }

    /// Record that (d, d + offset) holds for all d > 3.
    fn add_diagonal(&mut self, offset: i64) {
        if self.diagonals.insert(offset) {
            // Remove now-redundant individual points on this diagonal.
            self.points.retain(|&(d, r)| r as i64 - d as i64 != offset);
        }
    }

    /// Record that a specific (d, r) pair holds.
    fn add_point(&mut self, d: u64, r: u64) {
        let offset = r as i64 - d as i64;
        if !self.diagonals.contains(&offset) {
            self.points.insert((d, r));
        }
    }

    /// Check whether (d, r) is covered.
    fn contains(&self, d: u64, r: u64) -> bool {
        let offset = r as i64 - d as i64;
        self.diagonals.contains(&offset) || self.points.contains(&(d, r))
    }

    /// Check whether all offsets in -1..=max_offset have their diagonal.
    fn has_all_diagonals_up_to(&self, max_offset: i64) -> bool {
        (-1..=max_offset).all(|c| self.diagonals.contains(&c))
    }
}

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

/// Simple sieve of Eratosthenes returning all primes up to `limit`.
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

// Tables from ANTL/include/ANTL/quadratic/quadratic_order.hpp (OQvals[])
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

/// Given a (possibly non-fundamental) negative discriminant D, return
/// the fundamental discriminant D_K and the conductor f such that D = f² · D_K.
fn fundamental_discriminant(d: i64) -> (i64, u64) {
    let mut n = d.unsigned_abs();
    let mut f = 1u64;

    // Divide out all odd squared prime factors.
    let mut p = 3u64;
    while p * p <= n {
        while n % (p * p) == 0 {
            n /= p * p;
            f *= p;
        }
        p += 2;
    }

    // Divide out powers-of-4 from the even part.
    while n % 4 == 0 {
        let q = n / 4;
        // Keep the factor of 4 if removing it would leave a non-discriminant
        // (i.e. q ≢ 3 mod 4, meaning -q ≢ 1 mod 4).
        if q % 4 == 0 {
            n = q;
            f *= 2;
        } else if q % 4 == 3 {
            // -q ≡ 1 mod 4 and q is squarefree (odd squares already removed) → fundamental.
            n = q;
            f *= 2;
            break;
        } else {
            // q ≡ 1 or 2 mod 4 → need the factor of 4.
            break;
        }
    }

    (-(n as i64), f)
}

/// Port of `h_lower_bound` from clgrp.c.
///
/// Computes a conditional lower bound h* on the class number h(D),
/// satisfying h* <= h(D) <= 2*h*.
///
/// Algorithm from ANTL/src/L_function/L_function_long.cpp, approximateL1_impl(2).
/// The Kronecker symbol is evaluated at the fundamental discriminant D_K
/// (discriminant of the maximal order) so that (D_K/p) is the correct
/// primitive character, including at p = 2.
fn h_lower_bound(d: i64, primes: &[u64]) -> i32 {
    let d_abs = d.unsigned_abs();
    let (dk, f) = fundamental_discriminant(d);

    let r = ((d_abs as f64).log10() as usize) / 5;
    let r_c = r.min(C_TABLE.len() - 1);
    let r_q = r.min(Q_TABLE.len() - 1);
    let c = C_TABLE[r_c];
    let q = Q_TABLE[r_q];
    let q2 = q * 2;

    let dk_int = Integer::from(dk);

    // Compute partial product for p < Q.
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

    // Compute weighted partial products for Q <= p < 2Q.
    let mut wt = 1.0_f64;
    let i_start = q;
    for i in i_start..=p {
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

    // exp(E) ≈ L(1, χ_{D_K}). The naive estimate is L·√|D|/(√2·π), which
    // uses √|D| = f·√|D_K| and thus already absorbs the conductor f. For
    // non-maximal orders we must also include the conductor correction
    // Π_{p|f} (1 − (D_K/p)/p) that appears in the class number formula.
    let mut h =
        (e.exp() * (d_abs as f64).sqrt()) / (std::f64::consts::SQRT_2 * std::f64::consts::PI);

    if f > 1 {
        for &p in &prime_divisors(f) {
            let kron = (&dk_int).kronecker_symbol(Integer::from(p)) as f64;
            h *= 1.0 - kron / p as f64;
        }
    }

    let h_star = h as i32;
    if h_star < 5 { h_star + 1 } else { h_star }
}

fn bound(n: u64) -> i64 {
    // Compute bound = 4 * 2^n - 1 with overflow checks.
    let two_to_n = (2u128).checked_pow(n as u32).unwrap_or_else(|| {
        eprintln!("Overflow computing 2^n");
        std::process::exit(1);
    });
    let bound = 4u128
        .checked_mul(two_to_n)
        .and_then(|v| v.checked_sub(1))
        .unwrap_or_else(|| {
            eprintln!("Overflow computing 4*2^n - 1");
            std::process::exit(1);
        });
    if bound > i64::MAX as u128 {
        eprintln!("Bound {} exceeds i64 range", bound);
        std::process::exit(1);
    }
    bound as i64
}

fn element_order(d: i64, h: u64) -> u64 {
    let group = S64Group::new(d).unwrap();
    let f = group.prime_form(2i32).unwrap();

    let pdivs = prime_divisors(h);
    let mut n = h;
    for &p in &pdivs {
        let f_np = group.pow_u32(&f, (n / p) as u32);
        if group.is_id(&f_np) {
            n = n / p;
        }
    }
    n
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: split-order-clgrp <n>");
        std::process::exit(1);
    }
    let max_n: u64 = args[1].parse().unwrap_or_else(|_| {
        eprintln!("Invalid max order n: {}", args[2]);
        std::process::exit(1);
    });

    // Pre-compute primes for h_lower_bound. Need primes up to 2 * max(Q_TABLE).
    let sieve_limit = 2 * Q_TABLE[Q_TABLE.len() - 1] as usize + 100;
    let primes = sieve_primes(sieve_limit);

    let bound = bound(max_n);
    let discs = Discriminants::new()
        .take_while(|d| d.abs() <= bound)
        .map(|d| {
            let h_star = h_lower_bound(d, &primes);
            let result = compute_group_bjt(d, 1, h_star, 0).unwrap();
            (d, result.h, result.invariants)
        })
        .filter(|(_, _, invariants)| invariants.iter().any(|i| i % 4 == 0))
        .filter_map(|(d, h, invariants)| {
            let order = element_order(d, h as u64);
            (order <= max_n).then(|| (d, h, invariants, order))
        });
}
