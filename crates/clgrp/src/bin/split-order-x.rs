use clgrp::bjt::compute_group_bjt;
use malachite::base::num::arithmetic::traits::KroneckerSymbol;
use malachite::integer::Integer;
use qform::S64Group;

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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: split-order-clgrp <ell> <n>");
        std::process::exit(1);
    }
    let ell: u64 = args[1].parse().unwrap_or_else(|_| {
        eprintln!("Invalid prime ell: {}", args[1]);
        std::process::exit(1);
    });
    let n: u64 = args[2].parse().unwrap_or_else(|_| {
        eprintln!("Invalid order n: {}", args[2]);
        std::process::exit(1);
    });

    // Compute bound = 4 * ell^n - 1 with overflow checks.
    let ell_to_n = (ell as u128).checked_pow(n as u32).unwrap_or_else(|| {
        eprintln!("Overflow computing ell^n");
        std::process::exit(1);
    });
    let bound = 4u128
        .checked_mul(ell_to_n)
        .and_then(|v| v.checked_sub(1))
        .unwrap_or_else(|| {
            eprintln!("Overflow computing 4*ell^n - 1");
            std::process::exit(1);
        });
    if bound > i64::MAX as u128 {
        eprintln!("Bound {} exceeds i64 range", bound);
        std::process::exit(1);
    }
    let bound = bound as i64;

    // Pre-compute primes for h_lower_bound. Need primes up to 2 * max(Q_TABLE).
    let sieve_limit = 2 * Q_TABLE[Q_TABLE.len() - 1] as usize + 100;
    let primes = sieve_primes(sieve_limit);

    let pdivs = prime_divisors(n);

    let mut abs_d: i64 = 3;
    while abs_d <= bound {
        let mut d = -abs_d;

        let kron = (&Integer::from(d)).kronecker_symbol(Integer::from(ell));
        if kron != 1 {
            abs_d = next_abs_d(abs_d);
            continue;
        }

        let Some(group) = S64Group::new(d) else {
            abs_d = next_abs_d(abs_d);
            continue;
        };

        let Some(f) = group.prime_form(ell as i32) else {
            abs_d = next_abs_d(abs_d);
            continue;
        };

        // Check f^n == identity.
        let f_n = group.pow_u32(&f, n as u32);
        if !group.is_id(&f_n) {
            abs_d = next_abs_d(abs_d);
            continue;
        }

        // Check that order is exactly n: for each prime divisor p of n, f^(n/p) != identity.
        let mut exact = true;
        for &p in &pdivs {
            let f_np = group.pow_u32(&f, (n / p) as u32);
            if group.is_id(&f_np) {
                exact = false;
                break;
            }
        }

        if !exact {
            abs_d = next_abs_d(abs_d);
            continue;
        }

        let mut depth = 0;
        let h_star = h_lower_bound(d, &primes);

        let result = compute_group_bjt(d, 1, h_star, 0).unwrap();
        println!("{}\n\t{}:\t{}\t{:?}", d, depth, result.h, result.invariants);

        d *= 64;
        depth = 3;

        while depth < 10 {
            let h_star = h_lower_bound(d, &primes);
            let result = compute_group_bjt(d, 1, h_star, 0).unwrap();
            println!("\t{}:\t{}\t{:?}", depth, result.h, result.invariants);

            d *= 4;
            depth += 1;
        }

        abs_d = next_abs_d(abs_d);
    }
}

/// Advance abs_d to the next value where -abs_d ≡ 0 or 1 (mod 4).
fn next_abs_d(abs_d: i64) -> i64 {
    if abs_d % 4 == 3 { abs_d + 1 } else { abs_d + 3 }
}
