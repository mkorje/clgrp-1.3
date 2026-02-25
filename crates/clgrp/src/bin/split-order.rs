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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: split-order <ell> <n>");
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
    let bound = 4u128.checked_mul(ell_to_n).and_then(|v| v.checked_sub(1)).unwrap_or_else(|| {
        eprintln!("Overflow computing 4*ell^n - 1");
        std::process::exit(1);
    });
    // We need bound to fit in i64 magnitude (positive side) since discriminants are negative i64.
    if bound > i64::MAX as u128 {
        eprintln!("Bound {} exceeds i64 range", bound);
        std::process::exit(1);
    }
    let bound = bound as i64;

    let pdivs = prime_divisors(n);

    // Enumerate D from -3 down to -bound where D ≡ 0 or 1 (mod 4).
    // D ≡ 0 (mod 4): -4, -8, -12, -16, ...
    // D ≡ 1 (mod 4): -3, -7, -11, -15, ...
    // Interleaved in descending order: -3, -4, -7, -8, -11, -12, -15, -16, ...
    let mut abs_d: i64 = 3;
    while abs_d <= bound {
        // abs_d ≡ 3 (mod 4) gives D = -abs_d ≡ 1 (mod 4)
        // abs_d ≡ 0 (mod 4) gives D = -abs_d ≡ 0 (mod 4)
        let d = -abs_d;

        // Kronecker symbol (D/ell) must be 1 for ell to split.
        let kron = (&Integer::from(d)).kronecker_symbol(Integer::from(ell));
        if kron != 1 {
            // Advance to next valid abs_d.
            if abs_d % 4 == 3 {
                abs_d += 1; // 3 -> 4
            } else {
                abs_d += 3; // 4 -> 7
            }
            continue;
        }

        let Some(group) = S64Group::new(d) else {
            if abs_d % 4 == 3 {
                abs_d += 1;
            } else {
                abs_d += 3;
            }
            continue;
        };

        let Some(f) = group.prime_form(ell as i32) else {
            if abs_d % 4 == 3 {
                abs_d += 1;
            } else {
                abs_d += 3;
            }
            continue;
        };

        // Check f^n == identity.
        let f_n = group.pow_u32(&f, n as u32);
        if !group.is_id(&f_n) {
            if abs_d % 4 == 3 {
                abs_d += 1;
            } else {
                abs_d += 3;
            }
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

        if exact {
            println!("{}", d);
        }

        // Advance to next valid abs_d.
        if abs_d % 4 == 3 {
            abs_d += 1;
        } else {
            abs_d += 3;
        }
    }
}
