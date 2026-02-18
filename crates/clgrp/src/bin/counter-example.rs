use clgrp::bjt::compute_group_bjt;
use clgrp::{LmfdbClassGroupEntry, LmfdbFileSpec, stream_lmfdb_class_groups};
use malachite::base::num::arithmetic::traits::KroneckerSymbol;
use malachite::integer::Integer;
use std::fmt::Write;

fn remove_ell_fast_i32(mut x: i32, ell: i32) -> i32 {
    // Basic guard rails
    if x == 0 || ell == 0 || ell == 1 {
        return x;
    }
    if ell == -1 {
        return x.abs();
    }

    let ell = ell.abs(); // Valuation is defined for the prime's magnitude

    // 1. Successive Squaring Strategy
    let mut p = ell;
    while let Some(p_squared) = p.checked_mul(p) {
        if x % p_squared == 0 {
            p = p_squared;
        } else {
            break;
        }
    }

    // 2. Strip the large blocks
    while x % p == 0 {
        x /= p;
    }

    // 3. Clean up the remainder (in case p was a large power)
    while x % ell == 0 {
        x /= ell;
    }

    x
}

fn main() {
    let mut buffer = String::new();

    let k = 0;
    let r = 3;
    let m = 8;
    let spec = LmfdbFileSpec { k, r, m };
    let mut stream = stream_lmfdb_class_groups(spec).unwrap();

    let ell: i32 = 5;

    let mut count = 0;
    while let Some(Ok(LmfdbClassGroupEntry {
        discriminant,
        class_number,
        ..
    })) = stream.next()
    {
        if discriminant == -4 || discriminant == -3 {
            continue;
        }

        let kronecker = (&Integer::from(discriminant)).kronecker_symbol(Integer::from(ell));
        if kronecker != -1 {
            continue;
        }

        let mut h: i32 = class_number.try_into().unwrap();
        let mut d = discriminant;

        h *= (ell + 1) * ell;
        d *= (ell * ell * ell * ell) as i64;

        let init = remove_ell_fast_i32(h, ell);
        let h_star = h / init;
        if h_star == ell {
            continue;
        }

        writeln!(
            &mut buffer,
            "kronecker={:?},h={:?},d={:?},init={:?},h*={:?}",
            kronecker, h, d, init, h_star
        );

        let invariants = compute_group_bjt(d, init, h_star, ell).unwrap();
        let exponent = *invariants.invariants.iter().max().unwrap_or(&0);
        let mut depth = 2;
        writeln!(
            &mut buffer,
            "d={:?}, {:?} {:?}",
            depth, exponent, invariants.invariants
        );

        loop {
            h *= ell;
            d *= (ell * ell) as i64;
            depth = depth + 1;

            let init = remove_ell_fast_i32(h, ell);
            let invariants = compute_group_bjt(d, init, h / init, 0).unwrap();
            let new_exponent: i64 = *invariants.invariants.iter().max().unwrap_or(&0);
            writeln!(
                &mut buffer,
                " ={:?}, {:?}: {:?}",
                depth, new_exponent, invariants.invariants
            );
            if new_exponent > exponent {
                if depth >= 4 {
                    print!("{}", buffer);
                } else {
                    buffer.clear();
                }
                break;
            }
        }

        count += 1;

        if count >= 10000 {
            break;
        }
    }
}
