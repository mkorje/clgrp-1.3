use optarith::xgcd_binary_i64;
use qform::{S64Form, S64Group};
use std::collections::HashMap;
use std::error::Error;
use std::fmt::{Display, Formatter};

pub const MAX_RANK: usize = 10;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BjtResult {
    pub h: i64,
    pub invariants: Vec<i64>,
}

#[derive(Debug)]
pub enum BjtError {
    InvalidInput(&'static str),
    PrimeFormExhausted { prime_index: i32 },
    MaxRankExceeded { rank: usize, max_rank: usize },
    TableIndexOutOfBounds { index: usize, len: usize },
    Overflow(&'static str),
}

impl Display for BjtError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
            Self::PrimeFormExhausted { prime_index } => {
                write!(
                    f,
                    "unable to find next prime form after index {prime_index}"
                )
            }
            Self::MaxRankExceeded { rank, max_rank } => {
                write!(f, "rank {rank} exceeds MAX_RANK={max_rank}")
            }
            Self::TableIndexOutOfBounds { index, len } => {
                write!(f, "table index {index} out of bounds for length {len}")
            }
            Self::Overflow(msg) => write!(f, "arithmetic overflow: {msg}"),
        }
    }
}

impl Error for BjtError {}

#[derive(Debug, Clone)]
struct FormEntry {
    form: S64Form,
    value: Vec<i32>,
}

#[derive(Debug, Default)]
struct FormTable {
    entries: Vec<FormEntry>,
    by_ab: HashMap<(i32, i32), Vec<usize>>,
}

impl FormTable {
    fn size(&self) -> usize {
        self.entries.len()
    }

    fn get(&self, index: usize) -> Option<&FormEntry> {
        self.entries.get(index)
    }

    fn find(&self, form: &S64Form) -> Option<&FormEntry> {
        let key = (form.a, form.b);
        let indices = self.by_ab.get(&key)?;
        indices.iter().find_map(|&idx| self.entries.get(idx))
    }

    fn insert(&mut self, entry: FormEntry) {
        let idx = self.entries.len();
        let key = (entry.form.a, entry.form.b);
        self.entries.push(entry);
        self.by_ab.entry(key).or_default().push(idx);
    }

    fn empty(&mut self) {
        self.entries.clear();
        self.by_ab.clear();
    }

    fn delete_from(&mut self, index: usize) -> Result<(), BjtError> {
        let len = self.entries.len();
        if index >= len {
            return Err(BjtError::TableIndexOutOfBounds { index, len });
        }

        if index + 1 == len {
            let popped = self.entries.pop().expect("len > 0");
            let key = (popped.form.a, popped.form.b);
            if let Some(indices) = self.by_ab.get_mut(&key) {
                let _ = indices.pop();
                if indices.is_empty() {
                    self.by_ab.remove(&key);
                }
            }
            return Ok(());
        }

        // Rare path in this algorithm, but supported for correctness.
        self.entries.remove(index);
        self.rebuild_index();
        Ok(())
    }

    fn rebuild_index(&mut self) {
        self.by_ab.clear();
        for (idx, entry) in self.entries.iter().enumerate() {
            self.by_ab
                .entry((entry.form.a, entry.form.b))
                .or_default()
                .push(idx);
        }
    }
}

fn qform_pow_s32(group: &S64Group, base: &S64Form, exp: i32) -> Result<S64Form, BjtError> {
    if exp < 0 {
        let pos = exp
            .checked_neg()
            .ok_or(BjtError::Overflow("negating exponent"))?;
        let mut inv = *base;
        group.inverse(&mut inv);
        Ok(group.pow_u32(&inv, pos as u32))
    } else {
        Ok(group.pow_u32(base, exp as u32))
    }
}

fn next_form(
    group: &S64Group,
    init_pow: i32,
    mut prime_index: i32,
) -> Result<(i32, S64Form), BjtError> {
    let mut out = group.identity();
    while group.is_id(&out) {
        let search_index = prime_index
            .checked_add(1)
            .ok_or(BjtError::Overflow("prime_index + 1"))?;
        let Some((idx, prime_form)) = group.next_prime_form(search_index) else {
            return Err(BjtError::PrimeFormExhausted { prime_index });
        };
        prime_index = idx;
        out = qform_pow_s32(group, &prime_form, init_pow)?;
    }
    Ok((prime_index, out))
}

fn vec_set(v: &mut Vec<i32>, val: i32, j: usize) {
    if v.len() <= j {
        v.resize(j + 1, 0);
    }
    v[j] = val;
}

fn vec_add(dst: &mut [i64; MAX_RANK], v: &[i32], w: &[i32]) {
    let ws = w.len();
    let vs = v.len();

    if ws < vs {
        for k in 0..ws {
            dst[k] = i64::from(v[k]) + i64::from(w[k]);
        }
        for k in ws..vs {
            dst[k] = i64::from(v[k]);
        }
    } else {
        for k in 0..vs {
            dst[k] = i64::from(v[k]) + i64::from(w[k]);
        }
        for k in vs..ws {
            dst[k] = i64::from(w[k]);
        }
    }
}

fn as_i64(value: i128, context: &'static str) -> Result<i64, BjtError> {
    i64::try_from(value).map_err(|_| BjtError::Overflow(context))
}

pub fn smith_normal_form(matrix: &mut [[i64; MAX_RANK]], size: usize) -> Result<(), BjtError> {
    if size <= 1 {
        return Ok(());
    }

    let n = size;
    let mut i: isize = n as isize - 1;

    while i > 0 {
        let ii = i as usize;

        loop {
            let mut c = 0;

            for j in (0..ii).rev() {
                if matrix[ii][j] == 0 {
                    continue;
                }

                let (dd, u, v) = if matrix[ii][ii] != 0 && matrix[ii][j] % matrix[ii][ii] == 0 {
                    (matrix[ii][ii], 1_i64, 0_i64)
                } else {
                    let xg = xgcd_binary_i64(matrix[ii][ii], matrix[ii][j]);
                    (xg.gcd, xg.s, xg.t)
                };

                if dd == 0 {
                    return Err(BjtError::InvalidInput("zero gcd in smith_normal_form"));
                }

                let r = matrix[ii][ii] / dd;
                let b = matrix[ii][j] / dd;

                for row in matrix.iter_mut().take(n) {
                    let temp = as_i64(
                        i128::from(u) * i128::from(row[ii]) + i128::from(v) * i128::from(row[j]),
                        "SNF row transform temp",
                    )?;
                    let new_kj = as_i64(
                        i128::from(r) * i128::from(row[j]) - i128::from(b) * i128::from(row[ii]),
                        "SNF row transform entry",
                    )?;
                    row[j] = new_kj;
                    row[ii] = temp;
                }
            }

            for j in (0..ii).rev() {
                if matrix[j][ii] == 0 {
                    continue;
                }

                let (dd, u, v) = if matrix[ii][ii] != 0 && matrix[j][ii] % matrix[ii][ii] == 0 {
                    (matrix[ii][ii], 1_i64, 0_i64)
                } else {
                    let xg = xgcd_binary_i64(matrix[ii][ii], matrix[j][ii]);
                    (xg.gcd, xg.s, xg.t)
                };

                if dd == 0 {
                    return Err(BjtError::InvalidInput("zero gcd in smith_normal_form"));
                }

                let r = matrix[ii][ii] / dd;
                let b = matrix[j][ii] / dd;

                let (upper, lower) = matrix.split_at_mut(ii);
                let row_i = &mut lower[0];
                let row_j = &mut upper[j];

                for (i_entry, j_entry) in row_i.iter_mut().zip(row_j.iter_mut()).take(n) {
                    let temp = as_i64(
                        i128::from(u) * i128::from(*i_entry) + i128::from(v) * i128::from(*j_entry),
                        "SNF column transform temp",
                    )?;
                    let new_jk = as_i64(
                        i128::from(r) * i128::from(*j_entry) - i128::from(b) * i128::from(*i_entry),
                        "SNF column transform entry",
                    )?;
                    *j_entry = new_jk;
                    *i_entry = temp;
                }

                c += 1;
            }

            if c == 0 {
                break;
            }
        }

        let b = matrix[ii][ii];
        let mut done = true;

        if b != 0 {
            'scan: for k in 0..ii {
                for l in 0..ii {
                    if matrix[k][l] % b != 0 {
                        let (upper, lower) = matrix.split_at_mut(ii);
                        let row_k = &upper[k];
                        let row_i = &mut lower[0];
                        for (i_entry, k_entry) in row_i.iter_mut().zip(row_k.iter()).take(n) {
                            *i_entry = as_i64(
                                i128::from(*i_entry) + i128::from(*k_entry),
                                "SNF step 9 accumulation",
                            )?;
                        }
                        done = false;
                        break 'scan;
                    }
                }
            }
        }

        if done {
            i -= 1;
        }
    }

    Ok(())
}

// Port of `compute_group_bjt` from `clgrp.c`.
// Returns:
// - `h`: order of the subgroup generated by this stage (caller in C multiplies by `init_pow`)
// - `invariants`: SNF diagonal factors > 1 (or `[1]` for the trivial group)
pub fn compute_group_bjt(
    discriminant: i64,
    init_pow: i32,
    h_star: i32,
    _ell: i32,
) -> Result<BjtResult, BjtError> {
    if h_star < 1 {
        return Err(BjtError::InvalidInput("h_star must be >= 1"));
    }
    if init_pow < 1 {
        return Err(BjtError::InvalidInput("init_pow must be >= 1"));
    }

    let mut m = [[0_i64; MAX_RANK]; MAX_RANK];
    let group =
        S64Group::new(discriminant).ok_or(BjtError::InvalidInput("invalid discriminant"))?;

    let mut r_table = FormTable::default();
    let mut q_table = FormTable::default();

    let id_entry = FormEntry {
        form: group.identity(),
        value: vec![0],
    };
    r_table.insert(id_entry.clone());
    q_table.insert(id_entry);

    let omega = 2_i32;
    let mut h = 1_i64;
    let mut det = 1_i32;
    let mut j = 0_usize;
    let mut prime_index = -1_i32;

    while h < i64::from(h_star) {
        if j >= MAX_RANK {
            return Err(BjtError::MaxRankExceeded {
                rank: j,
                max_rank: MAX_RANK,
            });
        }

        let q_size = q_table.size();
        let r_prev_size = r_table.size();

        let (new_prime_index, g) = next_form(&group, init_pow, prime_index)?;
        prime_index = new_prime_index;

        let mut s = 1_i32;
        let mut y = omega;
        let mut u = omega;
        let mut b = group.pow_u32(&g, omega as u32);
        let mut c = b;

        if j > 0 {
            for i in 0..q_size {
                let it = q_table.get(i).ok_or(BjtError::TableIndexOutOfBounds {
                    index: i,
                    len: q_size,
                })?;
                let temp = group.compose(&it.form, &g);
                if let Some(e) = r_table.find(&temp) {
                    let n = e.value.len().min(it.value.len()).min(j);
                    let mut blocked = false;
                    for (k, row) in m.iter().enumerate().take(n) {
                        if i64::from(it.value[k]) + i64::from(e.value[k]) >= row[k] {
                            blocked = true;
                            break;
                        }
                    }
                    if !blocked {
                        m[j][j] = 1;
                        break;
                    }
                }
            }
        }

        while m[j][j] == 0 {
            let mut baby_found = false;
            for i in s..=u {
                let a = qform_pow_s32(&group, &g, -i)?;

                if a.a == 1 {
                    m[j][j] = i64::from(i);
                    break;
                }

                if s == 1 && i > 1 {
                    for t in 0..r_prev_size {
                        let it =
                            r_table
                                .get(t)
                                .cloned()
                                .ok_or(BjtError::TableIndexOutOfBounds {
                                    index: t,
                                    len: r_prev_size,
                                })?;
                        let temp = group.compose(&it.form, &a);

                        let found_unblocked = if let Some(e) = q_table.find(&temp) {
                            let n = e.value.len().min(it.value.len()).min(j);
                            let mut blocked = false;
                            for (k, row) in m.iter().enumerate().take(n) {
                                if i64::from(it.value[k]) + i64::from(e.value[k]) >= row[k] {
                                    blocked = true;
                                    break;
                                }
                            }

                            if !blocked {
                                vec_add(&mut m[j], &it.value, &e.value);
                                m[j][j] = i64::from(i);
                                true
                            } else {
                                false
                            }
                        } else {
                            false
                        };

                        if found_unblocked {
                            baby_found = true;
                            break;
                        }

                        let mut value = it.value.clone();
                        vec_set(&mut value, i, j);
                        r_table.insert(FormEntry { form: temp, value });
                    }
                } else {
                    for t in 0..r_prev_size {
                        let it =
                            r_table
                                .get(t)
                                .cloned()
                                .ok_or(BjtError::TableIndexOutOfBounds {
                                    index: t,
                                    len: r_prev_size,
                                })?;
                        let form = group.compose(&it.form, &a);
                        let mut value = it.value.clone();
                        vec_set(&mut value, i, j);
                        r_table.insert(FormEntry { form, value });
                    }
                }

                if baby_found {
                    break;
                }
            }

            while m[j][j] == 0 && i64::from(y) < i64::from(u) * i64::from(u) {
                let mut giant_found = false;
                for t in 0..q_size {
                    let it = q_table
                        .get(t)
                        .cloned()
                        .ok_or(BjtError::TableIndexOutOfBounds {
                            index: t,
                            len: q_size,
                        })?;
                    let temp = group.compose(&it.form, &b);

                    if let Some(e) = r_table.find(&temp) {
                        let n = it.value.len().min(e.value.len()).min(j);
                        let mut blocked = false;
                        for (k, row) in m.iter().enumerate().take(n) {
                            if i64::from(it.value[k]) + i64::from(e.value[k]) >= row[k] {
                                blocked = true;
                                break;
                            }
                        }

                        if !blocked {
                            vec_add(&mut m[j], &it.value, &e.value);
                            m[j][j] = m[j][j]
                                .checked_add(i64::from(y))
                                .ok_or(BjtError::Overflow("M[j][j] += y"))?;
                            if m[j][j] != 0 {
                                giant_found = true;
                                break;
                            }
                        }
                    }
                }

                if giant_found {
                    break;
                }

                y = y.checked_add(u).ok_or(BjtError::Overflow("y += u"))?;
                b = group.compose(&b, &c);
            }

            s = u.checked_add(1).ok_or(BjtError::Overflow("s = u + 1"))?;
            u = u.checked_mul(2).ok_or(BjtError::Overflow("u *= 2"))?;
            c = group.square(&c);
        }

        if m[j][j] > 1 {
            h = h
                .checked_mul(m[j][j])
                .ok_or(BjtError::Overflow("h *= M[j][j]"))?;

            if h < i64::from(h_star) {
                let q = (m[j][j] as f64).sqrt().ceil() as i32;
                det = det.checked_mul(q).ok_or(BjtError::Overflow("det *= q"))?;
                let r_size = r_table.size();
                let det_usize =
                    usize::try_from(det).map_err(|_| BjtError::Overflow("det to usize"))?;

                if r_size > det_usize {
                    for t in (det_usize..r_size).rev() {
                        r_table.delete_from(t)?;
                    }
                }

                let q_size = q_table.size();
                c = group.pow_u32(&g, q as u32);
                b = group.identity();

                let mut i = 1_i32;
                while i < q && i64::from(i) * i64::from(q) < m[j][j] {
                    b = group.compose(&b, &c);
                    for t in 0..q_size {
                        let it =
                            q_table
                                .get(t)
                                .cloned()
                                .ok_or(BjtError::TableIndexOutOfBounds {
                                    index: t,
                                    len: q_size,
                                })?;
                        let form = group.compose(&b, &it.form);
                        let mut value = it.value.clone();
                        let iq = i.checked_mul(q).ok_or(BjtError::Overflow("i * q"))?;
                        vec_set(&mut value, iq, j);
                        q_table.insert(FormEntry { form, value });
                    }
                    i = i.checked_add(1).ok_or(BjtError::Overflow("i += 1"))?;
                }
            }

            j += 1;
        } else {
            m[j][j] = 0;
        }
    }

    if h > 2 * i64::from(h_star) {
        return Err(BjtError::InvalidInput("h > 2 * h_star"));
    }

    let invariants = if h != 1 {
        let rank = j;
        smith_normal_form(&mut m, rank)?;
        let mut factors = Vec::new();
        for (i, row) in m.iter().enumerate().take(rank) {
            if row[i] > 1 {
                factors.push(row[i]);
            }
        }
        factors
    } else {
        vec![1]
    };

    r_table.empty();
    q_table.empty();

    Ok(BjtResult { h, invariants })
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::MultiGzDecoder;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::thread;
    use std::time::Instant;

    #[test]
    fn compute_group_bjt_trivial_bound() {
        let out = compute_group_bjt(-23, 1, 1, 0).expect("compute should succeed");
        assert_eq!(out.h, 1);
        assert_eq!(out.invariants, vec![1]);
    }

    #[test]
    fn compute_group_bjt_nontrivial_bound() {
        let out = compute_group_bjt(-23, 1, 3, 0).expect("compute should succeed");
        assert_eq!(out.h, 3);
        assert_eq!(out.invariants, vec![3]);
    }

    #[test]
    fn smith_normal_form_keeps_diagonal_matrix() {
        let mut m = [[0_i64; MAX_RANK]; MAX_RANK];
        m[0][0] = 2;
        m[1][1] = 6;
        m[2][2] = 12;
        smith_normal_form(&mut m, 3).expect("snf should succeed");

        // This port follows the C implementation ordering convention:
        // diagonal factors are produced from largest to smallest.
        assert_eq!(m[0][0], 12);
        assert_eq!(m[1][1], 6);
        assert_eq!(m[2][2], 2);
        assert_eq!(m[0][0] % m[1][1], 0);
        assert_eq!(m[1][1] % m[2][2], 0);
    }

    fn distinct_prime_divisors(mut n: i64) -> Vec<i64> {
        let mut factors = Vec::new();
        let mut p = 2_i64;
        while p * p <= n {
            if n % p == 0 {
                factors.push(p);
                while n % p == 0 {
                    n /= p;
                }
            }
            p += if p == 2 { 1 } else { 2 };
        }
        if n > 1 {
            factors.push(n);
        }
        factors
    }

    // Mirrors the C code in tabulate_bjt: product of primes with exponent exactly one.
    fn init_pow_from_h(h: i64) -> i32 {
        let factors = distinct_prime_divisors(h);
        let mut init_pow = 1_i64;
        let mut h_temp = h;

        for p in factors {
            h_temp /= p;
            if h_temp % p != 0 {
                init_pow *= p;
            }
        }

        i32::try_from(init_pow).expect("init_pow should fit in i32")
    }

    #[derive(Debug, Clone)]
    struct SampleRow {
        discriminant: i64,
        class_number: i64,
        invariants: Vec<i64>,
    }

    type SourceKey = (u32, u32, u32);
    type SourceCounts = HashMap<SourceKey, usize>;
    type SourceSamples = HashMap<SourceKey, Vec<SampleRow>>;

    fn fixture_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("testdata/lmfdb_bjt_local_samples.tsv.gz")
    }

    fn next_prng(state: &mut u64) -> u64 {
        let mut x = *state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *state = x;
        x
    }

    fn parse_fixture_samples_per_source(sample_per_source: usize) -> (SourceCounts, SourceSamples) {
        let path = fixture_path();
        let file = File::open(&path)
            .unwrap_or_else(|e| panic!("failed to open fixture {}: {e}", path.display()));
        let reader = BufReader::new(MultiGzDecoder::new(file));

        let mut counts: SourceCounts = HashMap::new();
        let mut samples: SourceSamples = HashMap::new();
        let mut rng_states: HashMap<SourceKey, u64> = HashMap::new();

        for line_result in reader.lines() {
            let line = line_result.expect("read fixture line");
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let cols: Vec<_> = line.split('\t').collect();
            assert_eq!(cols.len(), 9, "expected 9 tab-separated columns in {line}");

            let r = cols[1]
                .parse::<u32>()
                .unwrap_or_else(|_| panic!("invalid r in {line}"));
            let m = cols[2]
                .parse::<u32>()
                .unwrap_or_else(|_| panic!("invalid m in {line}"));
            let k = cols[3]
                .parse::<u32>()
                .unwrap_or_else(|_| panic!("invalid k in {line}"));
            let key = (r, m, k);

            let discriminant = cols[6]
                .parse::<i64>()
                .unwrap_or_else(|_| panic!("invalid discriminant in {line}"));
            let class_number = cols[7]
                .parse::<i64>()
                .unwrap_or_else(|_| panic!("invalid class number in {line}"));
            let invariants = cols[8]
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|x| {
                    x.parse::<i64>()
                        .unwrap_or_else(|_| panic!("invalid invariant in {line}"))
                })
                .collect::<Vec<_>>();

            let seen = counts.entry(key).or_insert(0);
            *seen += 1;

            let bucket = samples.entry(key).or_default();
            if bucket.len() < sample_per_source {
                bucket.push(SampleRow {
                    discriminant,
                    class_number,
                    invariants,
                });
            } else {
                let state = rng_states.entry(key).or_insert_with(|| {
                    ((u64::from(r) << 40) ^ (u64::from(m) << 20) ^ u64::from(k))
                        .wrapping_add(0x9E37_79B9_7F4A_7C15)
                });
                let j = (next_prng(state) % (*seen as u64)) as usize;
                if j < sample_per_source {
                    bucket[j] = SampleRow {
                        discriminant,
                        class_number,
                        invariants,
                    };
                }
            }
        }

        (counts, samples)
    }

    fn validate_sample_row(row: &SampleRow) {
        let expected_h = row.class_number;
        let init_pow = init_pow_from_h(expected_h);
        let h_star = i32::try_from(expected_h / i64::from(init_pow)).expect("h_star");

        let out = compute_group_bjt(row.discriminant, init_pow, h_star, 0)
            .unwrap_or_else(|e| panic!("compute_group_bjt failed for D={}: {e}", row.discriminant));

        let mut actual_invariants = out.invariants.clone();
        if let Some(first) = actual_invariants.first_mut() {
            *first *= i64::from(init_pow);
        }
        let actual_h = out.h * i64::from(init_pow);

        assert_eq!(
            actual_h, expected_h,
            "class number mismatch for D={}",
            row.discriminant
        );
        assert_eq!(
            actual_invariants, row.invariants,
            "invariants mismatch for D={}",
            row.discriminant
        );
    }

    #[test]
    fn compute_group_bjt_matches_local_fixture_samples() {
        let expected_sources = 32_usize;
        let expected_rows_per_source = 50_000_usize;
        let sample_per_source = expected_rows_per_source;

        let (counts, samples) = parse_fixture_samples_per_source(sample_per_source);
        assert_eq!(
            counts.len(),
            expected_sources,
            "unexpected source-file count in local fixture"
        );
        assert_eq!(
            samples.len(),
            expected_sources,
            "missing sample bucket for some source files"
        );

        for (key, count) in &counts {
            assert_eq!(
                *count, expected_rows_per_source,
                "fixture row count drift for source {:?}",
                key
            );
        }

        for rows in samples.values() {
            assert_eq!(
                rows.len(),
                sample_per_source,
                "sample bucket did not retain requested size"
            );
        }

        let mut all_rows = Vec::with_capacity(sample_per_source * expected_sources);
        for rows in samples.into_values() {
            all_rows.extend(rows);
        }

        let total_rows = all_rows.len();
        let progress_every = 100_000_usize;
        let processed = AtomicUsize::new(0);
        let started_at = Instant::now();

        let workers = thread::available_parallelism().map_or(1, |n| n.get());
        let chunk_size = all_rows.len().div_ceil(workers.max(1));

        thread::scope(|scope| {
            let mut handles = Vec::new();
            let processed_ref = &processed;
            for chunk in all_rows.chunks(chunk_size) {
                handles.push(scope.spawn(move || {
                    for row in chunk {
                        validate_sample_row(row);
                        let processed_now = processed_ref.fetch_add(1, Ordering::Relaxed) + 1;
                        if processed_now == total_rows || processed_now.is_multiple_of(progress_every)
                        {
                            let pct = (processed_now as f64 * 100.0) / (total_rows as f64);
                            let secs = started_at.elapsed().as_secs_f64().max(0.001);
                            let rate = (processed_now as f64) / secs;
                            eprintln!(
                                "[bjt-test] {processed_now}/{total_rows} ({pct:.2}%) rows, {rate:.0} rows/s"
                            );
                        }
                    }
                }));
            }
            for handle in handles {
                handle.join().expect("parallel worker panicked");
            }
        });
    }
}
