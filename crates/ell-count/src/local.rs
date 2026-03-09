use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Entry from a fundamental discriminant file (`cl{a}mod{m}/cl{a}mod{m}.{index}.gz`).
///
/// Line format: `dist h c1 c2 ... ct`
/// Statistics are precomputed during parsing to avoid allocating invariant factor vectors.
#[derive(Debug, Clone)]
pub struct FundamentalEntry {
    pub abs_discriminant: i64,
    pub ell_rank: usize,
    pub ell_exponent: u32,
}

/// Entry from an ell output file (`cl{a}mod{m}l{ell}/cl{a}mod{m}l{ell}.{index}.gz`).
///
/// Line format: `dist kron c1 c2 ... ct` (for split: `dist 1 0`)
/// Statistics are precomputed during parsing to avoid allocating invariant factor vectors.
#[derive(Debug, Clone)]
pub struct EllEntry {
    pub abs_discriminant: i64,
    pub kronecker: i8,
    pub ell_rank: usize,
}

/// Describes a set of local class group files for one congruence class.
#[derive(Debug, Clone)]
pub struct LocalSpec {
    pub folder: PathBuf,
    pub a: i32,
    pub m: i32,
    pub files: i64,
    pub d_max: i64,
}

impl LocalSpec {
    pub fn d_total(&self) -> i64 {
        self.d_max / (self.files * self.m as i64)
    }

    pub fn fundamental_path(&self, index: i64) -> PathBuf {
        self.folder
            .join(format!("cl{}mod{}", self.a, self.m))
            .join(format!("cl{}mod{}.{}.gz", self.a, self.m, index))
    }

    pub fn ell_path(&self, ell: u64, index: i64) -> PathBuf {
        self.folder
            .join(format!("cl{}mod{}l{}", self.a, self.m, ell))
            .join(format!("cl{}mod{}l{}.{}.gz", self.a, self.m, ell, index))
    }

    pub fn starting_abs_d(&self, index: i64) -> i64 {
        index * self.d_total() * self.m as i64 + self.a as i64
    }
}

type GzReader = BufReader<MultiGzDecoder<File>>;

fn open_gz(path: &Path) -> io::Result<GzReader> {
    Ok(BufReader::new(MultiGzDecoder::new(File::open(path)?)))
}

/// ℓ-adic valuation of n.
fn ell_val(mut n: u64, ell: u64) -> u32 {
    if n == 0 {
        return 0;
    }
    let mut v = 0;
    while n.is_multiple_of(ell) {
        n /= ell;
        v += 1;
    }
    v
}

/// Iterator over fundamental discriminant entries, yielded in order of
/// increasing |D|.
///
/// Transparently opens successive gzipped files across the given index
/// range and tracks the running discriminant via deltas.
pub struct FundamentalStream {
    spec: LocalSpec,
    ell: u64,
    file_idx: i64,
    end_idx: i64,
    reader: GzReader,
    d: i64,
    buf: String,
    done: bool,
}

impl FundamentalStream {
    /// Open a stream over a single file at the given index.
    pub fn open_file(spec: LocalSpec, ell: u64, index: i64) -> io::Result<Self> {
        let start = index;
        let end = index + 1;
        let d = spec.starting_abs_d(start);
        let reader = open_gz(&spec.fundamental_path(start))?;
        Ok(Self {
            spec,
            ell,
            file_idx: start,
            end_idx: end,
            reader,
            d,
            buf: String::new(),
            done: false,
        })
    }
}

impl Iterator for FundamentalStream {
    type Item = io::Result<FundamentalEntry>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            self.buf.clear();
            match self.reader.read_line(&mut self.buf) {
                Ok(0) => {
                    self.file_idx += 1;
                    if self.file_idx >= self.end_idx {
                        self.done = true;
                        return None;
                    }
                    self.d = self.spec.starting_abs_d(self.file_idx);
                    match open_gz(&self.spec.fundamental_path(self.file_idx)) {
                        Ok(r) => self.reader = r,
                        Err(e) => {
                            self.done = true;
                            return Some(Err(e));
                        }
                    }
                }
                Ok(_) => {
                    let mut iter = self.buf.split_whitespace();
                    let Some(dist) = iter.next().and_then(|s| s.parse::<i64>().ok()) else {
                        continue;
                    };
                    // skip class number
                    if iter.next().is_none() {
                        continue;
                    }

                    let mut ell_rank = 0usize;
                    let mut ell_exponent = 0u32;
                    for s in iter {
                        let Ok(c) = s.parse::<u64>() else {
                            continue;
                        };
                        if c > 0 {
                            let v = ell_val(c, self.ell);
                            if v > 0 {
                                ell_rank += 1;
                            }
                            if v > ell_exponent {
                                ell_exponent = v;
                            }
                        }
                    }

                    self.d += dist * self.spec.m as i64;

                    return Some(Ok(FundamentalEntry {
                        abs_discriminant: self.d,
                        ell_rank,
                        ell_exponent,
                    }));
                }
                Err(e) => {
                    self.done = true;
                    return Some(Err(e));
                }
            }
        }
    }
}

/// Iterator over ell file entries, yielded in order of increasing |D|.
///
/// Transparently opens successive gzipped files across the given index
/// range and tracks the running discriminant via deltas.
pub struct EllStream {
    spec: LocalSpec,
    ell: u64,
    file_idx: i64,
    end_idx: i64,
    reader: GzReader,
    d: i64,
    buf: String,
    done: bool,
}

impl EllStream {
    /// Open a stream over a single file at the given index.
    pub fn open_file(spec: LocalSpec, ell: u64, index: i64) -> io::Result<Self> {
        let start = index;
        let end = index + 1;
        let d = spec.starting_abs_d(start);
        let reader = open_gz(&spec.ell_path(ell, start))?;
        Ok(Self {
            spec,
            ell,
            file_idx: start,
            end_idx: end,
            reader,
            d,
            buf: String::new(),
            done: false,
        })
    }
}

impl Iterator for EllStream {
    type Item = io::Result<EllEntry>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            self.buf.clear();
            match self.reader.read_line(&mut self.buf) {
                Ok(0) => {
                    self.file_idx += 1;
                    if self.file_idx >= self.end_idx {
                        self.done = true;
                        return None;
                    }
                    self.d = self.spec.starting_abs_d(self.file_idx);
                    match open_gz(&self.spec.ell_path(self.ell, self.file_idx)) {
                        Ok(r) => self.reader = r,
                        Err(e) => {
                            self.done = true;
                            return Some(Err(e));
                        }
                    }
                }
                Ok(_) => {
                    let mut iter = self.buf.split_whitespace();
                    let Some(dist) = iter.next().and_then(|s| s.parse::<i64>().ok()) else {
                        continue;
                    };
                    let Some(kronecker) = iter.next().and_then(|s| s.parse::<i8>().ok()) else {
                        continue;
                    };

                    let mut ell_rank = 0usize;
                    for s in iter {
                        let Ok(c) = s.parse::<u64>() else {
                            continue;
                        };
                        if c > 0 && c.is_multiple_of(self.ell) {
                            ell_rank += 1;
                        }
                    }

                    self.d += dist * self.spec.m as i64;

                    return Some(Ok(EllEntry {
                        abs_discriminant: self.d,
                        kronecker,
                        ell_rank,
                    }));
                }
                Err(e) => {
                    self.done = true;
                    return Some(Err(e));
                }
            }
        }
    }
}
