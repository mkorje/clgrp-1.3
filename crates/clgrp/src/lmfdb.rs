use flate2::read::MultiGzDecoder;
use reqwest::blocking::Client;
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::io::{BufRead, BufReader};
use std::sync::mpsc::{self, Receiver, SyncSender, TryRecvError};

/// Stream class-group rows from one LMFDB `.gz` data file.
///
/// The download and gzip decoding happen on a background thread. Consumers can:
/// - call `try_next()` to poll without blocking, or
/// - iterate with `next()` to block until each parsed row is available.
///
/// ```no_run
/// use clgrp::{stream_lmfdb_class_groups, LmfdbFileSpec};
///
/// let spec = LmfdbFileSpec { k: 1, r: 4, m: 16 };
/// let mut stream = stream_lmfdb_class_groups(spec)?;
///
/// // Non-blocking poll
/// if let Some(row) = stream.try_next() {
///     let row = row?;
///     println!("first row d={} h={}", row.discriminant, row.class_number);
/// }
///
/// // Blocking iteration
/// for row in stream {
///     let row = row?;
///     // process row...
///     if row.line_index > 10 {
///         break;
///     }
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub const LMFDB_QICG_ENDPOINT: &str =
    "https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups";
const LMFDB_BUCKET_SIZE: i64 = 1_i64 << 28;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LmfdbFileSpec {
    pub k: u32,
    pub r: u32,
    pub m: u32,
}

impl LmfdbFileSpec {
    pub fn filename_base(self) -> String {
        format!("cl{}mod{}", self.r, self.m)
    }

    pub fn initial_discriminant(self) -> i64 {
        -(i64::from(self.k) * LMFDB_BUCKET_SIZE) - i64::from(self.r)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LmfdbClassGroupEntry {
    pub line_index: u64,
    pub delta: i64,
    pub discriminant: i64,
    pub class_number: u64,
    pub invariant_factors: Vec<u64>,
}

#[derive(Debug)]
pub enum LmfdbStreamError {
    InvalidSpec(&'static str),
    Http(reqwest::Error),
    Io(std::io::Error),
    Parse {
        line_index: u64,
        raw_line: String,
        message: &'static str,
    },
    Arithmetic {
        line_index: u64,
        message: &'static str,
    },
}

impl Display for LmfdbStreamError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidSpec(message) => write!(f, "invalid LMFDB file spec: {message}"),
            Self::Http(err) => write!(f, "HTTP error: {err}"),
            Self::Io(err) => write!(f, "I/O error: {err}"),
            Self::Parse {
                line_index,
                raw_line,
                message,
            } => {
                write!(
                    f,
                    "parse error on line {line_index}: {message}; raw line: {raw_line:?}"
                )
            }
            Self::Arithmetic {
                line_index,
                message,
            } => {
                write!(f, "arithmetic overflow on line {line_index}: {message}")
            }
        }
    }
}

impl Error for LmfdbStreamError {}

impl From<reqwest::Error> for LmfdbStreamError {
    fn from(value: reqwest::Error) -> Self {
        Self::Http(value)
    }
}

impl From<std::io::Error> for LmfdbStreamError {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}

#[derive(Debug)]
pub struct LmfdbClassGroupStream {
    rx: Receiver<Result<LmfdbClassGroupEntry, LmfdbStreamError>>,
}

impl LmfdbClassGroupStream {
    pub fn try_next(&mut self) -> Option<Result<LmfdbClassGroupEntry, LmfdbStreamError>> {
        match self.rx.try_recv() {
            Ok(item) => Some(item),
            Err(TryRecvError::Empty) => None,
            Err(TryRecvError::Disconnected) => None,
        }
    }
}

impl Iterator for LmfdbClassGroupStream {
    type Item = Result<LmfdbClassGroupEntry, LmfdbStreamError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.rx.recv().ok()
    }
}

pub fn lmfdb_download_url(spec: LmfdbFileSpec) -> String {
    format!(
        "{LMFDB_QICG_ENDPOINT}?filenamebase={}&k={}&Fetch=fetch",
        spec.filename_base(),
        spec.k
    )
}

pub fn stream_lmfdb_class_groups(
    spec: LmfdbFileSpec,
) -> Result<LmfdbClassGroupStream, LmfdbStreamError> {
    if spec.m == 0 {
        return Err(LmfdbStreamError::InvalidSpec("m must be non-zero"));
    }

    let client = Client::builder()
        .user_agent("clgrp-lmfdb-stream/0.1.0")
        .http1_only()
        .build()?;
    let (tx, rx) = mpsc::sync_channel(512);

    std::thread::spawn(move || {
        worker(spec, client, tx);
    });

    Ok(LmfdbClassGroupStream { rx })
}

fn worker(
    spec: LmfdbFileSpec,
    client: Client,
    tx: SyncSender<Result<LmfdbClassGroupEntry, LmfdbStreamError>>,
) {
    let url = lmfdb_download_url(spec);
    let response = match client.get(url).send() {
        Ok(response) => match response.error_for_status() {
            Ok(response) => response,
            Err(err) => {
                let _ = tx.send(Err(LmfdbStreamError::Http(err)));
                return;
            }
        },
        Err(err) => {
            let _ = tx.send(Err(LmfdbStreamError::Http(err)));
            return;
        }
    };

    let reader = BufReader::new(MultiGzDecoder::new(response));
    let mut discriminant = spec.initial_discriminant();

    for (idx, line_result) in reader.lines().enumerate() {
        let line_index = idx as u64 + 1;
        let raw_line = match line_result {
            Ok(line) => line,
            Err(err) => {
                let _ = tx.send(Err(LmfdbStreamError::Io(err)));
                return;
            }
        };

        let parsed = match parse_line(&raw_line, line_index) {
            Ok(parsed) => parsed,
            Err(err) => {
                let _ = tx.send(Err(err));
                return;
            }
        };

        let step = match i64::from(spec.m).checked_mul(parsed.delta) {
            Some(value) => value,
            None => {
                let _ = tx.send(Err(LmfdbStreamError::Arithmetic {
                    line_index,
                    message: "m * a overflow",
                }));
                return;
            }
        };

        discriminant = match discriminant.checked_sub(step) {
            Some(value) => value,
            None => {
                let _ = tx.send(Err(LmfdbStreamError::Arithmetic {
                    line_index,
                    message: "d_(i-1) - m*a overflow",
                }));
                return;
            }
        };

        let entry = LmfdbClassGroupEntry {
            line_index,
            delta: parsed.delta,
            discriminant,
            class_number: parsed.class_number,
            invariant_factors: parsed.invariant_factors,
        };

        if tx.send(Ok(entry)).is_err() {
            return;
        }
    }
}

struct ParsedLine {
    delta: i64,
    class_number: u64,
    invariant_factors: Vec<u64>,
}

fn parse_line(raw_line: &str, line_index: u64) -> Result<ParsedLine, LmfdbStreamError> {
    // LMFDB describes tab-separated columns: a \t b \t c1 c2 ... ct.
    // We parse that form first, then fall back to whitespace-only tokenization.
    let tab_columns: Vec<&str> = raw_line
        .split('\t')
        .map(str::trim)
        .filter(|col| !col.is_empty())
        .collect();

    if tab_columns.len() >= 2 {
        return parse_from_columns(raw_line, line_index, &tab_columns);
    }

    let ws_columns: Vec<&str> = raw_line.split_whitespace().collect();
    parse_from_columns(raw_line, line_index, &ws_columns)
}

fn parse_from_columns(
    raw_line: &str,
    line_index: u64,
    columns: &[&str],
) -> Result<ParsedLine, LmfdbStreamError> {
    let Some(delta_raw) = columns.first().copied() else {
        return Err(LmfdbStreamError::Parse {
            line_index,
            raw_line: raw_line.to_string(),
            message: "missing delta field",
        });
    };

    let Some(class_number_raw) = columns.get(1).copied() else {
        return Err(LmfdbStreamError::Parse {
            line_index,
            raw_line: raw_line.to_string(),
            message: "missing class number field",
        });
    };

    let delta = delta_raw
        .parse::<i64>()
        .map_err(|_| LmfdbStreamError::Parse {
            line_index,
            raw_line: raw_line.to_string(),
            message: "invalid integer in delta field",
        })?;

    let class_number = class_number_raw
        .parse::<u64>()
        .map_err(|_| LmfdbStreamError::Parse {
            line_index,
            raw_line: raw_line.to_string(),
            message: "invalid integer in class number field",
        })?;

    let mut invariant_factors = Vec::new();
    for column in columns.iter().skip(2) {
        for factor_raw in column.split_whitespace() {
            let factor = factor_raw
                .parse::<u64>()
                .map_err(|_| LmfdbStreamError::Parse {
                    line_index,
                    raw_line: raw_line.to_string(),
                    message: "invalid integer in invariant factors",
                })?;
            invariant_factors.push(factor);
        }
    }

    Ok(ParsedLine {
        delta,
        class_number,
        invariant_factors,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_line_matches_lmfdb_sample() {
        let first = parse_line("0 12160 380 4 4 2", 1).expect("first line should parse");
        assert_eq!(first.delta, 0);
        assert_eq!(first.class_number, 12160);
        assert_eq!(first.invariant_factors, vec![380, 4, 4, 2]);

        let second = parse_line("2 4392 2196 2", 2).expect("second line should parse");
        assert_eq!(second.delta, 2);
        assert_eq!(second.class_number, 4392);
        assert_eq!(second.invariant_factors, vec![2196, 2]);
    }

    #[test]
    fn parse_tab_separated_line() {
        let entry = parse_line("2\t4392\t2196 2", 1).expect("tab-separated line should parse");
        assert_eq!(entry.delta, 2);
        assert_eq!(entry.class_number, 4392);
        assert_eq!(entry.invariant_factors, vec![2196, 2]);
    }

    #[test]
    fn discriminant_progression_matches_lmfdb_sample() {
        let spec = LmfdbFileSpec { k: 1, r: 4, m: 16 };
        let mut d = spec.initial_discriminant();
        assert_eq!(d, -268_435_460);

        let first = parse_line("0 12160 380 4 4 2", 1).expect("first line should parse");
        d -= i64::from(spec.m) * first.delta;
        assert_eq!(d, -268_435_460);

        let second = parse_line("2 4392 2196 2", 2).expect("second line should parse");
        d -= i64::from(spec.m) * second.delta;
        assert_eq!(d, -268_435_492);
    }
}
