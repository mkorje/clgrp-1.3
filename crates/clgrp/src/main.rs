use clgrp::{LmfdbClassGroupEntry, LmfdbFileSpec, stream_lmfdb_class_groups};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use std::env;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::ExitCode;
use std::thread;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

fn main() -> ExitCode {
    let args: Vec<String> = env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("fetch-test") => match run_fetch_test(&args[2..]) {
            Ok(()) => ExitCode::SUCCESS,
            Err(msg) => {
                eprintln!("{msg}");
                ExitCode::from(1)
            }
        },
        Some("build-bjt-fixture") => match run_build_bjt_fixture(&args[2..]) {
            Ok(()) => ExitCode::SUCCESS,
            Err(msg) => {
                eprintln!("{msg}");
                ExitCode::from(1)
            }
        },
        Some("build-bjt-fixture-from-files") => {
            match run_build_bjt_fixture_from_files(&args[2..]) {
                Ok(()) => ExitCode::SUCCESS,
                Err(msg) => {
                    eprintln!("{msg}");
                    ExitCode::from(1)
                }
            }
        }
        _ => {
            print_usage(&args[0]);
            ExitCode::from(2)
        }
    }
}

fn run_fetch_test(args: &[String]) -> Result<(), String> {
    if args.len() < 3 {
        return Err(
            "fetch-test needs at least 3 positional args: <k> <r> <m> [--max N] [--poll-ms MS]"
                .to_string(),
        );
    }

    let k = parse_u32(&args[0], "k")?;
    let r = parse_u32(&args[1], "r")?;
    let m = parse_u32(&args[2], "m")?;

    let mut max_rows: usize = 20;
    let mut poll_ms: u64 = 25;
    let mut i = 3;
    while i < args.len() {
        match args[i].as_str() {
            "--max" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --max".to_string())?;
                max_rows = value
                    .parse::<usize>()
                    .map_err(|_| format!("invalid --max value: {value}"))?;
            }
            "--poll-ms" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --poll-ms".to_string())?;
                poll_ms = value
                    .parse::<u64>()
                    .map_err(|_| format!("invalid --poll-ms value: {value}"))?;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }

    let spec = LmfdbFileSpec { k, r, m };
    let mut stream =
        stream_lmfdb_class_groups(spec).map_err(|e| format!("failed to start stream: {e}"))?;

    println!(
        "stream started for cl{}mod{}.{} (printing up to {} rows)",
        r, m, k, max_rows
    );
    println!("line\tdelta\tdiscriminant\tclass_number\tinvariant_factors");

    let mut printed = 0usize;
    loop {
        if printed >= max_rows {
            break;
        }

        match stream.try_next() {
            Some(Ok(entry)) => {
                println!(
                    "{}\t{}\t{}\t{}\t{:?}",
                    entry.line_index,
                    entry.delta,
                    entry.discriminant,
                    entry.class_number,
                    entry.invariant_factors
                );
                printed += 1;
            }
            Some(Err(err)) => return Err(format!("stream failed: {err}")),
            None => {
                thread::sleep(Duration::from_millis(poll_ms));
            }
        }
    }

    Ok(())
}

fn parse_u32(raw: &str, name: &str) -> Result<u32, String> {
    raw.parse::<u32>()
        .map_err(|_| format!("invalid {name} value: {raw}"))
}

fn parse_u64(raw: &str, name: &str) -> Result<u64, String> {
    raw.parse::<u64>()
        .map_err(|_| format!("invalid {name} value: {raw}"))
}

fn parse_usize(raw: &str, name: &str) -> Result<usize, String> {
    raw.parse::<usize>()
        .map_err(|_| format!("invalid {name} value: {raw}"))
}

#[derive(Debug, Clone)]
struct FixtureRow {
    spec: LmfdbFileSpec,
    source_file: Option<String>,
    line_index: u64,
    delta: i64,
    discriminant: i64,
    class_number: u64,
    invariant_factors: Vec<u64>,
}

impl FixtureRow {
    fn from_entry(spec: LmfdbFileSpec, entry: LmfdbClassGroupEntry) -> Self {
        Self {
            spec,
            source_file: None,
            line_index: entry.line_index,
            delta: entry.delta,
            discriminant: entry.discriminant,
            class_number: entry.class_number,
            invariant_factors: entry.invariant_factors,
        }
    }
}

#[derive(Debug, Clone)]
struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    fn new(seed: u64) -> Self {
        // Avoid the all-zero absorbing state.
        let state = if seed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            seed
        };
        Self { state }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn gen_bounded(&mut self, upper_exclusive: u64) -> u64 {
        if upper_exclusive == 0 {
            return 0;
        }
        self.next_u64() % upper_exclusive
    }
}

fn reservoir_sample_file(
    spec: LmfdbFileSpec,
    sample_size: usize,
    scan_limit: usize,
    rng: &mut XorShift64,
) -> Result<Vec<LmfdbClassGroupEntry>, String> {
    let mut stream =
        stream_lmfdb_class_groups(spec).map_err(|e| format!("failed to start stream: {e}"))?;
    let mut sampled = Vec::with_capacity(sample_size);
    let mut seen = 0_u64;

    for next in &mut stream {
        if seen >= scan_limit as u64 {
            break;
        }
        let row = next.map_err(|e| format!("stream failed: {e}"))?;
        seen += 1;

        if sampled.len() < sample_size {
            sampled.push(row);
        } else {
            let j = rng.gen_bounded(seen);
            if j < sample_size as u64 {
                sampled[j as usize] = row;
            }
        }
    }

    if seen < sample_size as u64 {
        return Err(format!(
            "only saw {seen} rows (scan_limit={scan_limit}) for cl{}mod{}.{} but need {sample_size}",
            spec.r, spec.m, spec.k,
        ));
    }

    sampled.sort_by_key(|row| row.line_index);
    Ok(sampled)
}

fn write_fixture_file(
    path: &PathBuf,
    seed: u64,
    scan_limit: usize,
    k_min: u32,
    k_max: u32,
    k_files_per_class: usize,
    rows: &[FixtureRow],
) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "failed to create parent directory {}: {e}",
                parent.display()
            )
        })?;
    }

    let generated_unix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| format!("failed to get UNIX timestamp: {e}"))?
        .as_secs();

    let mut body = String::new();
    body.push_str("# LMFDB BJT fixture generated by `clgrp build-bjt-fixture`\n");
    body.push_str(&format!("# seed={seed}\n"));
    body.push_str(&format!("# scan_limit={scan_limit}\n"));
    body.push_str(&format!("# k_min={k_min}\n"));
    body.push_str(&format!("# k_max={k_max}\n"));
    body.push_str(&format!("# k_files_per_class={k_files_per_class}\n"));
    body.push_str(&format!("# generated_unix={generated_unix}\n"));
    body.push_str(
        "# source_file\tr\tm\tk\tline\tdelta\tdiscriminant\tclass_number\tinvariant_factors_csv\n",
    );

    for row in rows {
        let invariant_csv = row
            .invariant_factors
            .iter()
            .map(u64::to_string)
            .collect::<Vec<_>>()
            .join(",");
        let source_file = row.source_file.as_deref().unwrap_or("");
        body.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            source_file,
            row.spec.r,
            row.spec.m,
            row.spec.k,
            row.line_index,
            row.delta,
            row.discriminant,
            row.class_number,
            invariant_csv
        ));
    }

    fs::write(path, body).map_err(|e| format!("failed to write {}: {e}", path.display()))
}

fn run_build_bjt_fixture(args: &[String]) -> Result<(), String> {
    let mut samples_per_file = 1_000_usize;
    let mut scan_limit = 40_000_usize;
    let mut k_min = 0_u32;
    let mut k_max = 4_095_u32;
    let mut k_files_per_class = 8_usize;
    let mut seed = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| format!("failed to derive seed from system time: {e}"))?
        .as_nanos() as u64;
    let mut out_path = PathBuf::from("crates/clgrp/testdata/lmfdb_bjt_samples.tsv");

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--samples-per-file" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --samples-per-file".to_string())?;
                samples_per_file = parse_usize(value, "samples-per-file")?;
            }
            "--k-min" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --k-min".to_string())?;
                k_min = parse_u32(value, "k-min")?;
            }
            "--k-max" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --k-max".to_string())?;
                k_max = parse_u32(value, "k-max")?;
            }
            "--k-files-per-class" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --k-files-per-class".to_string())?;
                k_files_per_class = parse_usize(value, "k-files-per-class")?;
            }
            "--scan-limit" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --scan-limit".to_string())?;
                scan_limit = parse_usize(value, "scan-limit")?;
            }
            "--seed" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --seed".to_string())?;
                seed = parse_u64(value, "seed")?;
            }
            "--out" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --out".to_string())?;
                out_path = PathBuf::from(value);
            }
            other => {
                return Err(format!("unknown argument: {other}"));
            }
        }
        i += 1;
    }

    if k_min > k_max {
        return Err(format!("k-min ({k_min}) must be <= k-max ({k_max})"));
    }
    if k_files_per_class == 0 {
        return Err("k-files-per-class must be >= 1".to_string());
    }
    let k_space = (u64::from(k_max) - u64::from(k_min) + 1) as usize;
    if k_files_per_class > k_space {
        return Err(format!(
            "k-files-per-class ({k_files_per_class}) exceeds k-range size ({k_space})"
        ));
    }

    let base_specs: [(u32, u32); 4] = [(3, 8), (4, 16), (7, 8), (8, 16)];
    let mut rng = XorShift64::new(seed);
    let mut fixture_rows = Vec::with_capacity(samples_per_file * base_specs.len());

    for (r, m) in base_specs {
        let class_k_files = k_files_per_class;
        eprintln!(
            "Sampling {samples_per_file} rows for cl{r}mod{m} across {class_k_files} k-files in [{k_min}, {k_max}] (scan_limit={scan_limit})..."
        );

        let base_rows = samples_per_file / class_k_files;
        let extra_rows = samples_per_file % class_k_files;
        let range_size = u64::from(k_max) - u64::from(k_min) + 1;

        for bucket in 0..class_k_files {
            let start_off = (bucket as u64 * range_size) / class_k_files as u64;
            let end_off = ((bucket as u64 + 1) * range_size) / class_k_files as u64;
            let span = end_off.saturating_sub(start_off).max(1);
            let bucket_start = u64::from(k_min) + start_off;
            let rows_for_k = base_rows + usize::from(bucket < extra_rows);
            let k = (bucket_start + rng.gen_bounded(span)) as u32;
            let spec = LmfdbFileSpec { k, r, m };
            eprintln!("  sampling {rows_for_k} rows from cl{r}mod{m}.{k}");

            let attempts = 4_u8;
            let mut sampled = None;
            for attempt in 1..=attempts {
                match reservoir_sample_file(spec, rows_for_k, scan_limit, &mut rng) {
                    Ok(rows) => {
                        sampled = Some(rows);
                        break;
                    }
                    Err(err) => {
                        eprintln!(
                            "    attempt {attempt}/{attempts} failed for cl{r}mod{m}.{k}: {err}"
                        );
                        thread::sleep(Duration::from_millis(120));
                    }
                }
            }

            let rows = sampled.ok_or_else(|| {
                format!("unable to sample cl{r}mod{m}.{k} after multiple attempts")
            })?;
            for row in rows {
                fixture_rows.push(FixtureRow::from_entry(spec, row));
            }
            thread::sleep(Duration::from_millis(120));
        }
        eprintln!("  completed cl{r}mod{m}");
    }

    if fixture_rows.len() != samples_per_file * base_specs.len() {
        return Err(format!(
            "internal error: expected {} rows but produced {}",
            samples_per_file * base_specs.len(),
            fixture_rows.len()
        ));
    }

    fixture_rows.sort_by_key(|row| (row.spec.r, row.spec.m, row.spec.k, row.line_index));
    write_fixture_file(
        &out_path,
        seed,
        scan_limit,
        k_min,
        k_max,
        k_files_per_class,
        &fixture_rows,
    )?;
    println!(
        "Wrote {} sampled rows to {}",
        fixture_rows.len(),
        out_path.display()
    );
    Ok(())
}

fn parse_filename_spec(name: &str) -> Result<LmfdbFileSpec, String> {
    let base = name
        .strip_suffix(".gz")
        .ok_or_else(|| format!("expected .gz suffix in file name: {name}"))?;
    let (prefix, k_raw) = base
        .rsplit_once('.')
        .ok_or_else(|| format!("expected '<prefix>.<k>.gz' format in file name: {name}"))?;
    let k = k_raw
        .parse::<u32>()
        .map_err(|_| format!("invalid k in file name: {name}"))?;
    let prefix = prefix
        .strip_prefix("cl")
        .ok_or_else(|| format!("expected 'cl' prefix in file name: {name}"))?;
    let (r_raw, m_raw) = prefix
        .split_once("mod")
        .ok_or_else(|| format!("expected 'cl<r>mod<m>' prefix in file name: {name}"))?;
    let r = r_raw
        .parse::<u32>()
        .map_err(|_| format!("invalid r in file name: {name}"))?;
    let m = m_raw
        .parse::<u32>()
        .map_err(|_| format!("invalid m in file name: {name}"))?;
    Ok(LmfdbFileSpec { k, r, m })
}

fn parse_local_line(raw_line: &str) -> Result<(i64, u64, Vec<u64>), String> {
    let cols: Vec<&str> = raw_line.split_whitespace().collect();
    if cols.len() < 3 {
        return Err(format!("expected at least 3 columns, got: {raw_line}"));
    }

    let delta = cols[0]
        .parse::<i64>()
        .map_err(|_| format!("invalid delta: {raw_line}"))?;
    let class_number = cols[1]
        .parse::<u64>()
        .map_err(|_| format!("invalid class number: {raw_line}"))?;
    let mut invariants = Vec::with_capacity(cols.len() - 2);
    for col in cols.iter().skip(2) {
        invariants.push(
            col.parse::<u64>()
                .map_err(|_| format!("invalid invariant in line: {raw_line}"))?,
        );
    }

    Ok((delta, class_number, invariants))
}

fn reservoir_sample_local_file(
    path: &PathBuf,
    spec: LmfdbFileSpec,
    sample_size: usize,
    rng: &mut XorShift64,
) -> Result<Vec<FixtureRow>, String> {
    let file = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let mut discriminant = spec.initial_discriminant();
    let mut sampled = Vec::with_capacity(sample_size);
    let mut seen = 0_u64;
    let source_name = path
        .file_name()
        .and_then(|x| x.to_str())
        .ok_or_else(|| format!("non-utf8 file name: {}", path.display()))?
        .to_string();

    for (idx, line_result) in reader.lines().enumerate() {
        let raw_line = line_result.map_err(|e| format!("read {}: {e}", path.display()))?;
        if raw_line.trim().is_empty() {
            continue;
        }

        let line_index = idx as u64 + 1;
        let (delta, class_number, invariant_factors) = parse_local_line(&raw_line)?;
        let step = i64::from(spec.m)
            .checked_mul(delta)
            .ok_or_else(|| format!("overflow computing step in {}", path.display()))?;
        discriminant = discriminant
            .checked_sub(step)
            .ok_or_else(|| format!("overflow updating discriminant in {}", path.display()))?;

        let row = FixtureRow {
            spec,
            source_file: Some(source_name.clone()),
            line_index,
            delta,
            discriminant,
            class_number,
            invariant_factors,
        };

        seen += 1;
        if sampled.len() < sample_size {
            sampled.push(row);
        } else {
            let j = rng.gen_bounded(seen);
            if j < sample_size as u64 {
                sampled[j as usize] = row;
            }
        }
    }

    if seen < sample_size as u64 {
        return Err(format!(
            "{} has only {seen} rows; need {sample_size}",
            path.display()
        ));
    }

    sampled.sort_by_key(|row| row.line_index);
    Ok(sampled)
}

fn run_build_bjt_fixture_from_files(args: &[String]) -> Result<(), String> {
    let mut input_dir = PathBuf::from("files");
    let mut output_path = PathBuf::from("crates/clgrp/testdata/lmfdb_bjt_local_samples.tsv.gz");
    let mut sample_size = 50_000_usize;
    let mut seed = 20260217_u64;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--input-dir" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --input-dir".to_string())?;
                input_dir = PathBuf::from(value);
            }
            "--out" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --out".to_string())?;
                output_path = PathBuf::from(value);
            }
            "--samples-per-file" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --samples-per-file".to_string())?;
                sample_size = parse_usize(value, "samples-per-file")?;
            }
            "--seed" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "missing value for --seed".to_string())?;
                seed = parse_u64(value, "seed")?;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }

    let mut files = fs::read_dir(&input_dir)
        .map_err(|e| format!("read_dir {}: {e}", input_dir.display()))?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.extension().and_then(|x| x.to_str()) == Some("gz"))
        .collect::<Vec<_>>();
    files.sort();

    if files.is_empty() {
        return Err(format!("no .gz files found in {}", input_dir.display()));
    }

    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("create_dir_all {}: {e}", parent.display()))?;
    }
    let out_file =
        File::create(&output_path).map_err(|e| format!("create {}: {e}", output_path.display()))?;
    let mut out = GzEncoder::new(out_file, Compression::default());
    let generated_unix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| format!("failed to get UNIX timestamp: {e}"))?
        .as_secs();

    writeln!(
        out,
        "# LMFDB local BJT fixture generated from downloaded files"
    )
    .map_err(|e| format!("write fixture header: {e}"))?;
    writeln!(out, "# seed={seed}").map_err(|e| format!("write fixture header: {e}"))?;
    writeln!(out, "# samples_per_file={sample_size}")
        .map_err(|e| format!("write fixture header: {e}"))?;
    writeln!(out, "# source_dir={}", input_dir.display())
        .map_err(|e| format!("write fixture header: {e}"))?;
    writeln!(out, "# generated_unix={generated_unix}")
        .map_err(|e| format!("write fixture header: {e}"))?;
    writeln!(
        out,
        "# source_file\tr\tm\tk\tline\tdelta\tdiscriminant\tclass_number\tinvariant_factors_csv"
    )
    .map_err(|e| format!("write fixture header: {e}"))?;

    let mut rng = XorShift64::new(seed);
    for path in &files {
        let file_name = path
            .file_name()
            .and_then(|x| x.to_str())
            .ok_or_else(|| format!("non-utf8 file name: {}", path.display()))?;
        let spec = parse_filename_spec(file_name)?;
        eprintln!("sampling {} rows from {}", sample_size, file_name);

        let sampled = reservoir_sample_local_file(path, spec, sample_size, &mut rng)?;
        for row in sampled {
            let invariant_csv = row
                .invariant_factors
                .iter()
                .map(u64::to_string)
                .collect::<Vec<_>>()
                .join(",");
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.source_file.unwrap_or_default(),
                row.spec.r,
                row.spec.m,
                row.spec.k,
                row.line_index,
                row.delta,
                row.discriminant,
                row.class_number,
                invariant_csv
            )
            .map_err(|e| format!("write fixture row: {e}"))?;
        }
    }

    out.finish()
        .map_err(|e| format!("finalize {}: {e}", output_path.display()))?;

    println!(
        "Wrote sampled fixture from {} files to {}",
        files.len(),
        output_path.display()
    );
    Ok(())
}

fn print_usage(bin: &str) {
    eprintln!("Usage:");
    eprintln!("  {bin} fetch-test <k> <r> <m> [--max N] [--poll-ms MS]");
    eprintln!(
        "  {bin} build-bjt-fixture [--samples-per-file N] [--scan-limit L] [--k-min A] [--k-max B] [--k-files-per-class C] [--seed S] [--out PATH]"
    );
    eprintln!(
        "  {bin} build-bjt-fixture-from-files [--input-dir DIR] [--samples-per-file N] [--seed S] [--out PATH]"
    );
    eprintln!();
    eprintln!("Example:");
    eprintln!("  {bin} fetch-test 1 4 16 --max 10");
    eprintln!(
        "  {bin} build-bjt-fixture --samples-per-file 1000 --scan-limit 40000 --k-min 0 --k-max 4095 --k-files-per-class 8 --seed 20260217"
    );
    eprintln!(
        "  {bin} build-bjt-fixture-from-files --input-dir files --samples-per-file 50000 --seed 20260217 --out crates/clgrp/testdata/lmfdb_bjt_local_samples.tsv.gz"
    );
}
