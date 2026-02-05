# ℓ-adic Growth Analyzer

Analyzes class group tabulations to find discriminants where a ℤ/ℓ^N factor grows to ℤ/ℓ^(N+1) when passing from the maximal order to an order of index ℓ².

Automatically processes all four congruence classes for fundamental discriminants:
- 8 mod 16
- 4 mod 16
- 3 mod 8
- 7 mod 8

## Building

```bash
cargo build --release
```

## Usage

```bash
./target/release/ell_growth \
    --folder /path/to/data \
    --ell 3 \
    -D 100000000 \
    --files 100 \
    -N 1 \
    [--mode strict|any|net] \
    [--verbose]
```

### Arguments

| Arg | Description |
|-----|-------------|
| `--folder` | Base folder containing `cl[a]mod[m]/` and `cl[a]mod[m]l[ell]/` directories |
| `--ell` | Prime ℓ |
| `-D` | Maximum \|discriminant\| |
| `--files` | Number of files |
| `-N` | Target exponent: find ℓ^N → ℓ^(N+1) growth |
| `--mode` | Growth detection mode (see below) |
| `--verbose` | Print each matching discriminant |

### Growth Detection Modes

- **strict** (default): A ℓ^N factor disappears AND a ℓ^(N+1) factor appears
  - `count(ℓ^{N+1})_ell > count(ℓ^{N+1})_fund` AND `count(ℓ^N)_ell < count(ℓ^N)_fund`
  
- **any**: Fundamental has ℓ^N factor AND ell has ℓ^(N+1) factor
  - Less restrictive: just checks presence, not strict growth
  
- **net**: Net increase in ℓ^(N+1) factors
  - `count(ℓ^{N+1})_ell > count(ℓ^{N+1})_fund`

### Example

Find all discriminants where a ℤ/3ℤ factor grows to ℤ/9ℤ:

```bash
./target/release/ell_growth \
    --folder ./data \
    --ell 3 \
    -D 10000000 \
    --files 10 \
    -N 1
```

## File Formats

### Fundamental (cl[a]mod[m].[index].gz)

```
dist h c1 c2 ... ct
```

- `dist`: distance (D_{i+1} = D_i + dist × m)
- `h`: class number
- `c1...ct`: invariant factors (Smith Normal Form)

### Index ℓ² (cl[a]mod[m]l[ell].[index].gz)

```
dist kron c1 c2 ... ct
```

- `dist`: distance (same as fundamental)
- `kron`: Kronecker symbol (D₀/ℓ): -1 (inert), 0 (ramified), 1 (split)
- `c1...ct`: invariant factors of Cl(𝒪_f)

## Output

```
ℓ-adic growth analysis
======================
folder: "./data"
ℓ=3
D_max=10000000, files=10
Target: ℤ/3^1ℤ → ℤ/3^2ℤ growth
Detection mode: strict

Processing 8 mod 16 ...
Processing 4 mod 16 ...
Processing 3 mod 8 ...
Processing 7 mod 8 ...

Results by congruence class
===========================

8 mod 16:
  Total discriminants: 114197
  With ℤ/3^1ℤ factor: 3456 (3.0263%)
  With growth to ℤ/3^2ℤ: 234 (0.2049% of total, 6.7708% of those with factor)
    kron=-1 (inert): with_factor=1152, with_growth=78 (6.77%)
    kron= 0 (ramified): with_factor=1152, with_growth=78 (6.77%)
    kron= 1 (split): with_factor=1152, with_growth=78 (6.77%)

...

Grand Total (all congruence classes)
====================================
Total discriminants: 456789
With ℤ/3^1ℤ factor: 12345 (2.7022%)
With growth to ℤ/3^2ℤ: 789 (0.1727% of total, 6.3932% of those with factor)

Breakdown by Kronecker symbol (D₀/ℓ):
  kron=-1 (inert): with_factor=4115, with_growth=263 (6.39%)
  kron= 0 (ramified): with_factor=4115, with_growth=263 (6.39%)
  kron= 1 (split): with_factor=4115, with_growth=263 (6.39%)
```
