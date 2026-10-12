# dagchainer-rs

A Rust port of [DAGchainer](https://sourceforge.net/projects/dagchainer), identifying chains of gene pairs sharing conserved order between genomic regions by finding paths through a directed acyclic graph (DAG) of pairwise matches.

This is both a standalone tool and a library crate, built as part of [CRBHits](https://github.com/kullrich/CRBHits) — a toolkit for pairwise orthology prediction.

### Original project

- Original algorithm & C++ implementation: [DAGchainer](http://dagchainer.sourceforge.net/) by A. L. Delcher and B. Haas
- Original project admins: [Joshua Orvis](https://sourceforge.net/u/jorvis/), [Brian Haas](https://sourceforge.net/u/bhaas/)
- Original project page / manual: [https://dagchainer.sourceforge.net/](https://dagchainer.sourceforge.net/)

This port reimplements the core chaining algorithm (`dagchainer.cpp`), the orchestration pipeline (`run_DAG_chainer.pl`), and all twelve accessory Perl scripts natively in Rust — with **no dependency on Perl, no subprocess calls between pipeline stages, and no intermediate temp files**.

---

## Installation

### From source

```bash
git clone https://github.com/kullrich/dagchainer-rs.git
cd dagchainer-rs
cargo build --release
```

The two binaries are then available at:
```bash
target/release/dagchainer   # core chainer (drop-in CLI replacement for the original C++ tool)
target/release/dagtools     # pipeline + all accessory tools, as subcommands
```

#### Requirements
- Rust 1.85+ (install via [rustup][rustup.rs])
- Optional: gnuplot, if you want dagtools to-gnuplot to actually render plots (otherwise use --no-exec to just generate the data/script files)

---

## Usage

### 1. The core chainer: `dagchainer`

A direct, drop-in replacement for the original compiled C++ `dagchainer` binary. Takes a 4-column match file (`pair_id x y score`) and reports chains of co-linear matches under an affine gap penalty model.

```
dagchainer \
  -G 10000 \
  -O 0 \
  -E -3 \
  -S 30 \
  -D 200000 \
  -F matches.input \
  [-r]
```

|**Flag**|**Meaning**|
|----|-------|
|-G|Length (bp) of a single gap|
|-O|Gap open penalty (must be <= 0)|
|-E|Gap extension penalty (must be <= 0)|
|-S|Minimum alignment (chain) score to report|
|-D|Maximum distance (bp) allowed between two chained matches|
|-F|Input match file|
|-r|Reverse the second coordinate list (e.g. minus-strand synteny)|

### 2. The full pipeline: `dagtools pipeline`

Most users will want this instead — it's the Rust port of `run_DAG_chainer.pl`, taking a raw 9-column tabular match/BLAST-like file and producing full `.aligncoords` chain output directly, with **no subprocess calls, no temp files, and no Perl runtime required** (the original Perl wrapper shells out to the compiled C++ binary twice per molecule pair and parses its text output back; this port calls the chaining engine natively, in-process).

**Input format** (tab-separated, 9 columns):
```
mol_1  acc_1  end5_1  end3_1  mol_2  acc_2  end5_2  end3_2  e_value
```

```
dagtools pipeline \
  -i matches.tsv \
  -g 10000 \
  -o 0 \
  -e -3 \
  -D 200000 \
  -A 6 \
  -E 1e-5 \
  [-s] [-I] [-T] [-Z constant_score] [-x min_score]
```

|**Flag**|**Meaning**|
|----|-------|
|-i|Input match file|
|-g|Length (bp) of a single gap (default: 10000)|
|-o|Gap open penalty (default: 0)|
|-e|Gap extension penalty (default: -3)|
|-M|Maximum match score derived from e-value (default: 50)|
|-Z|Use a fixed match score instead of deriving one from e-value|
|-D|Maximum distance (bp) between chained matches (default: 200000)|
|-s|Include self-molecule (tandem) comparisons|
|-I|Ignore tandem/overlapping self-chain alignments (requires -s)|
|-T|Report only tandem alignments (implies -s)|
|-A|Minimum number of aligned pairs to report a chain (default: 6)|
|-E|Maximum e-value for an input match to be considered (default: 1e-5)|
|-x|Minimum alignment score (default: computed from -A and -M)|

Output is written to `<input>.aligncoords`, matching the original tool's naming convention, so it feeds directly into every accessory tool below.

### 3. Accessory tools: `dagtools <subcommand>`

All twelve original accessory Perl scripts are available as subcommands of the same `dagtools` binary:

```
dagtools limit-top-hits <N>                        # keep best N hits per accession
dagtools filter-repetitive <window_size>            # collapse clustered repetitive matches
dagtools find-tandem-dups -m <max_intervening>      # detect tandem gene duplicates
dagtools min-chain-len <N>                          # filter chains by minimum length
dagtools count-genes-min-len <N>                    # list unique genes in chains >= N pairs
dagtools to-chain-spans <file>                      # collapse each chain to its overall span
dagtools chain-length-ratio                          # compute per-chain length ratios (A vs B)
dagtools report-segment-dist <file>                  # midpoint distances within chains
dagtools extract-chromo-dups <matchFile> <diagFile>  # keep matches whose chromo-pair has a chain
dagtools separate-datasets <matchFile>               # split matches into per-molecule-pair files
dagtools sum-nr-content <chain_spans_file> [-v]       # non-redundant bp content covered by chains
dagtools to-gnuplot --genomeA A.fa --genomeB B.fa \
                    --aligncoords chains.aligncoords \
                    --prefix myplot [--no-exec]       # dot-plot visualization
```

Run `dagtools <subcommand> --help` for full per-tool options.

### 4. Example: end-to-end run

```
# 1. run the full chaining pipeline
dagtools pipeline -i Arabidopsis.Release5.matchList \
  -g 10000 -o 0 -e -3 -D 200000 -A 6 -E 1e-5 -s -I

# 2. inspect the resulting chains
less Arabidopsis.Release5.matchList.aligncoords

# 3. filter to chains with at least 4 aligned pairs
dagtools min-chain-len 4 < Arabidopsis.Release5.matchList.aligncoords > filtered.aligncoords

# 4. collapse to chain spans and sum non-redundant genome content
dagtools to-chain-spans filtered.aligncoords > chain_spans.tsv
dagtools sum-nr-content chain_spans.tsv

# 5. visualize
dagtools to-gnuplot --genomeA genomeA.fa --genomeB genomeB.fa \
  --aligncoords filtered.aligncoords --prefix plot
```

### 5. As a library

`dagchainer-rs` can also be used directly as a Rust dependency (e.g. from CRBHits-rs):

```
[dependencies]
dagchainer-rs = "0.1"
```

```
use dagchainer_rs::{chain_matches, MatchPoint, Params};
use dagchainer_rs::pipeline::{parse_input_lines, run_pipeline, PipelineParams};

// low-level: chain a pre-built list of match points directly
let matches = vec![
    MatchPoint { pair_id: 0, x: 10, y: 12, score: 50.0 },
    MatchPoint { pair_id: 1, x: 20, y: 22, score: 50.0 },
];
let params = Params {
    indel_score: -3.0,
    gap_open_penalty: 0.0,
    bp_gap_size: 10_000,
    min_alignment_score: 10.0,
    max_dist_between_matches: 200_000,
    reverse_order: false,
};
let alignments = chain_matches(&matches, &params)?;

// high-level: run the full pipeline over a raw tabular match file
let raw_matches = parse_input_lines(reader)?;
let mut output = Vec::new();
run_pipeline(&raw_matches, &PipelineParams::default(), &mut output)?;
```

---

## Numerical fidelity vs. the original DAGchainer (C++)

This port has been validated against the original compiled `dagchainer`
binary across a real-world dataset (Arabidopsis thaliana, intra- and
inter-chromosomal comparisons), and produces **identical chain output**
in the overwhelming majority of cases, including exact match on:

- Chain membership (which matches get grouped into which alignment)
- Chain scores, to one decimal place
- The final integer `dag_position_score` reported per aligned pair

### A rare, deliberate exception: near-tie precision

There is one known class of edge case where this port's output can
differ from the original C++ binary: **when two or more candidate
predecessor chains score within floating-point rounding error of each
other.**

The original `dagchainer.cpp` stores its dynamic-programming scores
using 32-bit `float`, with precision loss compounding at every step of
every chain. This Rust port instead uses 64-bit `f64` throughout, for
better numerical accuracy and headroom against overflow/precision loss
on long chains or extreme scoring parameters.

In the vast majority of cases this makes no observable difference — the
competing chain scores aren't close enough for precision to matter. But
in rare cases (most commonly: a query gene with two or more tandemly
duplicated, nearly-equally-good targets in the other genome), the two
implementations can select a **different winner** at what is, in
effect, a numerical coin-flip. A concrete example observed in real data:

```
# original (float32):
1  At1g35470  13054900  13051614  4  At4g09310  5902946  5905237  5.4e-81   210

# this port (f64):
1  At1g35470  13054900  13051614  4  At4g09340  5927594  5925012  2.9e-160  204
```

Both `At4g09310` and `At4g09340` are valid, closely-scoring tandem
duplicates on chromosome 4; the two implementations simply broke the
tie differently. Neither choice is "more correct" — the true underlying
chain scores genuinely differ by less than either floating-point
representation can reliably distinguish.

### Why this trade-off was chosen deliberately

Bit-for-bit reproducibility of the original's specific tie-breaking
behavior would require replicating its 32-bit float arithmetic exactly,
including narrowing every intermediate score back to `float` precision
at each DP step. This was prototyped and confirmed to resolve the
discrepancy above, but was ultimately **not adopted**, in favor of
keeping this implementation on `f64` throughout — prioritizing improved
numerical headroom over exact replication of a legacy precision
limitation.

**Practical implications for users:**

- If you are validating this port against existing, previously
  published results generated with the original C++ `dagchainer`
  binary, expect output to match exactly *except* in rare near-tie
  cases like the one above — these are not bugs, but a documented,
  intentional consequence of using higher-precision arithmetic.
- These discrepancies only ever affect **which specific tandem/near-
  equal match is chosen** in a tie; they do not affect overall chain
  detection, chain boundaries, or scores for the (overwhelming majority
  of) non-tied matches in a chain.
- If exact, bit-for-bit parity with the original binary is a hard
  requirement for your use case (e.g. reproducing a specific published
  analysis verbatim), please open an issue — a `--legacy-float32`
  compatibility mode could be added as an opt-in flag in a future
  release.

---

## Other documented deviations from the original Perl/C++ tools

A handful of other small, deliberate corrections were made relative to the original scripts during porting (full detail in each module's doc comments):
- `limit_to_top_N_hits.pl`'s secondary sort was inverted (kept worst e-values instead of best) — fixed to keep the best.
- Several scripts' `my $x = $ARGV[0] or die` pattern rejected an explicit `0` due to Perl truthiness — all numeric CLI arguments in this port accept `0` literally where it's a sensible value.
- `dagchainer_to_gnuplot.pl`'s `--accession` filter only checked one side of a chain (`accA`) despite its own help text describing both — fixed to check both sides.
- Non-deterministic Perl hash key iteration (randomized since Perl 5.18) in a few scripts' output ordering has been replaced with deterministic sorted ordering.

See individual module documentation (`cargo doc --open`) for full details, including cases where an apparent "bug" was deliberately **preserved** for pipeline compatibility rather than fixed (e.g. `dagchainer_to_gnuplot.pl`'s scaffolding population order).

---

## Manual

For algorithmic background and the original tool's manual, see [https://dagchainer.sourceforge.net/](https://dagchainer.sourceforge.net/)

---

## License

MIT

---

## References

Haas B. J., Delcher A. L., Wortman J. R., and Salzberg S. L. (2004). **DAGchainer: a tool for mining segmental genome duplications and synteny.** *Bioinformatics*, **20(18)**, 3643-3646. [https://doi.org/10.1093/bioinformatics/bth397](https://doi.org/10.1093/bioinformatics/bth397)
