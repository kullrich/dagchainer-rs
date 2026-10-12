//! Shared data types for visualization tools, replacing the
//! `SequenceFeature` / `DatabaseMatch` classes from the original
//! `JsyntenyView.java`.

/// A single match/alignment entry between two molecules, equivalent to
/// one `DatabaseMatch` in the Java original.
#[derive(Debug, Clone)]
pub struct MatchRecord {
    pub mol_a: String,
    pub acc_a: String,
    pub end5_a: i64,
    pub end3_a: i64,

    pub mol_b: String,
    pub acc_b: String,
    pub end5_b: i64,
    pub end3_b: i64,

    pub score: Option<f64>,

    /// Which chain/alignment ("diagonal" in the Java original) this match
    /// belongs to, if known. `None` = unassigned / not part of any chain.
    pub diagonal: Option<usize>,
}

impl MatchRecord {
    pub fn mid_a(&self) -> f64 {
        (self.end5_a + self.end3_a) as f64 / 2.0
    }

    pub fn mid_b(&self) -> f64 {
        (self.end5_b + self.end3_b) as f64 / 2.0
    }

    /// Canonical, order-independent key for grouping by molecule pair,
    /// mirroring `get_chromo_pair_key` in the Java original.
    pub fn molecule_pair_key(&self) -> (String, String) {
        if self.mol_a <= self.mol_b {
            (self.mol_a.clone(), self.mol_b.clone())
        } else {
            (self.mol_b.clone(), self.mol_a.clone())
        }
    }
}