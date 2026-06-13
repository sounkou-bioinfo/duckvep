//! VEP compute exposed as DuckDB functions, wrapping the fastVEP engine
//! (consequence / HGVS / ACMG) for exact parity. See DESIGN.md §3.

pub mod annotate;
mod engine;
