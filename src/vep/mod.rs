//! VEP compute exposed as DuckDB functions, wrapping the fastVEP engine
//! (consequence / HGVS / ACMG) for exact parity. See docs/DESIGN.md §3.

pub mod annotate;
pub mod consequence;
mod engine;
pub mod normalize;
pub mod tcache;
