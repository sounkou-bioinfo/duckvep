//! Small helpers for writing DuckDB output vectors.

use duckdb::core::{DataChunkHandle, Inserter};

/// Fill a `LIST<VARCHAR>` output column from a per-row slice of strings.
///
/// `rows` is the slice of source rows for this chunk; `get` extracts the list
/// of strings for one row. Writes the flattened child vector, the per-row
/// `(offset, length)` entries, and the child length.
pub(crate) fn fill_string_list<T>(
    output: &mut DataChunkHandle,
    col: usize,
    rows: &[T],
    get: impl Fn(&T) -> &[String],
) {
    let total: usize = rows.iter().map(|r| get(r).len()).sum();
    let mut list = output.list_vector(col);
    {
        let child = list.child(total.max(1));
        let mut off = 0usize;
        for r in rows {
            for (j, s) in get(r).iter().enumerate() {
                child.insert(off + j, s.as_str());
            }
            off += get(r).len();
        }
    }
    let mut off = 0usize;
    for (i, r) in rows.iter().enumerate() {
        let len = get(r).len();
        list.set_entry(i, off, len);
        off += len;
    }
    list.set_len(total);
}
