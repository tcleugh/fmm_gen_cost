# Contract: Turn Ban File

Optional input accompanying the network shapefile (see [network-shapefile.md](network-shapefile.md)). Lists ordered pairs of edges across which routing is forbidden — used to model "no left turn" / "no U-turn" restrictions.

## Format

Plain CSV, comma-delimited (`,`). UTF-8 / ASCII text. One header line followed by data rows.

```text
in_edge_id,out_edge_id
102,103
102,201
407,512
```

- The header line is **always skipped**; column names are ignored. The first column is interpreted as the in-edge ID, the second as the out-edge ID.
- Both IDs are parsed via `std::stoi` from `Network::read_turn_ban_file`, so they must fit in `int`. Mixed-line whitespace handling is implementation-defined; keep rows tight.
- Each row inserts one directed edge-pair into `Network::turn_bans` (a hash set keyed on `Turn{in_edge_index, out_edge_index}`).

## Lookup semantics

When `LinkGraph` is built (and during STMatch's `CompositeGraph` construction), each candidate arc `(in_edge, out_edge)` is filtered by `Network::is_turn_banned(in_e, out_e)`. Banned pairs are simply omitted from the adjacency — they are not routable, even via Dijkstra.

The pair is **directional**: banning `(102, 103)` does not ban `(103, 102)`. List both rows if both directions should be banned.

## "Disabled" sentinel

The `Network` constructor's `turn_ban_file` parameter has no default in the header; in CLI / XML configs the default value is the literal string `"NO_TURN_BANS"`. When `read_turn_ban_file` sees this value it returns immediately without touching the filesystem.

## Validation

`NetworkConfig::validate()` only verifies file existence + extension when `turn_ban_file != "NO_TURN_BANS"`. There is no per-row validation:

- Unparseable IDs cause `std::stoi` to throw / leave the values from the previous row.
- IDs that don't exist in the network throw from `get_edge_index` during the insert step.
- Lines without a comma are silently treated as a single-field row (the second field never gets assigned).

If you need stricter validation, do it upstream of the matcher.
