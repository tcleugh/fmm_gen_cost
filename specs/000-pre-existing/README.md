# 000 — Pre-Existing Inputs and Outputs

This folder documents the input/output contracts the FMM project inherited from upstream — formats and CLI flags that pre-date any spec-kit feature in this fork. They are written down here so that:

- Feature specs under `specs/NNN-feature-name/` can reference them by name instead of re-describing the same shapefile / CSV layouts in every contract.
- Anyone adding a new matcher or output column has a single canonical reference for the conventions the rest of the codebase already follows.

These are *descriptive* contracts (what the code does today), not *normative* ones (what some external spec requires). When the C++ implementation and this folder disagree, the C++ wins; please update the docs.

## Contracts

| File | Subject |
|---|---|
| [contracts/network-shapefile.md](contracts/network-shapefile.md) | Road network input (ESRI shapefile / GeoPackage / GeoJSON) |
| [contracts/turn-ban-file.md](contracts/turn-ban-file.md) | CSV of banned (in_edge, out_edge) transitions |
| [contracts/gps-trajectory.md](contracts/gps-trajectory.md) | GPS trajectory input (CSV-WKT, CSV-point, or OGR file) |
| [contracts/ubodt-file.md](contracts/ubodt-file.md) | Upper Bounded Origin-Destination Table consumed by FMM |
| [contracts/match-result-csv.md](contracts/match-result-csv.md) | Base CSV layout `CSVMatchResultWriter` emits for all matchers |
| [contracts/matcher-cli.md](contracts/matcher-cli.md) | Shared CLI flags (`NetworkConfig`, `GPSConfig`, `ResultConfig`) every matcher exposes |

## What lives where

The contracts in this folder describe behavior implemented by:

- `src/network/network.cpp` — network loading + turn-ban CSV reading.
- `src/config/{network,gps,result}_config.cpp` — CLI / XML parsing for the three shared sub-configs.
- `src/io/gps_reader.cpp` — GPS trajectory reading (three formats).
- `src/io/mm_writer.cpp` — base CSV result writer.
- `src/mm/fmm/ubodt.cpp` + `ubodt_gen_algorithm.cpp` — UBODT file format + generator.

When polymatch (specs/001) or any future feature adds new fields / formats, this folder stays untouched — the feature's own `contracts/` adds *delta* documentation.
