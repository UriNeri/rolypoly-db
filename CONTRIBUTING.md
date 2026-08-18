# Contributing to rolypoly-db

This repository builds external data consumed by RolyPoly. It follows the
RolyPoly project’s GPLv3 and scientific-provenance expectations, while favoring
direct construction workflows over application-style abstraction.

## Principles
Similar to RolyPoly's [contributing guidelines](https://github.com/UriNeri/rolypoly/blob/main/CONTRIBUTING.md), but slightly different. Briefly (especially if you are an LLM)..:

- Prefer one readable, reproducible, resumable construction workflow over
  interactive notebooks or layers of single-use wrappers.
- Bash is acceptable for orchestration when it makes parallelism, compression,
  or standard bioinformatics tools clearer.
- Keep paths, threads, memory, and temporary directories configurable.
- Use Polars for tabular files and SQLite queries, including its streaming and
  database-reading facilities. Do not write row-by-row parsers when an existing
  dependency already provides the operation.
- Reuse pinned RolyPoly utilities instead of copying FASTX, logging, download,
  profile-building, or command-runner implementations.
- Keep the primary Click orchestration in `build_data()` when extracting it
  would only create a pass-through function. Avoid underscore-prefixed helpers,
  decorative character separators, and functions used once in one fixed way.
  Put durable explanations in useful function docstrings.
- Never commit downloaded references, generated databases, or production logs.
- Record source versions, URLs, checksums, parameters, and output schemas.
- Do not silently discard unmatched identifiers or unresolved taxonomy. (this may be a bug or real 'deleted.dmp' or similar entries).
- Store useful intermediate work in an ignored, resumable work directory rather
  than `/tmp`; reserve `/tmp` for disposable data that cannot matter later.

## Environment and checks
`build_data.py` is a construction script: do not create unit tests or synthetic fixtures for its internal functions by default. Validation is only needed for the mature databases/formats, and that can be done by the regular rolypoly commands on the newly created DBs etc. Note that some reference downloads or full database builds can take long. only really worth doing if we think the remote changed compared to the local (or if no local, or if really needed).

## RolyPoly data contract
Generated bundles must follow `manifests/rolypoly-data.json`. Coordinate any
path or schema change with the RolyPoly runtime before publishing a bundle.

The production build may use a sibling checkout at `../rolypoly` and write to
`../rolypoly/data`, but committed workflows must not require that absolute local
layout.
