# rolypoly-db

Construction workflows for the external databases consumed by
[RolyPoly](https://github.com/UriNeri/rolypoly).

Generated databases, source downloads, and other large artifacts
are deliberately excluded from Git. RolyPoly remains the core/main suite's codebase, and
its `get-data` command installs bundles produced here.  
This repository is a little messy: some steps are still in notebooks or Bash
snippets, but it should work (TM...).

## Local layout

The sibling-repository layout is supported for development:

```text
rps/
├── rolypoly/
│   └── data/              # current production build/output tree
└── rolypoly-db/
```

The default data root is `../rolypoly/data`. Override it with
`--data-dir`, `ROLYPOLY_DATA`, or workflow-specific arguments.

## Environment

For local development, the sibling RolyPoly Pixi environment can run these
sources directly:

```bash
pixi shell -e dev --manifest-path=../rolypoly/pyproject.toml
PYTHONPATH=src python -m rolypoly_db.build_data --help
```

Alternatively, use this repo's manifest

```bash
pixi install
```

A separately pinned `rolypoly-tk` environment installed from Bioconda, PyPI,
or GitHub might also be used, unless the local repo of the -tk in the sister folder is more updated.

## Output contract

`manifests/rolypoly-data.json` records paths consumed by RolyPoly. Changes to
these paths must be coordinated with the RolyPoly runtime before publishing a
new data bundle.

### Rfam models and family typing

The `rfam` step pins covariance models and curated family metadata to Rfam 15.1.
Refresh the release in `download_and_extract_rfam` as one unit; do not mix
`CURRENT` models with an older family table. The existing bundle inclusion
patterns include `profiles/rfam_family.tsv.gz` and its provenance sidecar
`profiles/rfam_family.tsv.build.json`.

The headered metadata table contains `rfam_id`, `profile_name`,
`rfam_description`, `rfam_type` (Rfam's curated TP categories), `rfam_release`,
and `metadata_source`. Provenance records source URLs, SHA-256 checksums,
row count, and columns. No GO assignments are inferred. RolyPoly uses this
accession-keyed table to produce SO feature types and separate report classes.
Older bundles remain usable with conservative runtime fallbacks.
