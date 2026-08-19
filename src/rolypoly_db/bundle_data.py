"""Create the distributable RolyPoly data archive from the data contract."""

import fnmatch
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

from rich_click import command, option


def manifest_default_path() -> Path:
    return Path(__file__).resolve().parents[2] / "manifests" / "rolypoly-data.json"


def normalize_relative(path: Path) -> str:
    return path.as_posix()


def pattern_matches(path: str, pattern: str) -> bool:
    if pattern.endswith("/**"):
        return path == pattern[:-3] or path.startswith(pattern[:-2])
    return fnmatch.fnmatch(path, pattern)


def collect_bundle_files(data_dir: Path, manifest_path: Path) -> list[str]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    bundle = manifest.get("bundle")
    if not bundle:
        raise ValueError(f"{manifest_path} does not define a bundle section")

    include_patterns = bundle.get("include", [])
    include_small = bundle.get("include_small", {})
    include_small_max_bytes = include_small.get("max_bytes", 0)
    include_small_patterns = include_small.get("patterns", [])
    exclude_patterns = bundle.get("exclude", [])
    if not include_patterns:
        raise ValueError(f"{manifest_path} bundle section has no include patterns")

    included: set[str] = set()
    unmatched_patterns: list[str] = []
    for pattern in include_patterns:
        matches = [
            normalize_relative(path.relative_to(data_dir))
            for path in data_dir.glob(pattern)
            if path.is_file()
        ]
        if not matches:
            unmatched_patterns.append(pattern)
        included.update(matches)

    for pattern in include_small_patterns:
        matches = [
            normalize_relative(path.relative_to(data_dir))
            for path in data_dir.glob(pattern)
            if path.is_file() and path.stat().st_size <= include_small_max_bytes
        ]
        included.update(matches)

    files = [
        path
        for path in sorted(included)
        if not any(pattern_matches(path, pattern) for pattern in exclude_patterns)
    ]

    if not files:
        raise ValueError("Bundle include/exclude rules produced no files")
    if unmatched_patterns:
        print("Patterns with no local matches:")
        for pattern in unmatched_patterns:
            print(f"  {pattern}")
    return files


def write_file_list(files: list[str]) -> Path:
    handle = tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", delete=False, prefix="rolypoly-data-files-", suffix=".txt"
    )
    with handle:
        handle.write("\n".join(files))
        handle.write("\n")
    return Path(handle.name)


def create_archive(data_dir: Path, output: Path, files: list[str], threads: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    file_list = write_file_list(files)
    tar_cmd = [
        "tar",
        "--create",
        "--file",
        "-",
        "--directory",
        str(data_dir),
        "--files-from",
        str(file_list),
    ]

    try:
        if output.name.endswith((".tar.gz", ".tgz")):
            compressor = "pigz" if shutil.which("pigz") else "gzip"
            compressor_cmd = (
                [compressor, "-p", str(max(1, threads)), "-c"]
                if compressor == "pigz"
                else [compressor, "-c"]
            )
            with output.open("wb") as archive:
                tar_proc = subprocess.Popen(tar_cmd, stdout=subprocess.PIPE)
                if tar_proc.stdout is None:
                    raise RuntimeError("Could not open tar stdout")
                gzip_proc = subprocess.Popen(
                    compressor_cmd, stdin=tar_proc.stdout, stdout=archive
                )
                tar_proc.stdout.close()
                gzip_code = gzip_proc.wait()
                tar_code = tar_proc.wait()
            if tar_code or gzip_code:
                raise RuntimeError(
                    f"Archive creation failed: tar={tar_code}, compressor={gzip_code}"
                )
        else:
            with output.open("wb") as archive:
                result = subprocess.run(tar_cmd, stdout=archive, check=False)
            if result.returncode:
                raise RuntimeError(f"Archive creation failed: tar={result.returncode}")
    finally:
        file_list.unlink(missing_ok=True)


def write_sha256(archive: Path) -> Path:
    """Write a portable SHA-256 sidecar using the archive basename."""
    digest = hashlib.sha256()
    with archive.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    checksum = archive.with_name(f"{archive.name}.sha256")
    checksum.write_text(f"{digest.hexdigest()}  {archive.name}\n", encoding="utf-8")
    return checksum


def copy_bundle_files(data_dir: Path, stage_dir: Path, files: list[str]) -> None:
    """Copy the manifest-selected file set to a clean staging directory."""
    if stage_dir.exists() and any(stage_dir.iterdir()):
        raise FileExistsError(f"Staging directory is not empty: {stage_dir}")
    stage_dir.mkdir(parents=True, exist_ok=True)
    for relative_path in files:
        source = data_dir / relative_path
        destination = stage_dir / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


@command()
@option("--data-dir", required=True, type=Path, help="RolyPoly data directory.")
@option(
    "--manifest",
    "manifest_path",
    default=manifest_default_path,
    type=Path,
    show_default=True,
    help="Data contract manifest with bundle include/exclude rules.",
)
@option(
    "--output",
    type=Path,
    default=Path("data.tar.gz"),
    show_default=True,
    help="Output archive path.",
)
@option("--threads", default=4, show_default=True, help="pigz threads, if available.")
@option("--dry-run", is_flag=True, help="Print the archive file list and exit.")
@option(
    "--stage-dir",
    type=Path,
    help="Copy the selected bundle contents to this empty directory instead of archiving.",
)
def bundle_data(
    data_dir: Path,
    manifest_path: Path,
    output: Path,
    threads: int,
    dry_run,
    stage_dir: Path | None,
):
    """Create an archive or staging directory containing distributable runtime data."""
    data_dir = data_dir.resolve()
    manifest_path = manifest_path.resolve()
    output = output.resolve()

    files = collect_bundle_files(data_dir, manifest_path)
    total_size = sum((data_dir / path).stat().st_size for path in files)
    print(f"Bundle files: {len(files):,}")
    print(f"Uncompressed bytes: {total_size:,}")

    if dry_run:
        for path in files:
            print(path)
        return

    if stage_dir is not None:
        copy_bundle_files(data_dir, stage_dir.resolve(), files)
        print(f"Copied selected data files to {stage_dir.resolve()}")
        return

    create_archive(data_dir, output, files, threads)
    print(f"Wrote {output}")
    checksum = write_sha256(output)
    print(f"Wrote {checksum}")


if __name__ == "__main__":
    bundle_data()
