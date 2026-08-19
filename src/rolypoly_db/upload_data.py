"""Upload a bundled RolyPoly data archive to Zenodo or NERSC."""

import json
import os
import subprocess
from pathlib import Path
from typing import Any

import requests
from rich_click import argument, group, option


ZENODO_API_URL = "https://zenodo.org/api"
ZENODO_SANDBOX_API_URL = "https://sandbox.zenodo.org/api"
ZENODO_DEPOSITION_ID = "21639934"
DATA_ARCHIVE_NAME = "data.tar.gz"
DATA_CHECKSUM_NAME = "data.tar.gz.sha256"


def require_existing_file(path: Path) -> Path:
    """Return an existing file path or raise a clear error."""
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"File does not exist: {path}")
    return path


def require_release_name(path: Path, expected_name: str) -> Path:
    """Require the stable public filename used by RolyPoly get-data."""
    path = require_existing_file(path)
    if path.name != expected_name:
        raise ValueError(f"Expected {expected_name}, got {path.name}")
    return path


def zenodo_headers(token: str) -> dict[str, str]:
    """Return Zenodo authorization headers."""
    return {"Authorization": f"Bearer {token}"}


def request_json(
    method: str,
    url: str,
    token: str,
    *,
    json_payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Send one Zenodo API request and return JSON."""
    response = requests.request(
        method,
        url,
        headers=zenodo_headers(token),
        json=json_payload,
        timeout=120,
    )
    response.raise_for_status()
    return response.json()


def load_metadata(metadata: Path | None, title: str | None) -> dict[str, Any]:
    """Load Zenodo metadata JSON or create a minimal draft metadata object."""
    if metadata is not None:
        return json.loads(metadata.read_text(encoding="utf-8"))
    if not title:
        raise ValueError("Either --metadata or --title is required for Zenodo upload")
    return {
        "metadata": {
            "title": title,
            "upload_type": "dataset",
            "description": "RolyPoly external runtime data bundle.",
            "creators": [{"name": "Neri, Uri"}],
            "access_right": "open",
        }
    }


def create_or_update_deposition(
    api_url: str,
    token: str,
    metadata_payload: dict[str, Any] | None,
    deposition_id: str | None,
    new_version: bool,
) -> dict[str, Any]:
    """Create, update, or fork a Zenodo draft deposition."""
    if deposition_id and new_version:
        draft = request_json(
            "POST",
            f"{api_url}/deposit/depositions/{deposition_id}/actions/newversion",
            token,
        )
        latest_draft = draft["links"]["latest_draft"]
        deposition = request_json("GET", latest_draft, token)
    elif deposition_id:
        deposition = request_json(
            "GET", f"{api_url}/deposit/depositions/{deposition_id}", token
        )
    else:
        deposition = request_json(
            "POST", f"{api_url}/deposit/depositions", token, json_payload={}
        )

    if metadata_payload is None:
        return deposition

    return request_json(
        "PUT",
        f"{api_url}/deposit/depositions/{deposition['id']}",
        token,
        json_payload=metadata_payload,
    )


def upload_file_to_bucket(
    bucket_url: str, token: str, file_path: Path, remote_name: str
) -> dict[str, Any]:
    """Upload one file to a Zenodo bucket."""
    with file_path.open("rb") as handle:
        response = requests.put(
            f"{bucket_url}/{remote_name}",
            data=handle,
            headers=zenodo_headers(token),
            timeout=None,
        )
    response.raise_for_status()
    return response.json()


@group()
def upload_data():
    """Upload bundled RolyPoly data archives."""


@upload_data.command("zenodo")
@option("--archive", required=True, type=Path, help="Data archive to upload.")
@option("--checksum", type=Path, help="Optional checksum file to upload too.")
@option("--metadata", type=Path, help="Zenodo deposition metadata JSON.")
@option("--title", help="Minimal metadata title when --metadata is omitted.")
@option(
    "--deposition-id",
    default=lambda: os.environ.get("ZENODO_DEPOSITION_ID", ZENODO_DEPOSITION_ID),
    show_default=ZENODO_DEPOSITION_ID,
    help="Existing latest Zenodo deposition ID.",
)
@option(
    "--token",
    default=lambda: os.environ.get("ZENODO_ACCESS_TOKEN"),
    help="Zenodo token. Defaults to ZENODO_ACCESS_TOKEN.",
)
@option("--new-version", is_flag=True, help="Create a new version from deposition ID.")
@option("--publish", is_flag=True, help="Publish the draft after upload.")
@option("--sandbox", is_flag=True, help="Use sandbox.zenodo.org.")
@option("--dry-run", is_flag=True, help="Print intended actions without uploading.")
def upload_zenodo(
    archive: Path,
    checksum: Path | None,
    metadata: Path | None,
    title: str | None,
    deposition_id: str | None,
    token: str | None,
    new_version: bool,
    publish: bool,
    sandbox: bool,
    dry_run: bool,
):
    """Upload an archive to Zenodo using environment-held credentials."""
    archive = require_release_name(archive, DATA_ARCHIVE_NAME)
    checksum = (
        require_release_name(checksum, DATA_CHECKSUM_NAME) if checksum else None
    )
    metadata_payload = (
        load_metadata(metadata, title) if metadata is not None or title else None
    )
    api_url = ZENODO_SANDBOX_API_URL if sandbox else ZENODO_API_URL

    planned_files = [archive.name]
    if checksum:
        planned_files.append(checksum.name)
    print(f"Zenodo API: {api_url}")
    print(f"Deposition ID: {deposition_id or '<new draft>'}")
    print(f"Files: {', '.join(planned_files)}")
    print(f"Update metadata: {metadata_payload is not None}")
    print(f"Publish: {publish}")
    if dry_run:
        return
    if not token:
        raise ValueError("Set ZENODO_ACCESS_TOKEN or pass --token")
    if new_version and not deposition_id:
        raise ValueError("--new-version requires --deposition-id or ZENODO_DEPOSITION_ID")

    deposition = create_or_update_deposition(
        api_url, token, metadata_payload, deposition_id, new_version
    )
    bucket_url = deposition["links"]["bucket"]
    upload_file_to_bucket(bucket_url, token, archive, archive.name)
    if checksum:
        upload_file_to_bucket(bucket_url, token, checksum, checksum.name)

    if publish:
        request_json(
            "POST",
            f"{api_url}/deposit/depositions/{deposition['id']}/actions/publish",
            token,
        )
    print(f"Zenodo deposition ID: {deposition['id']}")


@upload_data.command("zenodo-drafts")
@option(
    "--token",
    default=lambda: os.environ.get("ZENODO_ACCESS_TOKEN"),
    help="Zenodo token. Defaults to ZENODO_ACCESS_TOKEN.",
)
@option("--sandbox", is_flag=True, help="Use sandbox.zenodo.org.")
def list_zenodo_drafts(token: str | None, sandbox: bool):
    """List Zenodo draft depositions without uploading or publishing anything."""
    if not token:
        raise ValueError("Set ZENODO_ACCESS_TOKEN or pass --token")
    api_url = ZENODO_SANDBOX_API_URL if sandbox else ZENODO_API_URL
    response = requests.get(
        f"{api_url}/deposit/depositions",
        headers=zenodo_headers(token),
        params={"status": "draft", "sort": "mostrecent", "size": 100},
        timeout=120,
    )
    response.raise_for_status()
    drafts = response.json()
    if not drafts:
        print("No draft depositions found.")
        return
    for draft in drafts:
        metadata = draft.get("metadata", {})
        title = metadata.get("title", "<untitled>")
        modified = draft.get("modified", "<unknown modified date>")
        print(f"{draft['id']}\t{modified}\t{title}")


@upload_data.command("nersc")
@argument("archive", type=Path)
@option(
    "--remote-dir",
    default=lambda: os.environ.get("NERSC_REMOTE_DIR"),
    help="Remote directory. Defaults to NERSC_REMOTE_DIR.",
)
@option(
    "--user",
    default=lambda: os.environ.get("NERSC_USER"),
    help="NERSC username. Defaults to NERSC_USER.",
)
@option(
    "--host",
    default=lambda: os.environ.get("NERSC_HOST", "dtn01.nersc.gov"),
    show_default=True,
    help="NERSC DTN host.",
)
@option("--checksum", type=Path, help="Optional checksum file to upload too.")
@option("--dry-run", is_flag=True, help="Print the scp command without running it.")
def upload_nersc(
    archive: Path,
    remote_dir: str | None,
    user: str | None,
    host: str,
    checksum: Path | None,
    dry_run: bool,
):
    """Upload an archive to a NERSC DTN with scp."""
    archive = require_release_name(archive, DATA_ARCHIVE_NAME)
    checksum = (
        require_release_name(checksum, DATA_CHECKSUM_NAME) if checksum else None
    )
    if not remote_dir:
        raise ValueError("Set NERSC_REMOTE_DIR or pass --remote-dir")
    if not user:
        raise ValueError("Set NERSC_USER or pass --user")

    destination = f"{user}@{host}:{remote_dir.rstrip('/')}/"
    cmd = ["scp", str(archive)]
    if checksum:
        cmd.append(str(checksum))
    cmd.append(destination)
    print(" ".join(cmd))
    if dry_run:
        return
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    upload_data()
