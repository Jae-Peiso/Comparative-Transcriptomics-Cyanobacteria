from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any, Mapping, Sequence


def save_and_push(
    repo_path: str | Path | None = None,
    message: str | None = None,
    files: Sequence[str] | None = None,
    remote: str = "origin",
    branch: str | None = None,
    dry_run: bool = False,
    config: Mapping[str, Any] | None = None,
) -> dict[str, object]:
    """Stage files, commit them, and push the current branch.

    This helper accepts either explicit arguments or a compact config mapping,
    which makes it convenient to call from a notebook as:

    save_and_push(current_nb)

    where current_nb is a dict-like object with keys such as repo_path,
    message, files, remote, branch, and dry_run.
    """
    if config is not None:
        repo_path = config.get("repo_path", repo_path)
        message = config.get("message", message)
        files = config.get("files", files)
        remote = config.get("remote", remote)
        branch = config.get("branch", branch)
        dry_run = config.get("dry_run", dry_run)

    if repo_path is None:
        raise ValueError("repo_path is required")
    if message is None:
        raise ValueError("message is required")

    repo = Path(repo_path).resolve()
    if not (repo / ".git").exists():
        raise FileNotFoundError(f"{repo} is not a Git repository")

    if branch is None:
        branch = _get_current_branch(repo)

    if files is None:
        stage_args = ["git", "-C", str(repo), "add", "-A"]
    else:
        stage_args = ["git", "-C", str(repo), "add", *[str(Path(p)) for p in files]]

    subprocess.run(stage_args, check=True, capture_output=True, text=True)

    status = subprocess.run(
        ["git", "-C", str(repo), "status", "--porcelain"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()

    if not status:
        return {"committed": False, "pushed": False, "branch": branch, "message": message}

    if dry_run:
        return {"committed": False, "pushed": False, "branch": branch, "message": message, "dry_run": True}

    subprocess.run(
        ["git", "-C", str(repo), "commit", "-m", message],
        check=True,
        capture_output=True,
        text=True,
    )

    try:
        subprocess.run(
            ["git", "-C", str(repo), "push", remote, branch],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        # Commit already succeeded above -- surface git's own explanation (e.g. "Updates
        # were rejected because the remote contains work that you do not have locally")
        # instead of a bare traceback with no diagnostic info.
        raise RuntimeError(
            f"Committed locally, but push to {remote}/{branch} failed:\n{exc.stderr}"
        ) from exc

    return {"committed": True, "pushed": True, "branch": branch, "message": message}


def _get_current_branch(repo: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo), "branch", "--show-current"],
        check=True,
        capture_output=True,
        text=True,
    )
    branch = result.stdout.strip()
    if not branch:
        raise RuntimeError("Could not determine current branch")
    return branch
