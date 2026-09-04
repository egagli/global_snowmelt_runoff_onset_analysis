"""Provenance stamps for the results tables the notebooks write.

Every table under ``analyses/<unit>/results/<version>/`` carries four underscore-prefixed columns
(underscore so they never collide with a data column), appended right before the explicit
``to_csv`` in the notebook::

    metrics_df.assign(**results.provenance(config)).to_csv(path, index=False)

``_version``          the dataset version (``config.version``)
``_git_sha``          short SHA of the production package (the side-by-side ``global_snowmelt_runoff_onset`` clone)
``_analysis_git_sha`` short SHA of this repository
``_written_at``       UTC time, ISO 8601
"""

import subprocess
from datetime import datetime, timezone

from gsro_analysis import paths

PRODUCTION_REPO = paths.ROOT.parent / "global_snowmelt_runoff_onset"


def git_short_sha(repo_dir):
    """Short SHA of ``repo_dir``'s HEAD (None outside a git checkout)."""
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True,
                             timeout=5, cwd=repo_dir)
        return out.stdout.strip() or None
    except Exception:  # noqa: BLE001
        return None


def analysis_git_sha():
    """Short SHA of this repository's HEAD."""
    return git_short_sha(paths.ROOT)


def provenance(config):
    """The four provenance columns as a dict, for ``df.assign(**provenance(config))``."""
    return {"_version": config.version,
            "_git_sha": git_short_sha(PRODUCTION_REPO),
            "_analysis_git_sha": analysis_git_sha(),
            "_written_at": datetime.now(timezone.utc).isoformat(timespec="seconds")}
