"""Provenance-stamped results tables: the production package's
``save_result_table`` plus this repository's commit.

Upstream ``global_snowmelt_runoff_onset.results.save_result_table`` appends
``_version``, ``_git_sha`` and ``_written_at`` — and its ``_git_sha`` is the
commit of the PRODUCTION package (it resolves ``git rev-parse`` in its own
directory), not of the code that computed the numbers. Every table written
from here therefore also carries ``_analysis_git_sha``, the short SHA of
this repository's HEAD. The upstream column is left as it is (it records
which production code the pipeline ran against); nothing upstream changes.
"""

import subprocess

from gsro_analysis import paths


def analysis_git_sha():
    """Short SHA of this repository's HEAD (None outside a git checkout)."""
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=5,
                             cwd=paths.ROOT)
        return out.stdout.strip() or None
    except Exception:  # noqa: BLE001
        return None


def save_result_table(df, name, results_dir, version=None):
    """Write ``df`` as ``<results_dir>[/<version>]/<name>.csv`` with both
    provenance stamps. ``results_dir`` is normally ``paths.resultsdir(topic,
    config.version)``, which already carries the version segment — then
    leave ``version`` None."""
    from global_snowmelt_runoff_onset.results import save_result_table as upstream

    stamped = df.copy()
    stamped["_analysis_git_sha"] = analysis_git_sha()
    return upstream(stamped, name, version=version, results_dir=results_dir)
