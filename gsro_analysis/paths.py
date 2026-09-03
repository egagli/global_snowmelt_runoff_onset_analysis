"""Repo-root-anchored, version-aware paths for the analysis notebooks.

Notebooks live in per-topic folders under ``analyses/`` but the shared
inputs stay under ``data/`` and pipeline outputs under
``aggregated_results/<version>/``. Import from here instead of writing
cwd-relative string literals, so a notebook resolves the same paths no
matter where it is run from — and always says which dataset version it
is working against (``config.version``, e.g. ``'v10'`` / ``'v11'``)::

    from gsro_analysis import paths

    ds = xr.open_dataset(paths.aggregate('mountain_ranges', config.version))
    metrics = pd.read_csv(paths.resultsdir('mountain_ranges', config.version)
                          / 'mountain_range_metrics.csv')
    f.savefig(paths.figdir('mountain_ranges', config.version) / 'global_mad.png',
              dpi=300)

There are deliberately no default versions or filter tags: every artifact
names the dataset version it came from and the pixel filter it was built
with (see ``gsro_analysis.aggregate.FILTERS``).
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# shared inputs (all under data/ — gitignored, see data/README.md)
DATA = ROOT / "data"
GEOMETRIES = DATA / "geometries"
ERA5 = DATA / "era5"

# outputs and working dirs. Two environment overrides exist for running the
# notebooks against a test set or a different disk without touching the
# tracked tree: GSRO_AGGREGATED_ROOT (pipeline outputs, default
# aggregated_results/) and GSRO_OUTPUT_ROOT (figures/results, default
# analyses/ — figdir/resultsdir keep the analyses/<topic>/... layout under it).
AGGREGATED = Path(os.environ.get("GSRO_AGGREGATED_ROOT", ROOT / "aggregated_results"))
ANALYSES = ROOT / "analyses"
OUTPUT_ROOT = Path(os.environ.get("GSRO_OUTPUT_ROOT", ANALYSES))
PIPELINE = ROOT / "pipeline"
SCRATCH = ROOT / "scratch"
LOGS = ROOT / "logs"


def figdir(topic, version, *subdirs):
    """Version-scoped figure directory of an analysis folder, created if
    missing. ``topic`` may be nested, e.g.
    ``figdir('case_studies/sierra_nevada', 'v10')``.

    >>> figdir('mountain_ranges', 'v10', 'triplets', 'pngs')
    PosixPath('.../analyses/mountain_ranges/figures/v10/triplets/pngs')
    """
    path = OUTPUT_ROOT / topic / "figures" / version
    if subdirs:
        path = path.joinpath(*subdirs)
    path.mkdir(parents=True, exist_ok=True)
    return path


def resultsdir(topic, version):
    """Version-scoped results directory for
    ``gsro_analysis.results.save_result_table`` — pass it as ``results_dir``
    (the version segment is this path's job, so leave ``version=None``).
    The per-range metrics table lives at
    ``resultsdir('mountain_ranges', version) / 'mountain_range_metrics.csv'``
    (``stats.range_metrics_gdf`` joins it to the GMBA polygons on read)."""
    path = OUTPUT_ROOT / topic / "results" / version
    path.mkdir(parents=True, exist_ok=True)
    return path


def aggregate(group, version, filter_tag="fcf_lte_50"):
    """Combined aggregate output for a whole group:
    ``aggregated_results/<version>/<group>/all_<group>_<filter_tag>.nc``.
    ``group`` is ``mountain_ranges``, ``river_basins`` or ``continents``."""
    return AGGREGATED / version / group / f"all_{group}_{filter_tag}.nc"


def aggregate_dir(version, group=None):
    """Directory holding a version's aggregate outputs (optionally one
    group's), created if missing — where the pipeline writes per-unit files."""
    path = AGGREGATED / version if group is None else AGGREGATED / version / group
    path.mkdir(parents=True, exist_ok=True)
    return path


def partials_cache(version):
    """Local cache of the fleet's per-tile partial-sum parquets for a
    version (downloaded by pipeline/scripts/reduce_partials.py)."""
    path = AGGREGATED / version / "partials"
    path.mkdir(parents=True, exist_ok=True)
    return path


def era5_zonal(version, unit_type):
    """Per-unit ERA5-Land anomaly zonal means (pipeline/scripts/era5_zonal.py):
    ``aggregated_results/<version>/era5_zonal/era5_anomaly_<unit_type>.nc``."""
    path = AGGREGATED / version / "era5_zonal" / f"era5_anomaly_{unit_type}.nc"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def gtopo30_histogram():
    """Static GTOPO30 latitude x elevation land-pixel histogram (input,
    grid/version independent): ``data/gtopo30_lat_elev_histogram.nc``."""
    return DATA / "gtopo30_lat_elev_histogram.nc"


def logfile(name="analysis.log"):
    """Path inside the repo-root ``logs/`` directory, created if missing."""
    LOGS.mkdir(parents=True, exist_ok=True)
    return LOGS / name
