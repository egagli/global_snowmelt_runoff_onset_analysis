"""Repo-root-anchored, version-aware paths for the notebooks and scripts.

Layout (since 2026-09-04)::

    analyses/<unit>/                           one folder per aggregation unit: continents, mountain_ranges, river_basins
        0_aggregate_by_<unit>.ipynb            the fleet's partial sums -> the unit's cube(s), zonal means, metrics table
        data/geometries/                       the unit's polygons, downloaded once (gitignored except small tracked tables)
        data/aggregation/<version>/            the cubes and zonal means (gitignored, regenerable)
        figures/<version>/  results/<version>/ tracked per dataset version
    partials/<version>/tile_RRR_CCC.parquet    local cache of the fleet product (gitignored)
    data/                                      shared inputs: the hillshade basemap and its source zip (gitignored)

Import from here instead of writing cwd-relative string literals, so a notebook resolves the same
paths no matter where it is run from, and always says which dataset version it is working against
(``config.version``, e.g. ``'v10'``)::

    from gsro_analysis import paths

    mountain_ranges_ds = xr.open_dataset(paths.aggregation_dir('mountain_ranges', config.version)
                                         / 'all_mountain_ranges_fcf_lte_50.nc')
    metrics_df = pd.read_csv(paths.resultsdir('mountain_ranges', config.version) / 'mountain_range_metrics.csv')
    f.savefig(paths.figdir('mountain_ranges', config.version) / 'global_mad_by_mountain_range.png', dpi=300)

Two environment overrides exist for running against a test set without touching the tree (the CI
smoke test): ``GSRO_PARTIALS_ROOT`` (the partials cache, default ``partials/``) and
``GSRO_OUTPUT_ROOT`` (everything a notebook writes: figures, results and ``data/aggregation``, which keep
the ``analyses/<unit>/...`` layout under it; default ``analyses/``). There are deliberately no default
versions or filter tags: every artifact names the dataset version it came from and the pixel filter it
was built with (``gsro_analysis.aggregate.FILTERS``).
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DATA = ROOT / "data"                # shared inputs (hillshade + its source zip); everything else lives next to its analysis
ANALYSES = ROOT / "analyses"
PIPELINE = ROOT / "pipeline"
SCRATCH = ROOT / "scratch"
LOGS = ROOT / "logs"

PARTIALS = Path(os.environ.get("GSRO_PARTIALS_ROOT", ROOT / "partials"))
OUTPUT_ROOT = Path(os.environ.get("GSRO_OUTPUT_ROOT", ANALYSES))

UNITS = ("continents", "mountain_ranges", "river_basins")   # the aggregation units = the analyses/ folders with a 0_aggregate notebook


def geometries(unit):
    """The unit's polygon folder ``analyses/<unit>/data/geometries/`` (created if missing): the
    downloaded GMBA / BasinATLAS / continents sources and, for river basins, the tracked
    level-6 population table."""
    path = ANALYSES / unit / "data" / "geometries"
    path.mkdir(parents=True, exist_ok=True)
    return path


def aggregation_dir(unit, version):
    """``analyses/<unit>/data/aggregation/<version>/`` (created if missing): the unit's cubes
    ``all_<unit>_<filter>.nc`` and its ERA5-Land zonal means ``era5_anomaly_<unit>.nc``, written by
    ``0_aggregate_by_<unit>.ipynb``."""
    path = OUTPUT_ROOT / unit / "data" / "aggregation" / version
    path.mkdir(parents=True, exist_ok=True)
    return path


def figdir(topic, version, *subdirs):
    """Version-scoped figure directory of an analysis folder, created if missing. ``topic`` may be
    nested, e.g. ``figdir('case_studies/sierra_nevada', 'v10')``.

    >>> figdir('mountain_ranges', 'v10', 'triplets')
    PosixPath('.../analyses/mountain_ranges/figures/v10/triplets')
    """
    path = OUTPUT_ROOT / topic / "figures" / version
    if subdirs:
        path = path.joinpath(*subdirs)
    path.mkdir(parents=True, exist_ok=True)
    return path


def resultsdir(topic, version):
    """Version-scoped results directory ``analyses/<topic>/results/<version>/`` (created if
    missing): the provenance-stamped tables, e.g. ``mountain_range_metrics.csv``."""
    path = OUTPUT_ROOT / topic / "results" / version
    path.mkdir(parents=True, exist_ok=True)
    return path


def partials_cache(version):
    """Local cache of the fleet's per-tile partial-sum parquets: ``partials/<version>/`` (created
    if missing), filled by the aggregation notebooks from Azure."""
    path = PARTIALS / version
    path.mkdir(parents=True, exist_ok=True)
    return path


def gtopo30_histogram():
    """The static GTOPO30 latitude x elevation land-pixel histogram (tracked, 40 KB, grid- and
    version-independent): ``analyses/continents/data/gtopo30_lat_elev_histogram.nc``, rebuilt by
    ``pipeline/scripts/get_gtopo30_histogram.py`` (``pixi run gtopo30``, Earth Engine)."""
    return ANALYSES / "continents" / "data" / "gtopo30_lat_elev_histogram.nc"


def hillshade():
    """The global hillshade basemap ``data/global_hillshade_robinson.tif`` (Natural Earth, World
    Robinson, 1 km; gitignored), built by ``pipeline/scripts/get_hillshade.py`` (``pixi run hillshade``)."""
    return DATA / "global_hillshade_robinson.tif"


def logfile(name="analysis.log"):
    """Path inside the repo-root ``logs/`` directory, created if missing."""
    LOGS.mkdir(parents=True, exist_ok=True)
    return LOGS / name
