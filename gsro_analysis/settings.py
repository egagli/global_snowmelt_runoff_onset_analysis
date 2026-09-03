"""Single source of the external-service identifiers and data-source URLs
that used to be hardcoded across the notebooks.

No cloud-cluster settings exist on purpose: batch compute is local
execution for dry runs and aggregation, and the GitHub Actions fleet pattern
(the production repo's icechunk fleet, adapted) for full tile campaigns.
"""

# Earth Engine (needed only for CHILI and the GTOPO30 histogram).
# The service-account key lives in the production clone:
# ../global_snowmelt_runoff_onset/config/ee_key.json
EE_PROJECT = "egagli-data-access"   # for INTERACTIVE auth only; the service
                                    # account carries its own GCP project
EE_HIGHVOLUME_URL = "https://earthengine-highvolume.googleapis.com"


def zarr_store_exists(fs, path):
    """Skip-if-exists check that works for both zarr formats: v3 stores mark
    their root with zarr.json, v2-consolidated with .zmetadata."""
    return fs.exists(f"{path}/zarr.json") or fs.exists(f"{path}/.zmetadata")


def ee_key_path():
    """The EE service-account key in the side-by-side production clone."""
    from gsro_analysis import paths
    return paths.ROOT.parent / "global_snowmelt_runoff_onset" / "config" / "ee_key.json"


def initialize_earthengine(key_file=None):
    """Initialize Earth Engine, service-account first.

    Verified working HEADLESS (2026-08-24): the key's client_email is read
    from the json itself and the service account carries its own GCP
    project, so no ``project=`` is passed — forcing :data:`EE_PROJECT`
    would 403 a service account that isn't registered in it. Falls back to
    interactive credentials (+ EE_PROJECT) only when no key file is found.
    """
    import json

    import ee

    key_file = key_file or ee_key_path()
    try:
        with open(key_file, "rb") as f:
            raw = f.read()
        email = json.loads(raw)["client_email"]
    except (FileNotFoundError, KeyError):
        ee.Initialize(project=EE_PROJECT, opt_url=EE_HIGHVOLUME_URL)
        return "interactive"
    credentials = ee.ServiceAccountCredentials(email, str(key_file))
    ee.Initialize(credentials, opt_url=EE_HIGHVOLUME_URL)
    # easysnowdata's @requires_earthengine checks EARTHENGINE_TOKEN or the
    # interactive credentials FILE — never the live session. On a CI runner
    # neither exists (fleet smoke run 2026-08-25 failed in get_chili with a
    # fully initialized session), so export the token in the form easysnowdata
    # documents (base64 of the service-account json).
    import base64
    import os
    os.environ.setdefault("EARTHENGINE_TOKEN",
                          base64.b64encode(raw).decode())
    return email


def fresh_blob_fs(config):
    """A NEW Azure filesystem instance, bypassing fsspec's instance cache.

    After minutes of Earth-Engine/odc-heavy work the process-shared
    ``config.azure_blob_fs`` can list EMPTY (zarr group reads come back
    with 0 data vars while a fresh process sees everything; teardown logs
    'Loop is not running'). Seen 2026-08-25 in the ancillary worker and in
    ERA5 Acquire run 32885451802 (a complete store failed its own post-write
    check). Anything that must be TRUSTED — verification reads, opens inside
    long-lived workers — goes through a fresh instance."""
    import adlfs
    return adlfs.AzureBlobFileSystem(
        account_name=config.azure_storage_account,
        sas_token=config.sas_token,
        skip_instance_cache=True)


def verify_in_subprocess(config, module, func, *args, timeout=900):
    """Run ``gsro_analysis.<module>.<func>(Config(<config file>), *args)`` in
    a FRESH interpreter and return its truthiness.

    Post-write verification reads through fsspec, whose event loop can be
    dead after 30+ minutes of Earth-Engine/dask work ('Loop is not running'
    at teardown): three complete stores failed their own in-process check
    on 2026-08-25 (a year store, then the anomaly store twice) while a fresh
    process read them fine. A fresh interpreter has a fresh loop, so the
    builders call this as the fallback before declaring a write failed.
    ``args`` must be JSON-serializable."""
    import json
    import subprocess
    import sys
    code = ("import json, sys\n"
            "from global_snowmelt_runoff_onset.config import Config\n"
            f"from gsro_analysis import {module}\n"
            f"ok = {module}.{func}(Config(sys.argv[1]), *json.loads(sys.argv[2]))\n"
            "sys.exit(0 if ok else 1)")
    result = subprocess.run([sys.executable, '-c', code,
                             str(config.config_file_path), json.dumps(list(args))],
                            capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0 and result.stderr:
        print(result.stderr.strip().splitlines()[-1][:300])
    return result.returncode == 0


# The ONE place the dataset version is chosen: every notebook and script
# builds its Config through load_config(). Override per shell with
# GSRO_CONFIG=config/global_config_v11.txt. Code assumes >= v10 (icechunk,
# versioned Azure layouts) - there is no v9 compatibility.
CONFIG_FILE = "config/global_config_v10.txt"


def load_config(config_file=None):
    """Production ``Config`` for the current dataset version
    (``config_file`` > ``$GSRO_CONFIG`` > :data:`CONFIG_FILE`); the path is
    resolved inside the side-by-side production clone."""
    import os
    from global_snowmelt_runoff_onset.config import Config
    return Config(config_file or os.environ.get("GSRO_CONFIG", CONFIG_FILE))


# Azure prefixes for analysis artifacts (container-relative, account from
# the production Config). <version> is config.version; the ancillary store
# is keyed by GRID generation, not dataset version — it only changes when
# the grid does.
ANCILLARY_PREFIX = "snowmelt/analysis/ancillary"                # + /<grid>/tile_RRR_CCC.zarr
PARTIALS_PREFIX = "snowmelt/analysis/partials"                  # + /<version>/tile_RRR_CCC.parquet  (the fleet product)
ANALYSIS_PARQUET_PREFIX = "snowmelt/analysis/parquets"          # + /<version>/tile_RRR_CCC.parquet  (opt-in pixel tables)
ERA5_DATA_PREFIX = "snowmelt/analysis/era5_data"                # + /<version>
AGGREGATED_PREFIX = "snowmelt/analysis/aggregated"              # + /<version>/  (mirror of aggregated_results/<version>)


# v9 per-tile parquets: kept ONLY as the validation reference for the
# tile-for-tile comparison in pipeline/pipeline.ipynb (v9 tile (r, c) == v10 tile (r+2, c))
V9_TILE_PARQUET_PREFIX = "snowmelt/analysis/parquets/tiles/v9"

# external vector/raster sources
GMBA_URL = ("https://data.earthenv.org/mountains/standard/"
            "GMBA_Inventory_v2.0_standard_300.zip")
CONTINENTS_URL = ("https://pubs.usgs.gov/of/2006/1187/basemaps/continents/"
                  "continents.zip")
# BasinATLAS v1.0 (HydroATLAS, figshare article 9890531, file 20082137).
# The direct figshare.com/ndownloader URL sits behind an AWS WAF bot
# challenge (returns 202 to non-browser clients) — the API endpoint
# redirects cleanly to a presigned S3 URL instead.
BASIN_ATLAS_URL = "https://api.figshare.com/v2/file/download/20082137"
BASIN_ATLAS_MD5 = "69af94baee68da5a3f80f09e7b85bd04"   # from the figshare API
# The basin id layer the fleet stores: HydroBASINS LEVEL 6 (decided
# 2026-09-02; ~16,400 basins, six-digit PFAF_ID, int32). Pfafstetter codes
# are prefix-hierarchical, so every coarser level is derived at reduce time
# (level j = PFAF_ID // 10**(6 - j)); the default `river_basins` cube is level
# 5 and the analyses, era5_zonal.py and the population join use level-5
# polygons (basin_atlas_layer(5)). Changing this after a campaign means a
# fleet pass (datacube.refresh_unit_layers + re-map).
BASIN_ATLAS_STORED_LEVEL = 6
BASIN_ATLAS_LAYER = "BasinATLAS_v10_lev06"


def basin_atlas_layer(level):
    """gdb layer name of HydroBASINS level ``level`` (1-12) in the cached
    BasinATLAS v1.0 file."""
    return f"BasinATLAS_v10_lev{int(level):02d}"


def basin_atlas_gdb():
    """The cached BasinATLAS gdb zip (downloaded once, md5-verified)."""
    return cached_source(BASIN_ATLAS_URL, filename='BasinATLAS_Data_v10.gdb.zip',
                        expected_md5=BASIN_ATLAS_MD5)


def cached_source(url, filename=None, max_retries=6, backoff=20,
                  expected_md5=None):
    """Download ``url`` once into ``data/geometries/sources/`` and return
    the local Path.

    The per-tile ancillary build must NEVER hit these vector sources over
    the network per tile: BasinATLAS alone is ~2.5 GB, and the direct
    figshare download URL sits behind an AWS WAF bot challenge (HTTP 202
    to non-browser clients — use the api.figshare.com endpoint instead).
    Retries transient statuses (202/429/5xx) with backoff; anything else
    raises. ``expected_md5`` verifies the streamed bytes.
    """
    import hashlib
    import time
    import urllib.error
    import urllib.request

    from gsro_analysis import paths

    dest_dir = paths.GEOMETRIES / "sources"
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / (filename or url.rstrip("/").split("/")[-1])
    if dest.exists() and dest.stat().st_size > 0:
        return dest

    tmp = dest.with_name(dest.name + ".part")
    last_status = None
    for attempt in range(1, max_retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "gsro-analysis"})
            with urllib.request.urlopen(req, timeout=300) as r:
                if r.status == 200:
                    md5 = hashlib.md5()
                    with open(tmp, "wb") as f:
                        while chunk := r.read(1 << 22):
                            f.write(chunk)
                            md5.update(chunk)
                    if expected_md5 and md5.hexdigest() != expected_md5:
                        tmp.unlink()
                        raise RuntimeError(
                            f"md5 mismatch for {url}: got {md5.hexdigest()}, "
                            f"expected {expected_md5}")
                    tmp.rename(dest)
                    return dest
                last_status = r.status  # e.g. a WAF challenge 202
        except urllib.error.HTTPError as e:
            last_status = e.code
            if e.code not in (202, 408, 429, 500, 502, 503, 504):
                raise
        except (urllib.error.URLError, TimeoutError) as e:
            last_status = str(e)
        print(f"{url} -> {last_status}; retry {attempt}/{max_retries} in {backoff}s")
        time.sleep(backoff)
    raise RuntimeError(f"failed to download {url} after {max_retries} attempts "
                       f"(last status: {last_status})")
# account is uwcryo (the old snowmelt.blob... hostname is gone — NXDOMAIN);
# the blob is anonymously readable, no SAS needed
FOREST_COVER_FRACTION_URL = (
    "https://uwcryo.blob.core.windows.net/snowmelt/eric/"
    "forest_cover_fraction/PROBAV_LC100_global_v3.0.1_2019-nrt_"
    "Tree-CoverFraction-layer_EPSG-4326.tif")
# Sturm & Liston (2021) seasonal snow classification, 300 m (NSIDC-0768),
# public mirror — the ERA5 zonal join's seasonal-snow weight (class 4 =
# ephemeral is excluded, matching the pixel filters' snow_classification != 4)
SNOW_CLASS_URL = (
    "https://uwcryo.blob.core.windows.net/snowmelt/eric/"
    "snow_classification/SnowClass_GL_300m_10.0arcsec_2021_v01.0.tif")
