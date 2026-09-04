"""Build data/global_hillshade_robinson.tif — the grey hillshade basemap under the two inset world
maps (gsro_analysis.world_maps), the river-basin maps (analyses/river_basins/basin_onset.ipynb) and the
Sierra Nevada case study (analyses/case_studies/sierra_nevada/sierra_nevada.ipynb).

Source: Natural Earth, "Gray Earth with Shaded Relief, Hypsography, Ocean Bottom, and Drainages"
(10 m raster GRAY_HR_SR_OB_DR, public domain; settings.HILLSHADE_URL), downloaded once into
data/sources/. Recipe: gdalwarp to World Robinson (ESRI:54030) at 1 km, average resampling, Byte,
Cloud-Optimized GeoTIFF with LZW — the recipe of the production repo's
visualize/data/download_and_preprocess_hillshade.ipynb, unchanged, so the product is the same file.
Value 0 = outside the ellipse (no data), 1-255 = grey.

No credentials. The zip is a few hundred MB and the warp takes a few minutes; the output (~240 MB) is
gitignored (*.tif). The script skips when the output exists (--force rebuilds).

    pixi run hillshade                       # = python pipeline/scripts/get_hillshade.py
    pixi run hillshade --out /tmp/test.tif   # build somewhere else
"""

import argparse
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

from gsro_analysis import paths, settings

GDALWARP_ARGS = ['-t_srs', 'ESRI:54030',             # World Robinson
                 '-te_srs', 'EPSG:4326', '-te', '-180', '-90', '180', '90',
                 '-tr', '1000', '1000',              # 1 km pixels
                 '-r', 'average',
                 '-ot', 'Byte',
                 '-of', 'COG', '-co', 'COMPRESS=LZW', '-co', 'BIGTIFF=IF_SAFER']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', default=str(paths.hillshade()), help='output GeoTIFF (default: data/global_hillshade_robinson.tif)')
    p.add_argument('--force', action='store_true', help='rebuild even if the output exists')
    return p.parse_args()


def main():
    args = parse_args()
    out = Path(args.out)
    if out.exists() and not args.force:
        print(f"{out} exists ({out.stat().st_size / 1e6:.0f} MB); --force rebuilds it")
        return
    gdalwarp = shutil.which('gdalwarp')
    if gdalwarp is None:
        sys.exit("gdalwarp not found on PATH: run through `pixi run` (libgdal is in the environment)")

    # 1. the Natural Earth zip, downloaded once
    zip_path = settings.cached_source(settings.HILLSHADE_URL, dest_dir=paths.DATA / 'sources')
    print(f"source zip: {zip_path} ({zip_path.stat().st_size / 1e6:.0f} MB)")

    # 2. the GeoTIFF inside it, extracted next to the zip
    extract_dir = zip_path.with_suffix('')
    with zipfile.ZipFile(zip_path) as zf:
        tif_member = next(n for n in zf.namelist() if n.lower().endswith('.tif'))
        src = extract_dir / tif_member
        if not src.exists():
            zf.extract(tif_member, extract_dir)
    print(f"source raster: {src} ({src.stat().st_size / 1e6:.0f} MB)")

    # 3. reproject to World Robinson at 1 km (the production recipe), atomically
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.stem + '.part.tif')
    cmd = [gdalwarp, '-overwrite', *GDALWARP_ARGS, str(src), str(tmp)]
    print('$', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    tmp.replace(out)
    print(f"wrote {out} ({out.stat().st_size / 1e6:.0f} MB)")


if __name__ == '__main__':
    main()
