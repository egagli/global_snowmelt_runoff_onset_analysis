"""Build analyses/continents/data/gtopo30_lat_elev_histogram.nc — the land-pixel count per continent x
1 degree latitude x 100 m elevation from GTOPO30 (USGS/GTOPO30 on Earth Engine, reduced at 1 km).

It is the grey background of the continental latitude-elevation panels (all land, not just the
mapped pixels), merged into the continents cube as `dem_pixel_count` by
analyses/continents/0_aggregate_by_continent.ipynb. Grid- and version-independent, so the 40 KB file
is TRACKED; rebuild it only if the latitude/elevation bins in gsro_analysis.aggregate change. The bins
are the cube's (aggregate.LAT_EDGES, aggregate.DEM_EDGES); the continents are the USGS polygons
(settings.CONTINENTS_URL, read from the web) with Australia folded into Oceania.

Needs the Earth Engine service key; a few minutes.

    pixi run gtopo30                        # = python pipeline/scripts/get_gtopo30_histogram.py
"""

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import xarray as xr

from gsro_analysis import aggregate, paths, settings

GTOPO30_ASSET = "USGS/GTOPO30"
SCALE_M = 1000


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', default=str(paths.gtopo30_histogram()))
    p.add_argument('--force', action='store_true', help='rebuild even if the output exists')
    p.add_argument('--ee-key', help='Earth Engine service-account key json (default: the production clone\'s config/ee_key.json)')
    return p.parse_args()


def continent_histogram(dem, continent_geometry, lat_bins, elevation_bins):
    """(latitude bin, elevation bin) land-pixel counts inside one continent polygon set."""
    import ee
    continent_fc = ee.FeatureCollection(continent_geometry.__geo_interface__)
    lat_binned = ee.Image.pixelLonLat().select('latitude').add(90).floor().int()      # 1 degree bins from -90
    dem_binned = dem.divide(100).floor().int()                                         # 100 m bins from 0
    binned = ee.Image.cat([ee.Image.constant(1), lat_binned, dem_binned]).rename(['count', 'lat_bin', 'dem_bin'])
    reducer = (ee.Reducer.sum()
               .group(groupField=1, groupName='lat_bin')
               .group(groupField=2, groupName='dem_bin'))
    result = binned.reduceRegion(reducer=reducer, geometry=continent_fc.geometry(),
                                 scale=SCALE_M, maxPixels=1e12, bestEffort=True).getInfo()
    counts = np.zeros((len(lat_bins), len(elevation_bins)))
    for dem_group in result.get('groups', []):
        dem_idx = int(dem_group.get('dem_bin', -1))
        if not 0 <= dem_idx < len(elevation_bins):
            continue
        for lat_group in dem_group.get('groups', []):
            lat_idx = int(lat_group.get('lat_bin', -1))
            if 0 <= lat_idx < len(lat_bins):
                counts[lat_idx, dem_idx] = lat_group.get('sum', 0)
    return counts


def main():
    args = parse_args()
    out = Path(args.out)
    if out.exists() and not args.force:
        print(f"{out} exists ({out.stat().st_size / 1e3:.0f} KB); --force rebuilds it")
        return
    import ee
    print(f"Earth Engine initialized as {settings.initialize_earthengine(key_file=args.ee_key)}", flush=True)

    continents_gdf = gpd.read_file('zip+' + settings.CONTINENTS_URL)
    dem = ee.Image(GTOPO30_ASSET)
    lat_bins = aggregate.bin_centers('latitude')
    elevation_bins = aggregate.bin_centers('elevation')
    counts = np.zeros((len(aggregate.CONTINENT_NAMES), len(lat_bins), len(elevation_bins)), dtype=np.int64)
    for i, continent in enumerate(aggregate.CONTINENT_NAMES):
        members = ['Oceania', 'Australia'] if continent == 'Oceania' else [continent]   # Australia -> Oceania
        geometry = continents_gdf[continents_gdf['CONTINENT'].isin(members)]
        print(f"GTOPO30 histogram: {continent} ({len(geometry)} polygons)", flush=True)
        counts[i] = continent_histogram(dem, geometry, lat_bins, elevation_bins)

    ds = xr.Dataset({'pixel_count': (('continent', 'latitude', 'elevation'), counts)},
                    coords={'continent': list(aggregate.CONTINENT_NAMES), 'latitude': lat_bins, 'elevation': elevation_bins})
    ds['latitude'].attrs['units'] = 'degrees'
    ds['elevation'].attrs['units'] = 'meters'
    ds['pixel_count'].attrs['description'] = f'Count of land pixels in each latitude-elevation bin ({GTOPO30_ASSET}, {SCALE_M} m scale)'
    out.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out, encoding={'pixel_count': {'zlib': True, 'complevel': 4}})
    print(f"wrote {out} ({out.stat().st_size / 1e3:.0f} KB), {int(counts.sum()):,} land pixels")


if __name__ == '__main__':
    main()
