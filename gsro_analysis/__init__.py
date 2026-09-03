"""
gsro_analysis — analysis-side helpers for the global snowmelt runoff onset
dataset: paths (version-aware), settings (external-service IDs), datacube
(per-tile UTM ancillary construction), aggregate (the map/reduce engine:
partial sums -> cubes), stats, plotting, colorbars, world_maps, and
results (provenance-stamped tables).

Submodules import lazily (PEP 562): ``import gsro_analysis`` (or importing
``paths``/``plotting``) does not drag in the heavy stack (xdem,
earthengine-api, odc.stac, dask, …), and merely importing the package
never reconfigures the root logger. File logging into the repo's logs/ is
opt-in: call ``gsro_analysis.datacube.setup_logging()`` from entrypoints
that want it.
"""

import importlib

_SUBMODULES = ('aggregate', 'colorbars', 'datacube', 'era5', 'paths', 'plotting',
               'results', 'settings', 'stats', 'world_maps')


def __getattr__(name):
    if name in _SUBMODULES:
        return importlib.import_module(f'{__name__}.{name}')
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')


def __dir__():
    return sorted(set(globals()) | set(_SUBMODULES))
