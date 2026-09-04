"""The two inset world maps, drawn natively in matplotlib.

    Figure A  spring-temperature sensitivity choropleth + per-range scatterplot insets   plot_temperature_sensitivity_map
    Figure B  elevation lapse-rate choropleth + per-range polar-triplet insets          plot_lapse_rate_triplet_map

Each map is assembled by its own notebook, which also builds the component parts and documents
every knob (label placement, adding or hiding ranges, text style, page, colorbar and legend
positions): ``analyses/mountain_ranges/topography_triplet_composite_figure.ipynb`` (B) and
``spring_temperature_sensitivity_composite_figure.ipynb`` (A).

Inputs (all pipeline products or tracked files; nothing figure-specific except the label table):

* ``stats_gdf``: GMBA polygons + per-range metrics, one row per range, any CRS (the GMBA
  inventory merged on ``GMBA_V2_ID`` with ``results/<version>/mountain_range_metrics.csv``, the
  table ``0_aggregate_by_mountain_range.ipynb`` writes). Columns used: ``GMBA_V2_ID``, ``MapName``,
  ``Level_04``, ``anomaly_slope``, ``anomaly_corr``, ``snowmelt_lapse_rate_per_100m``,
  ``snowmelt_lapse_rate_n``.
* ``analyses/mountain_ranges/label_layout.csv``: the curated display flags and the hand-placed
  label anchors in World Robinson metres (``display_map``, ``display_label``,
  ``show_topo``/``show_anom``, ``topo_x/y``, ``anom_x/y``). An anchor is the BOTTOM-LEFT corner
  of the label block. The only input that is neither data nor code.
* the per-range sweep PNGs (``figures/<version>/{triplets,anomaly_scatterplots}/``), the
  triplet legend PNG and the hillshade GeoTIFF in ``data/`` (``pipeline/scripts/get_hillshade.py``). The 60° x 30° graticule is
  generated here.

Page model: a figure of the page size (mm); a map axes in projected metres over the map item;
one page-sized overlay axes in page millimetres (origin top-left) that carries callouts, inset
images, label text, legend and frames; colorbar axes at the picture positions.

Provenance. The published figures were laid out in a QGIS project; this module reproduces its
pages from the same inputs (validated side by side, 2026-09-01) and then departs from it where
legibility asked for it. The project, its spec and the transfer material are archived in the
private repo ``recreate_global_snowmelt_runoff_onset_analysis_QGIS_figures_in_mpl``. Numbers
that came out of that reverse engineering and are NOT obvious from the code:

* Inset width: QGIS rendered the HTML ``<img width=120*90e6/@map_scale px>`` at 1 px = 1 pt,
  so the insets are 46.5 mm (A) and 53.8 mm (B) wide (``INSET_WIDTH_PX``).
* The image sits 3 pt above the anchor (``HTML_MARGIN_PT``); the name is one text line above
  the image, centred on it.
* Figure A names are bold and coloured with the feature's own fill colour; the ColorBrewer
  ramps are continuous; hillshade 0 = outside the ellipse = no-data; the graticule is 0.45 pt.
* The ``*_QGIS`` styles and ``top_margin_mm=None`` reproduce the project's look; the default
  styles are larger and bolder, the page top is cut (A) or extended (B) to ``TOP_MARGIN_MM``
  above the highest label, and B's colorbar sits 11.5 mm further left than in the project so
  its label clears the legend.

:func:`layout_report` lists collisions (name vs name, name vs inset, colorbar vs legend, off
page) after any change; the ``plot_*`` functions warn when it finds one.

Usage::

    from gsro_analysis import world_maps
    metrics_df = pd.read_csv(paths.resultsdir('mountain_ranges', config.version) / 'mountain_range_metrics.csv')
    gmba_gdf = gpd.read_file('zip+' + settings.GMBA_URL)
    stats_gdf = gmba_gdf[['GMBA_V2_ID', 'MapName', 'Level_04', 'geometry']].merge(metrics_df.drop(columns=['name']), on='GMBA_V2_ID')
    fig, placed = world_maps.plot_lapse_rate_triplet_map(
        stats_gdf, config.version,
        out=paths.figdir('mountain_ranges', config.version) / 'global_lapse_rates_with_triplets_map.png')
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, replace
from pathlib import Path

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from matplotlib import patheffects
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.patches import Rectangle
from PIL import Image
from rasterio.enums import Resampling
from shapely.geometry import Point, box

from gsro_analysis import colorbars, paths

CRS = 'ESRI:54030'                       # World Robinson, +proj=robin +datum=WGS84 +units=m
MM_PER_PT = 25.4 / 72
LAYOUT_CSV = paths.ANALYSES / 'mountain_ranges' / 'label_layout.csv'
HILLSHADE = paths.hillshade()               # built by pipeline/scripts/get_hillshade.py
HILLSHADE_STRETCH = (1, 231)             # single-band grey, linear 1 -> 231 (black -> white)
GRATICULE_DEG = (60, 30)                 # graticule cell size (lon, lat) in degrees
FONT_FAMILIES = ('Arial', 'Liberation Sans', 'DejaVu Sans')   # first one found is used
TEXT_COLOR = '#323232'

# QGIS HTML label: <div style="margin-bottom: 3px">MapName</div><br/><img width=120*90e6/@map_scale>
INSET_WIDTH_PX = 120 * 90e6              # / map scale -> CSS px, rendered by QGIS at 1 px = 1 pt
HTML_MARGIN_PT = 3                       # gap text -> image, and image -> anchor (measured 1.06 mm)
FONT_ASCENT, FONT_DESCENT = 0.905, 0.212  # Arial, in em (used for the label-block rectangle)

# z-order inside the page overlay axes
Z_CALLOUT, Z_INSET, Z_TEXT, Z_LEGEND, Z_FRAME = 1, 2, 3, 4, 5


# --------------------------------------------------------------------------- pages
@dataclass(frozen=True)
class Page:
    """A QGIS layout page: size in mm, the map item (full width, from the top), its extent."""
    width: float
    height: float
    map_height: float
    extent: tuple                       # (x0, x1, y0, y1) in Robinson metres

    @property
    def scale(self):
        """Map scale denominator (1 : scale)."""
        x0, x1, _, _ = self.extent
        return (x1 - x0) / (self.width / 1000)

    @property
    def m_per_mm(self):
        """Map metres per page millimetre (≈ 82 km on A, 71 km on B)."""
        return self.scale / 1000

    def to_page(self, x, y):
        """Robinson metres -> page millimetres (origin top-left, y down)."""
        x0, x1, y0, y1 = self.extent
        return ((np.asarray(x) - x0) / (x1 - x0) * self.width,
                (y1 - np.asarray(y)) / (y1 - y0) * self.map_height)

    def from_page(self, x_mm, y_mm):
        """Page millimetres -> Robinson metres: the inverse of :meth:`to_page`, for editing
        anchors in ``label_layout.csv`` ("move it 10 mm right" = ``+10 * page.m_per_mm``)."""
        x0, x1, y0, y1 = self.extent
        return (x0 + np.asarray(x_mm) / self.width * (x1 - x0),
                y1 - np.asarray(y_mm) / self.map_height * (y1 - y0))

    def trim_top(self, trim_mm):
        """The same page with ``trim_mm`` removed from the top (negative = added): the map
        extent, page and map-item heights change together, so nothing on the map moves
        relative to the map."""
        x0, x1, y0, y1 = self.extent
        return Page(self.width, self.height - trim_mm, self.map_height - trim_mm,
                    (x0, x1, y0, y1 - trim_mm * self.m_per_mm))

    def mm_to_fig(self, x, y, w, h):
        """Page-mm rectangle (top-left origin) -> matplotlib figure-fraction rect."""
        return [x / self.width, 1 - (y + h) / self.height, w / self.width, h / self.height]


PAGES = {
    'temperature_sensitivity': Page(480, 300, 300,
                                    (-19467691.66, 19904028.44, -10706262.74, 13901062.32)),
    'lapse_rates_with_triplets': Page(480, 313.254, 289.201,
                                      (-16933759.31, 17077012.83, -9555248.82, 10936277.14)),
}

# layout pictures (spec §1.6; page mm, top-left origin): the item rect, its frame width (None =
# no frame) and whether the embedded artwork was opaque (figure A's colorbar PNG has a white
# background; figure B's colorbar SVG rendered transparent, the legend PNG carries its own white)
PICTURES = {
    'temperature_sensitivity': dict(
        colorbar=dict(rect=(194.126, 246.279, 226.386, 50.335), frame_mm=1.0, background='white')),
    'lapse_rates_with_triplets': dict(
        # x was 119.456 in the project, where the bar's label ran under the legend frame;
        # moved 11.5 mm left (2026-09-01) so the label clears the legend by ~4 mm
        colorbar=dict(rect=(108.0, 276.201, 200.173, 37.052), frame_mm=None, background=None),
        legend=dict(rect=(312.320, 238.064, 166.575, 74.267), frame_mm=1.0, background=None)),
}


# --------------------------------------------------------------------------- data
def inset_stem(map_name, level_04):
    """Inset filename stem, the project's rule: Level_04 instead of MapName for the three
    Andes cordilleras; spaces and hyphens -> ``_``; parentheses removed."""
    name = level_04 if any(k in map_name for k in ('Cordillera Occidental', 'Cordillera Central',
                                                   'Cordillera Oriental')) else map_name
    return name.replace(' ', '_').replace('-', '_').replace('(', '').replace(')', '')


def prepare(stats_gdf, layout_csv=LAYOUT_CSV, simplify_m=1000):
    """Join the per-range stats with the label table, derive the display flags the way the
    archived QGIS project did, project to Robinson and simplify (1000 m, invisible at
    ~7 km per output pixel; the project's layer files used the same tolerance)."""
    layout = pd.read_csv(layout_csv)
    gdf = stats_gdf.merge(layout.drop(columns=['MapName', 'note'], errors='ignore'),
                          on='GMBA_V2_ID', how='left')
    gdf['display_map'] = np.where(gdf['display_map'] == 0, 0, 1)        # curated: 0 for four ranges
    gdf['display_label'] = gdf['display_label'].fillna(0).astype(int)   # curated: the 41 ranges
    for col in ('show_topo', 'show_anom'):                              # the .qgd show overrides
        gdf[col] = gdf[col].fillna(1).astype(int)
    gdf = gdf.to_crs(CRS)
    if simplify_m:
        gdf['geometry'] = gdf.geometry.simplify(simplify_m, preserve_topology=True)
    return gdf


def sensitivity_fill(gdf):
    """Figure A fill: YlOrRd at clip(-anomaly_slope/12, 0, 1); alpha interp(anomaly_corr,
    [-1, 0], [1, 76/255]) (clamped). Returns (rgba (n, 4), drawn mask)."""
    drawn = (gdf['display_map'] == 1) & gdf['anomaly_slope'].notna()
    t = np.clip(-gdf['anomaly_slope'].to_numpy(dtype=float) / 12, 0, 1)
    rgba = mpl.colormaps['YlOrRd'](np.nan_to_num(t))
    rgba[:, 3] = np.interp(np.nan_to_num(gdf['anomaly_corr'].to_numpy(dtype=float)),
                           [-1, 0], [1, 76 / 255])
    rgba[~drawn.to_numpy()] = 0
    return rgba, drawn.to_numpy()


def lapse_rate_fill(gdf):
    """Figure B fill: YlGnBu at clip(rate/8, 0, 1), opaque, drawn when n >= 10."""
    drawn = (gdf['display_map'] == 1) & (gdf['snowmelt_lapse_rate_n'].fillna(0) >= 10)
    t = np.clip(gdf['snowmelt_lapse_rate_per_100m'].to_numpy(dtype=float) / 8, 0, 1)
    rgba = mpl.colormaps['YlGnBu'](np.nan_to_num(t))
    rgba[:, 3] = 1
    rgba[~drawn.to_numpy()] = 0
    return rgba, drawn.to_numpy()


# --------------------------------------------------------------------------- base map
def _read_hillshade(path, decimation):
    """Decimated hillshade (average resampling) with 0 = outside the ellipse masked; the mask
    comes from a nearest-neighbour read so averaging does not leave a dark rim."""
    if not Path(path).exists():
        raise FileNotFoundError(f"{path}: no hillshade basemap; build it with `pixi run hillshade` "
                                "(pipeline/scripts/get_hillshade.py, no credentials)")
    with rasterio.open(path) as src:
        out_shape = (src.height // decimation, src.width // decimation)
        values = src.read(1, out_shape=out_shape, resampling=Resampling.average)
        nodata = src.read(1, out_shape=out_shape, resampling=Resampling.nearest) == 0
        b = src.bounds
    return np.ma.masked_array(values, mask=nodata), (b.left, b.right, b.bottom, b.top)


def graticule(cell_deg=GRATICULE_DEG, densify_deg=0.5):
    """The graticule as projected cell boundaries: a grid of ``cell_deg`` (lon, lat) cells over
    the globe (60 x 30 -> 36 cells), each boundary densified before projection so the curves
    are smooth. Cells rather than lines so shared edges are drawn twice, which is how the
    original polygon-grid layer looked at 50 % opacity."""
    lon_step, lat_step = cell_deg
    cells = [box(lon, lat, lon + lon_step, lat + lat_step)
             for lon in range(-180, 180, lon_step) for lat in range(-90, 90, lat_step)]
    return gpd.GeoSeries(cells, crs='EPSG:4326').boundary.segmentize(densify_deg).to_crs(CRS)


def _base_map(fig, page, hillshade=HILLSHADE, graticule_deg=GRATICULE_DEG, hillshade_decimation=4):
    """Map axes over the page's map item: hillshade (grey 1-231) and the graticule (white,
    50 %, 0.45 pt dashed). Returns the axes (projected metres)."""
    ax = fig.add_axes(page.mm_to_fig(0, 0, page.width, page.map_height), zorder=0)
    x0, x1, y0, y1 = page.extent
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect('auto')     # the extent has the map item's aspect to < 0.01 %, so 1 m = 1 m
    ax.axis('off')
    values, extent = _read_hillshade(hillshade, hillshade_decimation)
    ax.imshow(values, extent=extent, cmap='gray', vmin=HILLSHADE_STRETCH[0],
              vmax=HILLSHADE_STRETCH[1], interpolation='antialiased', aspect='auto', zorder=1)
    if graticule_deg:
        graticule(graticule_deg).plot(ax=ax, color='white', alpha=0.5, linewidth=0.45,
                                      linestyle=(0, (4, 2)), zorder=2)
    return ax


def _choropleth(ax, gdf, rgba, drawn):
    gdf[drawn].plot(ax=ax, color=list(rgba[drawn]), edgecolor='none', linewidth=0, zorder=3)


# --------------------------------------------------------------------------- labels
_FONT = None


def font_family():
    """First of FONT_FAMILIES matplotlib can find (Arial is what the QGIS project used)."""
    global _FONT
    if _FONT is None:
        for family in FONT_FAMILIES:
            try:
                findfont(FontProperties(family=family), fallback_to_default=False)
                _FONT = family
                break
            except ValueError:
                continue
        else:
            _FONT = 'sans-serif'
    return _FONT


@dataclass(frozen=True)
class LabelStyle:
    """How a range name and its callout are drawn. ``STYLE_A_QGIS`` / ``STYLE_B_QGIS`` reproduce
    the project; ``STYLE_A`` / ``STYLE_B`` (the defaults) are the more legible variants: larger
    text, no alpha on the name, a light halo outside the buffer. Pass a
    ``dataclasses.replace(STYLE_B, size_pt=13)`` as ``style=`` to tune."""
    size_pt: float
    weight: str
    buffer_mm: float
    buffer_color: str
    shadow: bool
    callout_lw_mm: float
    callout_origin: str          # 'centroid' (label block centre) or 'exterior' (nearest boundary point)
    callout_alpha: float         # symbol opacity (multiplies the colour's alpha)
    colored_text: bool           # text colour = feature fill colour (A) or TEXT_COLOR (B)
    colored_callout: bool
    gap_pt: float = HTML_MARGIN_PT   # text baseline block -> image top (QGIS: the 3 px div margin)
    opaque_text: bool = False        # drop the fill colour's alpha for the name (legibility)
    halo_mm: float = 0.0             # extra light halo outside the buffer (0 = none)
    halo_color: str = 'white'
    inset_width_mm: float | None = None   # None = the project's rule, 120*90e6/scale pt (46.5 / 53.8 mm)

    def text_color(self, fill):
        if not self.colored_text:
            return TEXT_COLOR
        return (*fill[:3], 1.0) if self.opaque_text else tuple(fill)

    def path_effects(self):
        buffer_lw = 2 * self.buffer_mm / MM_PER_PT       # stroke is centred on the outline
        effects = []
        if self.shadow:
            # QGIS: black, 70 %, offset 1.0 mm at 135 deg (lower right), blur 1.5 mm; the blur is
            # approximated by four stroked shadows of growing width and falling opacity
            offset = (0.7071 / MM_PER_PT, -0.7071 / MM_PER_PT)
            for extra_mm, alpha in ((0.0, 0.2), (0.5, 0.15), (1.0, 0.11), (1.5, 0.08)):
                effects.append(patheffects.SimplePatchShadow(
                    offset=offset, shadow_rgbFace='black', alpha=alpha,
                    foreground='black', linewidth=buffer_lw + extra_mm / MM_PER_PT))
        if self.halo_mm:
            effects.append(patheffects.Stroke(linewidth=buffer_lw + 2 * self.halo_mm / MM_PER_PT,
                                              foreground=self.halo_color))
        effects.append(patheffects.Stroke(linewidth=buffer_lw, foreground=self.buffer_color))
        effects.append(patheffects.Normal())
        return effects


def _load_inset(path, width_px):
    """Inset PNG resized (Lanczos, alpha-aware) to its on-page pixel width at the output dpi."""
    img = Image.open(path).convert('RGBA')
    height_px = max(1, round(width_px * img.height / img.width))
    return np.asarray(img.resize((round(width_px), height_px), Image.LANCZOS))


_Y_COLS = ('anchor_y_mm', 'inset_top', 'inset_bottom', 'baseline', 'block_top', 'block_bottom')


def _label_geometry(page, labels, rgba, style, inset_dir, suffix, keep_on_page=False):
    """The label block of every labelled range, in page mm, without drawing anything: anchor
    (bottom-left of the block), inset rect (its bottom 3 pt above the anchor), text baseline,
    block rect, polygon centroid, fill colour. ``keep_on_page`` shifts a block down when it
    would poke above the page (the fallback when the page is not resized to the labels)."""
    width_mm = style.inset_width_mm or INSET_WIDTH_PX / page.scale * MM_PER_PT
    margin_mm = HTML_MARGIN_PT * MM_PER_PT
    size_mm = style.size_pt * MM_PER_PT
    rows = []
    for i, row in zip(np.flatnonzero(labels['_labelled'].to_numpy()), labels[labels['_labelled']].itertuples()):
        ax_mm, ay_mm = page.to_page(row.anchor_x, row.anchor_y)
        png = Path(inset_dir) / f'{inset_stem(row.MapName, row.Level_04)}{suffix}.png'
        if png.exists():
            with Image.open(png) as im:
                aspect = im.height / im.width
        else:
            warnings.warn(f'{row.MapName}: no inset {png.name} in {inset_dir}; drawing text and callout only')
            aspect = 0.35
        height_mm = width_mm * aspect
        bottom = ay_mm - margin_mm
        top = bottom - height_mm
        baseline = top - style.gap_pt * MM_PER_PT - FONT_DESCENT * size_mm
        block_top = baseline - FONT_ASCENT * size_mm - style.buffer_mm - style.halo_mm
        shift = max(0.0, -block_top) if keep_on_page else 0.0
        cx, cy = page.to_page(row.geometry.centroid.x, row.geometry.centroid.y)
        fill = tuple(float(v) for v in rgba[i])
        rows.append(dict(GMBA_V2_ID=row.GMBA_V2_ID, MapName=row.MapName, png=str(png), has_png=png.exists(),
                         anchor_x_mm=ax_mm, anchor_y_mm=ay_mm + shift, shift_down_mm=shift,
                         inset_left=ax_mm, inset_right=ax_mm + width_mm, inset_top=top + shift,
                         inset_bottom=bottom + shift, inset_w_mm=width_mm, inset_h_mm=height_mm,
                         baseline=baseline + shift, block_left=ax_mm, block_right=ax_mm + width_mm,
                         block_top=block_top + shift, block_bottom=ay_mm + shift,
                         centroid_x=float(cx), centroid_y=float(cy),
                         fill_r=fill[0], fill_g=fill[1], fill_b=fill[2], fill_a=fill[3]))
    return pd.DataFrame(rows)


def _draw_labels(page_ax, geo, style, dpi):
    """Draw the callouts, inset images and names described by :func:`_label_geometry`; adds
    the callout origin to the table (page mm)."""
    font = font_family()
    origins = []
    for row in geo.itertuples():
        block = box(row.block_left, row.block_top, row.block_right, row.block_bottom)
        target = Point(row.centroid_x, row.centroid_y)
        if style.callout_origin == 'exterior':                  # QGIS "point on exterior" = shortest line
            origin = block.exterior.interpolate(block.exterior.project(target))
        else:
            origin = block.centroid
        fill = (row.fill_r, row.fill_g, row.fill_b, row.fill_a)
        if style.colored_callout:
            callout_color = (*fill[:3], fill[3] * style.callout_alpha)
        else:
            callout_color = (0, 0, 0, style.callout_alpha)
        page_ax.plot([origin.x, row.centroid_x], [origin.y, row.centroid_y], color=callout_color,
                     linewidth=style.callout_lw_mm / MM_PER_PT, solid_capstyle='round',
                     zorder=Z_CALLOUT, clip_on=False)
        if row.has_png:
            img = _load_inset(row.png, row.inset_w_mm / 25.4 * dpi)
            page_ax.imshow(img, extent=(row.inset_left, row.inset_right, row.inset_bottom, row.inset_top),
                           aspect='auto', interpolation='antialiased', zorder=Z_INSET)
        page_ax.text((row.inset_left + row.inset_right) / 2, row.baseline, row.MapName, ha='center',
                     va='baseline', fontsize=style.size_pt, fontfamily=font, fontweight=style.weight,
                     color=style.text_color(fill), path_effects=style.path_effects(),
                     zorder=Z_TEXT, clip_on=False)
        origins.append((origin.x, origin.y))
    geo = geo.copy()
    geo[['callout_x0', 'callout_y0']] = origins
    geo['callout_x1'] = geo['centroid_x']
    geo['callout_y1'] = geo['centroid_y']
    return geo


def layout_report(fig, placed):
    """Collisions on a drawn map, from the RENDERED text extents: name vs name, name vs another
    range's inset, inset vs inset, anything off the page, and the colorbar's full extent
    (ticks + label) against the other pictures and the insets. ``placed`` is the table a
    ``plot_*`` function returned (its ``.attrs`` carry the page, pictures and style)."""
    page, pictures, style = (placed.attrs[k] for k in ('page', 'pictures', 'style'))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    page_ax = fig.axes[0]
    to_mm = page_ax.transData.inverted()

    def mm_box(bb, pad=0.0):
        (x0, y0), (x1, y1) = to_mm.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
        return box(min(x0, x1) - pad, min(y0, y1) - pad, max(x0, x1) + pad, max(y0, y1) + pad)

    pad = style.buffer_mm     # the visible edge of a name; the lighter halo outside it is not counted
    texts = {t.get_text(): mm_box(t.get_window_extent(renderer), pad) for t in page_ax.texts}
    insets = {r.MapName: box(r.inset_left, r.inset_top, r.inset_right, r.inset_bottom) for r in placed.itertuples()}
    page_box = box(0, 0, page.width, page.height)
    report = dict(
        text_text=[(a, b) for a in texts for b in texts if a < b and texts[a].intersects(texts[b])],
        text_inset=[(a, b) for a in texts for b in insets if a != b and texts[a].intersects(insets[b])],
        inset_inset=[(a, b) for a in insets for b in insets if a < b and insets[a].intersects(insets[b])],
        off_page=sorted({n for n, g in list(texts.items()) + list(insets.items()) if not page_box.contains(g)}),
        # None when fewer than two blocks are drawn (a partially processed dataset version)
        min_text_inset_clearance_mm=(round(min(d), 2) if (d := [texts[a].distance(insets[b]) for a in texts for b in insets if a != b]) else None),
    )
    cax = placed.attrs.get('colorbar_ax')
    if cax is not None:
        cbb = mm_box(cax.get_tightbbox(renderer))
        report['colorbar_bbox_mm'] = tuple(round(v, 1) for v in cbb.bounds)
        report['colorbar_overlaps'] = ([n for n, p in pictures.items() if n != 'colorbar' and cbb.intersects(box(
            p['rect'][0], p['rect'][1], p['rect'][0] + p['rect'][2], p['rect'][1] + p['rect'][3]))]
            + [n for n, g in insets.items() if cbb.intersects(g)])
        if not page_box.contains(cbb):
            report['off_page'].append('colorbar')
    return report


# --------------------------------------------------------------------------- layout pictures
def _fit_top_left(x, y, w, h, aspect):
    """QGIS picture item, 'zoom' mode, top-left anchor: the largest (w', h') of the given
    aspect (w/h) inside the frame; returns (x, y, w', h')."""
    if w / h > aspect:
        w = h * aspect
    else:
        h = w / aspect
    return x, y, w, h


def _frame(page_ax, picture, fitted):
    """The picture item's frame (1 mm black, on the item rect) and, for opaque artwork, the
    white ground of the fitted image area."""
    x, y, w, h = picture['rect']
    if picture.get('background'):
        page_ax.add_patch(Rectangle(fitted[:2], fitted[2], fitted[3], facecolor=picture['background'],
                                    edgecolor='none', zorder=Z_LEGEND - 0.5))
    if picture.get('frame_mm'):
        page_ax.add_patch(Rectangle((x, y), w, h, fill=False, edgecolor='black',
                                    linewidth=picture['frame_mm'] / MM_PER_PT, zorder=Z_FRAME, clip_on=False))


def _colorbar(fig, page_ax, page, preset, picture):
    """Draw a ``gsro_analysis.colorbars`` preset natively into the page at the position the
    QGIS layout gave its picture item. The layout embedded the preset's standalone export
    (8 x 3 in figure, tight bbox, 0.1 in pad) fitted into the frame with its aspect kept, so
    the standalone geometry is measured once and everything (bar axes, fonts) is scaled by the
    same factor. The bar axes is taken at its ORIGINAL position: the builder's ``extend``
    triangles shrink whatever axes they are drawn into, in the standalone and here alike."""
    standalone = preset()
    standalone.canvas.draw()
    tight = standalone.get_tightbbox(standalone.canvas.get_renderer()).padded(0.1)   # inches
    fw, fh = standalone.get_size_inches()
    b = standalone.axes[0].get_position(original=True)                                # figure fraction
    bar = (b.x0 * fw, b.y1 * fh, b.width * fw, b.height * fh)                         # inches: x0, y1 (top), w, h
    plt.close(standalone)
    fitted = _fit_top_left(*picture['rect'], tight.width / tight.height)
    x, y, w, h = fitted
    s = w / tight.width                                     # page mm per standalone inch
    f = s / 25.4                                            # font scale (25.4 mm per standalone inch)
    rect = (x + (bar[0] - tight.x0) * s, y + (tight.y1 - bar[1]) * s, bar[2] * s, bar[3] * s)
    _frame(page_ax, picture, fitted)
    cax = fig.add_axes(page.mm_to_fig(*rect), zorder=2)
    preset(ax=cax, label_fontsize=15 * f, tick_labelsize=15 * f, text_fontsize=18 * f)
    return cax


def _legend_picture(page_ax, path, picture):
    img = Image.open(path).convert('RGBA')
    fitted = _fit_top_left(*picture['rect'], img.width / img.height)
    x, y, w, h = fitted
    page_ax.imshow(np.asarray(img), extent=(x, x + w, y + h, y), aspect='auto',
                   interpolation='antialiased', zorder=Z_LEGEND)
    _frame(page_ax, picture, fitted)


# --------------------------------------------------------------------------- the two figures
def _new_page(page, dpi):
    fig = plt.figure(figsize=(page.width / 25.4, page.height / 25.4), dpi=100, facecolor='white')
    page_ax = fig.add_axes([0, 0, 1, 1], zorder=1)         # page overlay, millimetres, y down
    page_ax.set_xlim(0, page.width)
    page_ax.set_ylim(page.height, 0)
    page_ax.axis('off')
    page_ax.patch.set_visible(False)
    return fig, page_ax


def _draw(kind, stats_gdf, version, style, fill_fn, show_col, anchor_cols, inset_sub, suffix,
          layout_csv, hillshade, graticule_deg, hillshade_decimation, out, dpi, top_margin_mm,
          pictures=None, legend_png=None):
    page = PAGES[kind]
    pictures = {k: dict(v) for k, v in (pictures or PICTURES[kind]).items()}
    gdf = prepare(stats_gdf, layout_csv)
    rgba, drawn = fill_fn(gdf)
    labels = gdf.copy()
    labels['anchor_x'] = labels[anchor_cols[0]]
    labels['anchor_y'] = labels[anchor_cols[1]]
    labels['_labelled'] = ((labels['display_map'] == 1) & (labels['display_label'] == 1)
                           & (labels[show_col] == 1) & labels['anchor_x'].notna())
    inset_dir = paths.figdir('mountain_ranges', version, inset_sub)

    # page height: either the project's page (blocks that would poke above it are shifted
    # down) or resized so the highest label block sits top_margin_mm below the top edge (the
    # pictures move with the map; nothing else changes)
    geo = _label_geometry(page, labels, rgba, style, inset_dir, suffix, keep_on_page=top_margin_mm is None)
    if top_margin_mm is not None and len(geo):
        trim = float(geo['block_top'].min() - top_margin_mm)
        page = page.trim_top(trim)
        for pic in pictures.values():
            x, y, w, h = pic['rect']
            pic['rect'] = (x, y - trim, w, h)
        geo = _label_geometry(page, labels, rgba, style, inset_dir, suffix)

    fig, page_ax = _new_page(page, dpi)
    ax = _base_map(fig, page, hillshade, graticule_deg, hillshade_decimation)
    _choropleth(ax, gdf, rgba, drawn)
    placed = _draw_labels(page_ax, geo, style, dpi)

    preset = {'temperature_sensitivity': colorbars.temperature_sensitivity,
              'lapse_rates_with_triplets': colorbars.elevation_delay}[kind]
    cax = _colorbar(fig, page_ax, page, preset, pictures['colorbar'])
    if 'legend' in pictures:
        legend_png = legend_png or paths.figdir('mountain_ranges', version) / 'polar_triplet_legend.png'
        if Path(legend_png).exists():
            _legend_picture(page_ax, legend_png, pictures['legend'])
        else:
            warnings.warn(f'no triplet legend at {legend_png}; run the legend cell of topography_triplet_composite_figure.ipynb')

    placed.attrs.update(kind=kind, page=page, pictures=pictures, style=style, colorbar_ax=cax)
    report = layout_report(fig, placed)
    placed.attrs['layout_report'] = report
    problems = {k: v for k, v in report.items() if k in ('text_text', 'text_inset', 'inset_inset', 'off_page',
                                                          'colorbar_overlaps') and v}
    if problems:
        warnings.warn(f'{kind}: layout collisions {problems} (see world_maps.layout_report)')
    if out is not None:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=dpi, facecolor='white')
    return fig, placed


# The project's label styles. Figure A's font size is 300 000 map metres (RenderMetersInMapUnits):
# 3.66 mm = 10.4 pt on the page; its name colour is the feature fill INCLUDING its alpha.
STYLE_A_QGIS = LabelStyle(size_pt=300000 / PAGES['temperature_sensitivity'].scale * 1000 / MM_PER_PT,
                          weight='bold', buffer_mm=0.5, buffer_color='black', shadow=False,
                          callout_lw_mm=0.8, callout_origin='centroid', callout_alpha=0.502,
                          colored_text=True, colored_callout=True)
STYLE_B_QGIS = LabelStyle(size_pt=10, weight='normal', buffer_mm=1.0, buffer_color='#FAFAFA', shadow=True,
                          callout_lw_mm=0.4, callout_origin='exterior', callout_alpha=1.0,
                          colored_text=False, colored_callout=False)

# Legible defaults (2026-09-01): as large as the curated anchors allow without a name touching a
# neighbouring inset (measured with Arial metrics: 11 pt on A, where 12 pt collides in the
# Sierra Nevada row; 12 pt on B with the text gap halved). A: opaque name colours, thin black
# buffer, white halo so dark reds read over the dark ocean. B: bold, thin white buffer, no
# shadow (the project's 1 mm buffer + blurred shadow read as bubbles at 12 pt).
STYLE_A = replace(STYLE_A_QGIS, size_pt=11, opaque_text=True, buffer_mm=0.4, halo_mm=0.5)
STYLE_B = replace(STYLE_B_QGIS, size_pt=12, gap_pt=1.5, weight='bold', buffer_mm=0.6,
                  buffer_color='white', shadow=False)
TOP_MARGIN_MM = 3.0   # page top sits this far above the highest label block (None = project page)


def plot_temperature_sensitivity_map(stats_gdf, version, layout_csv=LAYOUT_CSV, hillshade=HILLSHADE,
                                     graticule_deg=GRATICULE_DEG, hillshade_decimation=4, out=None, dpi=300,
                                     style=None, top_margin_mm=TOP_MARGIN_MM, pictures=None):
    """Figure A: spring-temperature sensitivity choropleth (YlOrRd, alpha from
    the correlation) with a per-range scatterplot inset for the 30 curated ranges. Page 480 mm
    wide; 300 mm tall with ``top_margin_mm=None``, otherwise cut to the labels (≈ 287 mm).
    ``style``: a :class:`LabelStyle` (default ``STYLE_A``; ``STYLE_A_QGIS`` for the project's
    look); ``pictures``: a ``PICTURES``-shaped dict to move the colorbar. Returns
    ``(fig, placed)``; ``placed`` lists every inset drawn (page mm) and carries the layout
    report in ``placed.attrs``. The composite notebook documents every knob."""
    return _draw('temperature_sensitivity', stats_gdf, version, style or STYLE_A, sensitivity_fill,
                 'show_anom', ('anom_x', 'anom_y'), 'anomaly_scatterplots',
                 '_temperature_2m_spring_months_mean', layout_csv, hillshade, graticule_deg,
                 hillshade_decimation, out, dpi, top_margin_mm, pictures)


def plot_lapse_rate_triplet_map(stats_gdf, version, layout_csv=LAYOUT_CSV, hillshade=HILLSHADE,
                                graticule_deg=GRATICULE_DEG, hillshade_decimation=4, out=None, dpi=300,
                                legend_png=None, style=None, top_margin_mm=TOP_MARGIN_MM, pictures=None):
    """Figure B: elevation lapse-rate choropleth (YlGnBu 0-8 d/100 m) with a
    polar-triplet inset for the 40 curated ranges, the elevation-delay colorbar and the
    triplet legend. Page 480 mm wide; 313.254 mm tall with ``top_margin_mm=None``, otherwise
    sized to the labels (≈ 317 mm). ``style`` (default ``STYLE_B``; ``STYLE_B_QGIS`` for the
    project's look), ``pictures`` (colorbar and legend rects). Returns ``(fig, placed)``.
    The composite notebook documents every knob."""
    return _draw('lapse_rates_with_triplets', stats_gdf, version, style or STYLE_B, lapse_rate_fill,
                 'show_topo', ('topo_x', 'topo_y'), 'triplets', '', layout_csv, hillshade, graticule_deg,
                 hillshade_decimation, out, dpi, top_margin_mm, pictures, legend_png=legend_png)
