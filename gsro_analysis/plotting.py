"""Shared plotting helpers for the mountain-range analyses.

Extracted verbatim from ``mountain_ranges/geo_and_topo_analysis.ipynb``, where
each of these was defined inline — several of them two to five times over, in
successive iterations of the same figure. ``era5_analysis.ipynb`` carried its
own copies of ``count_valid_obs`` and ``wrap_labels``.

Two functions previously read notebook globals; those are now module-level
constants below, keeping the original lowercase names so the function bodies
are unchanged from the notebook versions:

* ``create_mini_polar_plot`` and ``plot_mountain_range_anomalies`` use
  ``major_tick_radii``, ``rgrid_vals``, ``rgrid_labels_blank``, ``label_angle``
  (defined identically in cells 21 and 27 of the source notebook).
* ``get_sorted_ranges`` closed over the notebook's ``mountains_ds`` while
  taking an unused ``df`` argument; it now uses the argument. Every call site
  passed ``mountains_ds`` as ``df``, so behaviour is unchanged.

Usage::

    from gsro_analysis import plotting

    lapse_rates_df = stats.lapse_rate_weighted_bins(
        mountains_ds, range_metadata)
"""

import textwrap

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ---------------------------------------------------------------------------
# polar-plot radial grid (source notebook cells 21 and 27)
# ---------------------------------------------------------------------------
major_tick_radii = [0,1000,2000,3000,4000,5000,6000,7000,8000,9000,]
rgrid_vals = [0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,]
rgrid_labels_blank = ['','','','','','','','','','','','','','','','','','','',]
rgrid_labels = ['0m','','1000m','','2000m','','3000m','','4000m','','5000m','','6000m','','7000m','','8000m','','9000m',]
label_angle = 225


# ---------------------------------------------------------------------------
# global triplet-map layout (source notebook cell 37)
# ---------------------------------------------------------------------------
MAP_EXTENT               = [-200, 200, -62, 86]
INSIDE_BASE_DIAMETER_DEG = 14.0
INSIDE_SPACING_FACTOR    = 0.75
OUTSIDE_SIZE             = 0.085
OUTSIDE_GAP              = 0.004
OUTSIDE_MARGIN           = 0.018
OUTSIDE_ORIENTATION      = 'auto'
OUTSIDE_MAX_PER_ROW      = 14
LABEL_FONT_SIZE          = 8.0
CONNECTOR_LW             = 0.85


def to_polar(ds):
    """The aggregate stores aspect in DEGREES; matplotlib polar axes want
    radians. Returns the dataset (or array) with the aspect coordinate in
    radians — call once after aggregate.open_aggregate / collapse.
    Also expects the bin MEANS as plain variables (the unified schema;
    the pre-2026 ``.sel(statistic='mean')`` is gone) and, for the triplets'
    third panel, ``runoff_onset_elev_relative`` computed by the caller via
    aggregate.elevation_relative."""
    return ds.assign_coords(aspect=np.deg2rad(ds["aspect"]))


def count_valid_obs(row):
    return row.notna().sum()


def wrap_labels(ax, width=30):
    labels = [textwrap.fill(label.get_text(), width=width) for label in ax.get_yticklabels()]
    ax.set_yticklabels(labels)


def get_sorted_ranges(ds, continents):
    """Range names in ``continents``, north to south. Boolean selection (not
    ``where(drop=True)``, which errors when no range matches — a partially
    processed version)."""
    sel = ds['mountain_range'][ds['continent'].isin(continents).values]
    if sel.size == 0:
        return np.array([], dtype=str)
    return ds.sel(mountain_range=sel).sortby('centroid_latitude', ascending=False)['mountain_range'].values


def create_gmba_exists_gdf(gmba_gdf, mountains_ds):
    """Create a filtered GMBA geodataframe containing only ranges that exist in mountains_ds"""
    ranges_in_metadata = mountains_ds['mountain_range'].values
    
    # Filter GMBA to only include ranges that exist in range_metadata
    mask = (gmba_gdf['MapName'].isin(ranges_in_metadata)) | (gmba_gdf['Level_04'].isin(ranges_in_metadata))
    gmba_exists_gdf = gmba_gdf[mask].copy()
    
    # Add a standardized name column for easier matching
    gmba_exists_gdf['standard_name'] = gmba_exists_gdf.apply(
        lambda row: row['MapName'] if row['MapName'] in ranges_in_metadata else row['Level_04'], 
        axis=1
    )
    
    return gmba_exists_gdf


def plot_anomaly_heatmap_panels(anomalies_df, water_years, min_valid_obs=5,
                                vmin=-30, vmax=30):
    """Three-panel per-range annual-anomaly heatmaps (Americas / Europe+
    Africa / Asia+Oceania), ranges sorted north to south, rows with
    <= ``min_valid_obs`` valid years dropped. ``anomalies_df`` comes from
    ``stats.build_anomalies_df``; ``water_years`` from the dataset.
    Extracted from the block duplicated across both mountain_ranges
    notebooks. Returns the figure (caller saves it)."""
    import seaborn as sns

    water_years = list(water_years)
    fig, axes = plt.subplots(1, 3, figsize=(20, 20))
    panels = [
        ('North America and South America', ['North America', 'South America']),
        ('Europe and Africa', ['Europe', 'Africa']),
        ('Asia and Oceania', ['Asia', 'Oceania', 'Australia']),
    ]
    for k, (ax, (title, continents)) in enumerate(zip(axes, panels)):
        parts = [anomalies_df[anomalies_df['continent'] == c]
                 .sort_values('latitude', ascending=False)
                 for c in continents]
        panel_df = pd.concat(parts)
        panel_df = panel_df[
            panel_df[water_years].apply(count_valid_obs, axis=1) > min_valid_obs]
        sns.heatmap(panel_df.set_index('name')[water_years], ax=ax,
                    cmap='RdBu', vmin=vmin, vmax=vmax,
                    cbar_kws={'label': 'Days'}, linewidths=1, square=True,
                    cbar=(k == len(panels) - 1))
        ax.set_title(title)
        ax.set_facecolor('darkgrey')
        ax.set_xlabel('Water Year')
        ax.set_ylabel('')
    fig.suptitle('Snowmelt runoff onset anomalies by mountain range')
    return fig


def create_mini_polar_plot(data, metric, cmap, vmin, vmax, dpi=160):
    """Return RGBA array of polar panel with hard circular alpha (no faint halo)."""
    fig_temp = plt.figure(figsize=(2.0, 2.0), dpi=dpi)
    ax_temp  = fig_temp.add_subplot(111, projection='polar')

    theta, r = np.meshgrid(data.aspect.values, data.elevation.values)
    ax_temp.pcolormesh(theta, r, data[metric].values,
                       cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')

    ax_temp.set_facecolor('darkgrey')
    ax_temp.set_theta_zero_location('N')
    ax_temp.set_theta_direction(-1)

    for radius in major_tick_radii:
        ax_temp.plot(np.linspace(0, 2*np.pi, 160),
                     np.ones(160)*radius, color='black', lw=1.1)

    ax_temp.set_thetagrids([0,45,90,135,180,225,270,315,360], labels=['']*9)
    ax_temp.set_rgrids(rgrid_vals, labels=rgrid_labels_blank,
                       angle=label_angle, fontsize=6)
    ax_temp.xaxis.grid(True, linestyle=':', color='gray', alpha=0.75, lw=0.9)
    [s.set_linewidth(1.1) for s in ax_temp.spines.values()]
    ax_temp.set_thetamin(0); ax_temp.set_thetamax(360)
    ax_temp.set_rlim(bottom=8000, top=0)
    ax_temp.set_xticks([]); ax_temp.set_yticks([])

    fig_temp.canvas.draw()
    buf = np.frombuffer(fig_temp.canvas.buffer_rgba(), dtype=np.uint8)
    buf = buf.reshape(fig_temp.canvas.get_width_height()[::-1] + (4,))

    # Hard circular mask (remove all outside pixels fully)
    h, w = buf.shape[:2]
    cx, cy = w//2, h//2
    rad = min(w, h)//2 - 8
    yy, xx = np.ogrid[:h, :w]
    mask = (xx-cx)**2 + (yy-cy)**2 > rad**2
    buf[mask] = 0

    plt.close(fig_temp)
    return buf


def plot_mountain_range_anomalies(mountains_ds, mountain_ranges, mountain_params):
    n_years = len(mountains_ds.water_year)
    n_cols = n_years + 1  # +1 for median column
    n_rows = len(mountain_ranges)

    # Increase left margin for mountain range labels and reduce right margin
    left_margin = 0.05
    right_margin = 0.0  # Reduced right margin to minimize empty space
    
    fig = plt.figure(figsize=(16, 1.5 * len(mountain_ranges)), dpi=300)
    
    # Base axes with adjusted margins and hidden border
    ax0 = plt.axes([0.0, 0, 1, 1])
    ax0.set_xticks([])
    ax0.set_yticks([])
    ax0.spines['top'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    
    # Adjust vertical line positions
    x_start = left_margin
    x_width = (1.0 - x_start - right_margin) / n_cols  # Using defined right margin
    ax0_adjustment = 0.0004#2
    ax0_offset = -0.0022
    #x_width = x_width - 0.02


    # Add alternating shading for anomaly columns
    for i in range(n_years):
        if i % 2 == 1:  # Even years get shading
            x_left = x_start + (i + 1) * (x_width+ax0_adjustment) + ax0_offset
            x_right = x_start + (i + 2) * (x_width+ax0_adjustment) + ax0_offset
            ax0.axvspan(x_left, x_right, color='gray', alpha=0.2, zorder=0)
    

    
    # Draw vertical lines with precise positioning
    for i in range(n_cols):
        x = x_start + (i) * (x_width+ax0_adjustment)  + ax0_offset
        if i == 0:
            # Solid line for first divider
            ax0.axvline(x=x, color='black', linewidth=1, zorder=1)
        elif i == 1:
            ax0.axvline(x=x, color='black', linewidth=1, zorder=1)
        else:
            # Dotted lines with adjusted spacing
            #ax0.axvline(x=x - 0.001, color='black', linestyle=':', linewidth=1, zorder=1)
            ax0.axvline(x=x, color='black', linestyle=':', linewidth=1, zorder=1)
    
    ax0.set_xlim([0,1])
    
    # Adjusted gridspec with better margins
    gs = fig.add_gridspec(n_rows, n_cols,
                         hspace=0.04,  # Increased spacing between rows
                         wspace=0.06,
                         left=left_margin,
                         right=(1.0 - right_margin),    # Using defined right margin
                         bottom=0.01,   # Added bottom margin
                         top=0.95)      # Adjusted top margin

    # Improved title positioning
    gs.figure.text(x_start + x_width/2, 0.98, '10-yr median \nrunoff onset ',
                  ha='center', va='center', fontsize=14)
    
    # Improved year labels
    for i, year in enumerate(mountains_ds.water_year.values):
        x_pos = x_start + (i + 1.5) * x_width
        gs.figure.text(x_pos, 0.98, f'WY {year}\n anomaly',
                      ha='center', va='center', fontsize=14)

    # Mountain range labels with smart text wrapping based on rendered width
    for i, location in enumerate(mountain_ranges):
        words = location.split()
        lines = []
        current_line = []
        current_width = 0
        target_width = gs.figure.get_figheight() * 0.13  # Adjust this factor to control wrap width
        
        # Create a dummy text to measure string widths
        test_text = gs.figure.text(0, 0, '', fontsize=10)
        
        for word in words:
            # Get the width of the current word
            test_text.set_text(word)
            word_width = test_text.get_window_extent(renderer=gs.figure.canvas.get_renderer()).width
            
            # Get the width of current line + new word
            test_line = ' '.join(current_line + [word])
            test_text.set_text(test_line)
            line_width = test_text.get_window_extent(renderer=gs.figure.canvas.get_renderer()).width
            
            if current_line and line_width > target_width:
                # Add current line to lines and start new line with current word
                lines.append(' '.join(current_line))
                current_line = [word]
                current_width = word_width
            else:
                current_line.append(word)
                current_width = line_width
        
        # Add the last line
        if current_line:
            lines.append(' '.join(current_line))
        
        # Remove the dummy text
        test_text.remove()
        
        # Join lines with newlines
        wrapped_text = '\n'.join(lines)
        
        gs.figure.text(left_margin/2-0.01, 0.89 - 0.95*(i)/n_rows, wrapped_text,
                      rotation=90, va='center', ha='center', fontsize=14)
        
        for j in range(n_cols):
            ax = gs.figure.add_subplot(gs[i, j], projection='polar')
            
            if j == 0:  # Median column
                range_data = mountains_ds.sel(mountain_range=location)
                median_data = range_data['runoff_onset'].median(dim='water_year').values
                theta, r = np.meshgrid(range_data.aspect.values, range_data.elevation.values)
                im_median = ax.pcolormesh(theta, r, median_data,
                                        vmin=80, vmax=300,
                                        cmap='viridis', shading='auto')
            else:
                year = mountains_ds.water_year.values[j-1]
                range_data = mountains_ds.sel(mountain_range=location)
                data_anomaly = range_data['runoff_onset_anomaly'].sel(water_year=year).values
                im_anomaly = ax.pcolormesh(theta, r, data_anomaly,
                                         vmin=-30, vmax=30,
                                         cmap='RdBu', shading='auto')

            params = mountain_params.get(location, {
                'bottom_r': 8000,
                'top_r': 0
            })
            
            ax.set_facecolor('darkgrey')
            ax.set_theta_zero_location('N')
            ax.set_theta_direction(-1)
            
            # Draw major tick radius circles
            for major_tick_radius in major_tick_radii:
                ax.plot(np.linspace(0, 2 * np.pi, 100),
                       np.ones(100) * major_tick_radius,
                       color='black', linestyle='-', alpha=1, linewidth=1)
                

            def roundup(x):
                return x if x % 1000 == 0 else x + 1000 - x % 1000
            
            bottom_r = roundup(mountains_ds['runoff_onset_median'].sel(mountain_range=location).where(lambda x: x > 0, drop=True).elevation.max().values)

            ax.scatter(x=0, y=bottom_r, marker='.', color='black', s=100)
            ax.set_xlabel('')
            ax.set_ylabel('')
            
            # Improved grid and label settings
            ax.set_thetagrids([0, 45, 90, 135, 180, 225, 270, 315],
                             labels=[''] * 8)
            ax.set_rgrids(rgrid_vals, labels=rgrid_labels_blank,
                         angle=45,  # Adjusted label angle
                         fontsize=7,
                         ha='center', va='center')
            
            # Enhanced grid appearance
            ax.xaxis.grid(True, which="both", linestyle=":", color='gray',
                         alpha=0.8, linewidth=0.8)
            ax.yaxis.grid(True, which="both", linestyle=":", color='black',
                         alpha=0.8, linewidth=0.8)
            
            # Hide the outer circle border and adjust other spines
            ax.spines['polar'].set_visible(False)
            
            # Adjust other visual elements
            for spine in ax.spines.values():
                if spine != ax.spines['polar']:
                    spine.set_linewidth(1.0)
                    
            ax.set_thetamin(0)
            ax.set_thetamax(360)
            ax.set_rlim(bottom=bottom_r, top=params['top_r'])

    return fig


def _edge_flags(lon, lat, extent):
    mn_lon, mx_lon, mn_lat, mx_lat = extent
    return {
        'top':    lat >  mx_lat,
        'bottom': lat <  mn_lat,
        'left':   lon <  mn_lon,
        'right':  lon >  mx_lon
    }


def _lonlat_to_axes_fraction(ax, lon, lat):
    """Map lon/lat to axes fraction (0-1)^2 given current projection/extent."""
    proj_pt = ax.projection.transform_point(lon, lat, ccrs.PlateCarree())
    disp    = ax.transData.transform(proj_pt)
    return ax.transAxes.inverted().transform(disp)


def _outside_triplet_centers(edge_flags, anchor_axes_xy):
    """
    Compute three centers (axes coords) for outside triplet based on anchor axes position.
    anchor_axes_xy is (ax_x, ax_y) of the *anchor* point projected into axes space (clamped inside 0..1).
    For top/bottom we optionally pack columns by snapping x to discrete slots to avoid crowding.
    """
    ax_x, ax_y = anchor_axes_xy
    # Decide primary edge (priority vertical -> horizontal)
    if edge_flags['top']:    edge = 'top'
    elif edge_flags['bottom']: edge = 'bottom'
    elif edge_flags['left']: edge = 'left'
    elif edge_flags['right']: edge = 'right'
    else: return None

    # Orientation policy
    if OUTSIDE_ORIENTATION == 'horizontal':
        orient = 'horizontal'
    elif OUTSIDE_ORIENTATION == 'vertical':
        orient = 'vertical'
    else:
        orient = 'horizontal' if edge in ('top','bottom') else 'vertical'

    # Base line position
    if edge == 'top':
        base_y = 1 + OUTSIDE_MARGIN
        base_x = np.clip(ax_x, 0.03, 0.97)
    elif edge == 'bottom':
        base_y = -OUTSIDE_MARGIN
        base_x = np.clip(ax_x, 0.03, 0.97)
    elif edge == 'left':
        base_x = -OUTSIDE_MARGIN
        base_y = np.clip(ax_y, 0.03, 0.97)
    else:
        base_x = 1 + OUTSIDE_MARGIN
        base_y = np.clip(ax_y, 0.03, 0.97)

    # Simple packing for top/bottom to reduce clustering
    if edge in ('top','bottom') and OUTSIDE_MAX_PER_ROW > 2:
        slot_width = 1.0 / OUTSIDE_MAX_PER_ROW
        slot_index = int(np.clip(np.round(base_x / slot_width), 0, OUTSIDE_MAX_PER_ROW-1))
        base_x = (slot_index + 0.5) * slot_width

    spacing = OUTSIDE_SIZE + OUTSIDE_GAP
    if orient == 'horizontal':
        centers = np.column_stack((
            [base_x - spacing, base_x, base_x + spacing],
            [base_y, base_y, base_y]
        ))
    else:  # vertical
        centers = np.column_stack((
            [base_x, base_x, base_x],
            [base_y + spacing, base_y, base_y - spacing]
        ))
    return centers, edge, orient


# Anchor (lon, lat) per display range for the exploratory Equal-Earth triplet map: anchors
# outside MAP_EXTENT place the triplet along the map edge, inside ones sit
# geographically. Recovered 2026-08-25 from the pre-2026 combined notebook
# (kept outside this repository) — the defining cell was lost in the 2026-08-24
# debris purge while geo_and_topo_analysis.ipynb kept calling it.
#
# NOTE (2026-09-01): this Equal-Earth map was never the polar-triplet world map that was published. The
# published label set and hand-placed anchors (41 ranges, Robinson metres; nine
# differ from the keys below each way) are in analyses/mountain_ranges/label_layout.csv,
# read by gsro_analysis.world_maps, which draws the actual figure
# (topography_triplet_composite_figure.ipynb). Kept as the exploratory analogue.
TRIPLET_MAP_ANCHORS = {
    # North America - Alaska/Arctic
    'Brooks Range': (-145, 80),              # Moved further west
    'Mackenzie Mountains': (-100, 70),       # Added Mackenzie Mountains
    'South-Central Alaska': (-165, 62),      # Moved west, lower
    'Saint Elias Mountains': (-125, 60),     # Positioned between Alaska ranges
    'Alaska Range': (-155, 58),              # Moved to avoid overlap
    'Kolyma Mountains': (155, 68),           # Moved to avoid Verkhoyansk overlap
    
    # North America - Western US/Canada
    'Coast Mountains': (-120, 52),           # Positioned in Pacific
    'Cascade Range': (-135, 35),             # Moved east slightly
    'Columbia Mountains': (-120, 52),        # Moved east
    'Sierra Nevada': (-135, 20),             # Moved west into Pacific
    'Great Basin Ranges': (-110, 25),        # Moved east
    'Appalachian Mountains': (-60, 25),      # Added Appalachian Mountains
    'Southern Rocky Mountains': (-95, 40),  # Added Southern Rocky Mountains

    # Europe/Nordic
    'Iceland': (-40, 75),                    # Moved further west
    'Northern Scandes': (0, 75),             # Moved north to separate from Southern
    'Southern Scandes': (0, 60),             # Kept in Scandinavia
    'European Alps': (15, 30),                # Moved west into Atlantic
    'Pyrenees': (-35, 50),                   # Moved west
    'Ural Mountains': (50, 75),              # Moved west slightly
    'South European Highlands': (-25, 35), # Added Southern European Highlands

    # Asia - Central/Western
    'Caucasus Mountains': (45, 60),          # Moved west slightly
    'Taurus Mountains': (25, 15),            # Moved west
    'Pamir Mountains': (70, 40),             # Slight adjustment
    'Hindu Kush': (65, 38),                  # Moved west slightly
    'Karakoram': (80, 15),                   # Moved south to separate from Hindu Kush
    'Himalaya': (95, 0),                    # Moved south
    'Tian Shan': (75, 60),                   # Moved east
    'Kunlun Mountains': (120, 30),            # Moved east
    'Mongolian Highlands': (130, 45),        # Added Mongolian Highlands
    'Central Siberian Plateau': (90, 75),   # Added Central Siberian Plateau
    'Zagros Mountains': (50, 10),             # Added Zagros Mountains

    # Asia - Far East
    'Verkhoyansk Range': (120, 70),          # Moved west to avoid Kolyma overlap
    'Kamchatka Peninsula': (160, 55),        # Added Kamchatka Peninsula

    # Africa
    'High Atlas Range': (-8, 15),            # Moved west into Atlantic
    
    # South America
    'Dry Andes': (-100, -20),                 # Moved east
    'Meseta Patagónica': (-100, -35),         # Moved east, better separation
    'Patagonian Andes': (-100, -50),          # Moved east, further south
    'Cordillera Occidental (Central Andes)': (-100, -5),  # Moved west
    
    # Oceania
    'Southern Alps': (165, -42),             # Moved west slightly
    'Kaikoura Ranges': (175, -38),           # Moved north to separate from Southern Alps
    'Southern Great Dividing Range': (155, -32)  # Moved west into ocean
}


def plot_mountain_ranges_on_map(mountains_ds, test_ranges_dict, gmba_exists_gdf):
    """
    Global map with 3-panel polar triplets:
      - If anchor point (in test_ranges_dict) lies outside MAP_EXTENT -> place triplet outside edge.
      - Otherwise place triplet geographically at anchor (size scaled by latitude).
    """
    map_projection  = ccrs.EqualEarth()
    data_projection = ccrs.PlateCarree()

    fig = plt.figure(figsize=(16, 9), dpi=300)
    ax  = fig.add_subplot(111, projection=map_projection)
    ax.set_position([0.055, 0.055, 0.87, 0.87])

    # Base map
    ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.75)
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.75)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.45)

    # Background ranges
    gmba_exists_gdf.plot(ax=ax, color='gray', alpha=0.25,
                         edgecolor='none', transform=data_projection)

    sel_mask = gmba_exists_gdf['standard_name'].isin(test_ranges_dict.keys())
    gmba_exists_gdf[sel_mask].boundary.plot(ax=ax, color='black',
                                            linewidth=0.35, alpha=0.55,
                                            transform=data_projection)

    metrics = [
        ('runoff_onset_median',       'viridis', 100, 300),
        ('runoff_onset_mad',          'Reds',    0,   30),
        ('runoff_onset_elev_relative','RdBu',   -15,  15),
    ]

    for mtn_name, anchor in test_ranges_dict.items():
        if (mtn_name not in mountains_ds.mountain_range.values or
            mtn_name not in mountains_ds['mountain_range'].values):
            continue

        ds_range = mountains_ds.sel(mountain_range=mtn_name)
        cen_lon = float(ds_range['centroid_longitude'].values)
        cen_lat = float(ds_range['centroid_latitude'].values)
        anchor_lon, anchor_lat = anchor

        flags  = _edge_flags(anchor_lon, anchor_lat, MAP_EXTENT)
        is_out = any(flags.values())

        if is_out:
            # Derive axes fraction from *anchor* (clip into current box for relative positioning)
            # If anchor is far away, first clamp to numeric outer shell so transform is stable.
            clamp_lon = np.clip(anchor_lon, MAP_EXTENT[0]-1, MAP_EXTENT[1]+1)
            clamp_lat = np.clip(anchor_lat, MAP_EXTENT[2]-1, MAP_EXTENT[3]+1)
            ax_frac = _lonlat_to_axes_fraction(ax, clamp_lon, clamp_lat)
            centers_info = _outside_triplet_centers(flags, ax_frac)
            if centers_info is None:
                continue
            centers_axes, edge, orient = centers_info

            inset_axes_list = []
            for (metric, cmap, vmin, vmax), (cx, cy) in zip(metrics, centers_axes):
                mini = create_mini_polar_plot(ds_range, metric, cmap, vmin, vmax)
                llx  = cx - OUTSIDE_SIZE/2
                lly  = cy - OUTSIDE_SIZE/2
                inner = inset_axes(ax,
                                   width=OUTSIDE_SIZE,
                                   height=OUTSIDE_SIZE,
                                   bbox_to_anchor=(llx, lly, OUTSIDE_SIZE, OUTSIDE_SIZE),
                                   bbox_transform=ax.transAxes,
                                   borderpad=0)
                inner.imshow(mini, aspect='equal', interpolation='nearest')
                inner.axis('off')
                inset_axes_list.append(inner)

            # Build connector: centroid -> middle outside circle center in data coords
            mid_ax = inset_axes_list[1]
            mid_disp = mid_ax.transAxes.transform((0.5, 0.5))
            mid_data = ax.transData.inverted().transform(mid_disp)
            # Draw line
            ax.plot([cen_lon, mid_data[0]], [cen_lat, mid_data[1]],
                    transform=data_projection, color='black',
                    lw=CONNECTOR_LW, zorder=9, clip_on=False)

            # Centroid marker
            ax.scatter(cen_lon, cen_lat, s=20, color='black',
                       edgecolor='white', linewidth=0.55,
                       transform=data_projection, zorder=10)

            # Label above/below/side of middle circle
            mid_cx, mid_cy = centers_axes[1]
            if edge == 'top':
                tx, ty, va = mid_cx, mid_cy + OUTSIDE_SIZE/2 + 0.006, 'bottom'
            elif edge == 'bottom':
                tx, ty, va = mid_cx, mid_cy - OUTSIDE_SIZE/2 - 0.006, 'top'
            elif edge == 'left':
                tx, ty, va = mid_cx, mid_cy - OUTSIDE_SIZE/2 - 0.01, 'top'
            else:  # right
                tx, ty, va = mid_cx, mid_cy - OUTSIDE_SIZE/2 - 0.01, 'top'

            ax.text(tx, ty, mtn_name,
                    transform=ax.transAxes,
                    ha='center', va=va,
                    fontsize=LABEL_FONT_SIZE, fontweight='bold',
                    bbox=dict(boxstyle="round,pad=0.25",
                              facecolor='white', edgecolor='black',
                              linewidth=0.4, alpha=0.92),
                    zorder=11, clip_on=False)

        else:
            # Inside placement (geographic)
            lat_scale = 1.0 + abs(cen_lat)/90.0 * 0.22
            diam      = INSIDE_BASE_DIAMETER_DEG * lat_scale
            spacing   = diam * INSIDE_SPACING_FACTOR
            centers_deg = [
                (anchor_lon - spacing, anchor_lat),
                (anchor_lon,           anchor_lat),
                (anchor_lon + spacing, anchor_lat)
            ]

            # Connector centroid -> anchor middle
            ax.plot([cen_lon, anchor_lon], [cen_lat, anchor_lat],
                    color='black', lw=CONNECTOR_LW, alpha=0.85,
                    transform=data_projection, zorder=6)

            ax.scatter(cen_lon, cen_lat, s=20, color='black',
                       edgecolor='white', linewidth=0.55,
                       transform=data_projection, zorder=7)

            # Projection size (horizontal) for consistent circle size in projection space
            x1, y1 = map_projection.transform_point(anchor_lon - diam/2, anchor_lat, data_projection)
            x2, y2 = map_projection.transform_point(anchor_lon + diam/2, anchor_lat, data_projection)
            proj_size = abs(x2 - x1)
            half = proj_size / 2

            for (metric, cmap, vmin, vmax), (lon_c, lat_c) in zip(metrics, centers_deg):
                mini = create_mini_polar_plot(ds_range, metric, cmap, vmin, vmax)
                px, py = map_projection.transform_point(lon_c, lat_c, data_projection)
                ax.imshow(mini,
                          extent=[px - half, px + half, py - half, py + half],
                          transform=map_projection, zorder=8,
                          interpolation='nearest', clip_on=False)

            # Label under triplet
            lx, ly = map_projection.transform_point(anchor_lon,
                                                    anchor_lat - diam*0.52,
                                                    data_projection)
            ax.text(lx, ly, mtn_name,
                    transform=map_projection,
                    ha='center', va='top',
                    fontsize=LABEL_FONT_SIZE, fontweight='bold',
                    bbox=dict(boxstyle="round,pad=0.22",
                              facecolor='white', edgecolor='black',
                              linewidth=0.4, alpha=0.92),
                    zorder=9, clip_on=False)

    ax.set_extent(MAP_EXTENT, crs=data_projection)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                      xlocs=[-180,-120,-60,0,60,120,180],
                      ylocs=[-60,-30,0,30,60],
                      linestyle='--', linewidth=0.5, alpha=0.6)
    gl.top_labels = False
    gl.right_labels = False
    return fig, ax
