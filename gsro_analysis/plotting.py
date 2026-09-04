"""Polar-panel helpers for the mountain-range notebooks.

Every polar panel in this repository draws the same thing: aspect around the circle (north up,
clockwise), elevation along the radius (0 m at the rim, the range's top elevation at the centre),
black rings every 1000 m, a dark-grey face. :func:`style_polar_axes` applies that look to an axes
after the data has been drawn with ``pcolormesh``; :func:`plot_mountain_range_anomalies` is the
one multi-panel polar figure kept as a function (a row per range, the median plus one anomaly
panel per water year, with the column dividers and labels that make it a table).

The cubes store aspect in DEGREES; polar axes want radians, so the notebooks convert once::

    mountain_ranges_ds = mountain_ranges_ds.assign_coords(aspect=np.deg2rad(mountain_ranges_ds['aspect']))
"""

import matplotlib.pyplot as plt
import numpy as np

# radial grid: elevation in metres (radius = elevation)
major_tick_radii = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000]
rgrid_vals = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000]
rgrid_labels_blank = [''] * len(rgrid_vals)
rgrid_labels = ['0m', '', '1000m', '', '2000m', '', '3000m', '', '4000m', '', '5000m', '', '6000m', '', '7000m', '',
                '8000m', '', '9000m']
label_angle = 225
ASPECT_TICKS = [0, 45, 90, 135, 180, 225, 270, 315, 360]
ASPECT_LABELS = ['N', '', 'E', '', 'S', '', 'W', '', 'N']


def style_polar_axes(ax, bottom_r=8000, top_r=0, aspect_labels=False, elevation_labels=False,
                     center_dot=True, ring_linewidth=1.0, spine_linewidth=1.2):
    """The shared look of a polar elevation x aspect panel, applied after the data is drawn.

    ``bottom_r`` is the elevation at the CENTRE of the circle and ``top_r`` the elevation at the
    rim (the radial axis is inverted so low elevations sit outside). ``aspect_labels`` writes
    N/E/S/W on the theta grid, ``elevation_labels`` the metre labels on the rings; both are off
    for the small multiples. ``center_dot`` marks the centre elevation with a black dot."""
    ax.set_facecolor('darkgrey')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    for radius in major_tick_radii:
        ax.plot(np.linspace(0, 2 * np.pi, 100), np.full(100, radius), color='black', linestyle='-',
                alpha=1, linewidth=ring_linewidth)
    if center_dot:
        ax.scatter(x=0, y=bottom_r, marker='.', color='black', s=100)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_thetagrids(ASPECT_TICKS, labels=ASPECT_LABELS if aspect_labels else [''] * len(ASPECT_TICKS))
    ax.set_rgrids(rgrid_vals, labels=rgrid_labels if elevation_labels else rgrid_labels_blank,
                  angle=label_angle, fontsize=7, ha='left', va='bottom', zorder=0)
    ax.xaxis.grid(True, which='both', linestyle=':', color='gray', alpha=0.8, linewidth=1)
    ax.yaxis.grid(True, which='both', linestyle=':', color='black', alpha=1, linewidth=1)
    for spine in ax.spines.values():
        spine.set_linewidth(spine_linewidth)
    ax.set_thetamin(0)
    ax.set_thetamax(360)
    ax.set_rlim(bottom=bottom_r, top=top_r)
    return ax


def plot_mountain_range_anomalies(mountain_ranges_ds, mountain_ranges, mountain_params):
    """One row per range: the median runoff onset (viridis, 80-300 DOWY) followed by one polar
    anomaly panel per water year (RdBu, +-30 days), with alternating column shading, dividers,
    column headers and wrapped range names. ``mountain_ranges_ds`` is the analyses' view of the
    cube with aspect in radians; ``mountain_params`` may give ``top_r`` per range (the centre
    elevation is rounded up to the next 1000 m above the range's highest bin with data)."""
    n_years = len(mountain_ranges_ds.water_year)
    n_cols = n_years + 1  # +1 for median column
    n_rows = len(mountain_ranges)

    left_margin = 0.05
    right_margin = 0.0

    fig = plt.figure(figsize=(16, 1.5 * len(mountain_ranges)), dpi=300)

    # base axes carrying the column shading and dividers
    ax0 = plt.axes([0.0, 0, 1, 1])
    ax0.set_xticks([])
    ax0.set_yticks([])
    for side in ('top', 'bottom', 'left', 'right'):
        ax0.spines[side].set_visible(False)

    x_start = left_margin
    x_width = (1.0 - x_start - right_margin) / n_cols
    ax0_adjustment = 0.0004
    ax0_offset = -0.0022

    for i in range(n_years):
        if i % 2 == 1:
            x_left = x_start + (i + 1) * (x_width + ax0_adjustment) + ax0_offset
            x_right = x_start + (i + 2) * (x_width + ax0_adjustment) + ax0_offset
            ax0.axvspan(x_left, x_right, color='gray', alpha=0.2, zorder=0)

    for i in range(n_cols):
        x = x_start + i * (x_width + ax0_adjustment) + ax0_offset
        ax0.axvline(x=x, color='black', linewidth=1, zorder=1, linestyle='-' if i <= 1 else ':')
    ax0.set_xlim([0, 1])

    gs = fig.add_gridspec(n_rows, n_cols, hspace=0.04, wspace=0.06, left=left_margin,
                          right=(1.0 - right_margin), bottom=0.01, top=0.95)

    gs.figure.text(x_start + x_width / 2, 0.98, f'{n_years}-yr median \nrunoff onset ',
                   ha='center', va='center', fontsize=14)
    for i, year in enumerate(mountain_ranges_ds.water_year.values):
        gs.figure.text(x_start + (i + 1.5) * x_width, 0.98, f'WY {year}\n anomaly',
                       ha='center', va='center', fontsize=14)

    def roundup(x):
        return x if x % 1000 == 0 else x + 1000 - x % 1000

    for i, location in enumerate(mountain_ranges):
        # range name wrapped to the row height, measured with the renderer
        words = location.split()
        lines, current_line = [], []
        target_width = gs.figure.get_figheight() * 0.13
        test_text = gs.figure.text(0, 0, '', fontsize=10)
        renderer = gs.figure.canvas.get_renderer()
        for word in words:
            test_text.set_text(' '.join(current_line + [word]))
            if current_line and test_text.get_window_extent(renderer=renderer).width > target_width:
                lines.append(' '.join(current_line))
                current_line = [word]
            else:
                current_line.append(word)
        if current_line:
            lines.append(' '.join(current_line))
        test_text.remove()
        gs.figure.text(left_margin / 2 - 0.01, 0.89 - 0.95 * i / n_rows, '\n'.join(lines),
                       rotation=90, va='center', ha='center', fontsize=14)

        range_ds = mountain_ranges_ds.sel(mountain_range=location)
        theta, r = np.meshgrid(range_ds.aspect.values, range_ds.elevation.values)
        bottom_r = roundup(range_ds['runoff_onset_median'].where(lambda x: x > 0, drop=True).elevation.max().values)
        top_r = mountain_params.get(location, {}).get('top_r', 0)

        for j in range(n_cols):
            ax = gs.figure.add_subplot(gs[i, j], projection='polar')
            if j == 0:
                ax.pcolormesh(theta, r, range_ds['runoff_onset'].median(dim='water_year').values,
                              vmin=80, vmax=300, cmap='viridis', shading='auto')
            else:
                year = mountain_ranges_ds.water_year.values[j - 1]
                ax.pcolormesh(theta, r, range_ds['runoff_onset_anomaly'].sel(water_year=year).values,
                              vmin=-30, vmax=30, cmap='RdBu', shading='auto')
            style_polar_axes(ax, bottom_r=bottom_r, top_r=top_r, ring_linewidth=1.0, spine_linewidth=1.0)
            ax.spines['polar'].set_visible(False)

    return fig
