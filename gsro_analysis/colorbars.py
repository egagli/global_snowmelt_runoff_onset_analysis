"""Annotated colorbars for the analysis figures.

Thin presets over the production repo's colorbar builders
(``global_snowmelt_runoff_onset.plot_utils.create_month_colorbar`` and
``create_diverging_colorbar``): those own the drawing (bar, ticks, month
boundaries, white-on-black end texts, auto-shrinking to the axes width); this
module owns the *analysis-specific* variables, ranges and wording — the spring
temperature sensitivity, the elevation delay, the sunny–shaded difference, and
the mean-anomaly / MAD / median-onset wording used on the figures.

Until 2026-09 these bars were rendered standalone by a loose ``create_colorbar.py``
outside the repository and pasted onto the notebook renders in slides. Now the notebooks
draw them in place (pass ``ax=`` and the ``EMBEDDED*`` size kwargs); the two
world maps (``gsro_analysis.world_maps``) call the presets with ``ax=`` too and
scale the fonts to the picture size.

Per production convention 3 no default label carries a year count: pass
``n_years=len(ds.water_year)`` and it is interpolated.

Usage::

    from gsro_analysis import colorbars
    cax = fig.add_axes([0.3, 0.0, 0.4, 0.03])
    colorbars.anomaly(cax, n_years=len(ds.water_year), **colorbars.EMBEDDED)
"""

from global_snowmelt_runoff_onset import plot_utils as pu

# Font sizes for a bar embedded in a multi-panel figure; the builders default to
# standalone 8 x 3 inch bars (label 15, ticks 15, end text 18).
EMBEDDED = dict(label_fontsize=10, tick_labelsize=9, text_fontsize=10)
EMBEDDED_MONTH = dict(label_fontsize=10, tick_labelsize=9, month_fontsize=8)


def _years(n_years, suffix='-year '):
    return f'{n_years}{suffix}' if n_years else ''


def median_onset(ax=None, n_years=None, dowy=(100, 300), hemisphere='both',
                 abbreviate=True, **kw):
    """viridis day-of-water-year bar with month boundaries and names (NH bold on
    top, SH italic below when ``hemisphere='both'``). ``abbreviate=False`` gives
    the full month names of the published standalone bar (needs an 8-inch-wide
    axes); embedded bars want the 3-letter form."""
    label = f'{_years(n_years)}median runoff onset date [day of water year]'
    return pu.create_month_colorbar(dowy[0], dowy[1], hemisphere=hemisphere,
                                    major_tick_spacing=50, minor_tick_spacing=10,
                                    ax=ax, label=label, abbreviate_month_names=abbreviate, **kw)


def anomaly(ax=None, n_years=None, what='Mountain range mean runoff onset anomaly', **kw):
    """RdBu ±30 days, 'earlier / later than the N-yr median' end texts (the regional anomaly panels' bar)."""
    n = f'{n_years}-YR ' if n_years else ''
    return pu.create_diverging_colorbar(-30, 30, cmap='RdBu', label=f'{what} [days]',
                                        ticks=[-30, -20, -10, 0, 10, 20, 30], minor_tick_spacing=5,
                                        left_text=f'◀ EARLIER THAN {n}MEDIAN',
                                        right_text=f'LATER THAN {n}MEDIAN ▶', ax=ax, **kw)


def mad(ax=None, n_years=None, **kw):
    """Reds 0–30 days, 'lower / higher variability'."""
    return pu.create_diverging_colorbar(0, 30, cmap='Reds',
                                        label=f'{_years(n_years)}median absolute deviation [days]',
                                        ticks=[0, 5, 10, 15, 20, 25, 30], minor_tick_spacing=2.5,
                                        left_text='lower variability', right_text='higher variability',
                                        ax=ax, **kw)


def sunny_shaded(ax=None, **kw):
    """PuOr ±30 days: sunny (CHILI warm) minus shaded (cool) onset timing."""
    return pu.create_diverging_colorbar(-30, 30, cmap='PuOr',
                                        label='Runoff onset timing difference between sunny and shaded areas [days]',
                                        ticks=[-30, -20, -10, 0, 10, 20, 30], minor_tick_spacing=5,
                                        left_text='sunny areas melt first', right_text='shaded areas melt first',
                                        ax=ax, **kw)


def temperature_sensitivity(ax=None, **kw):
    """YlOrRd (reversed) −12–0 days/°C: the temperature-sensitivity world map's choropleth ramp."""
    return pu.create_diverging_colorbar(-12, 0, cmap='YlOrRd_r',
                                        label='Runoff onset anomaly per 1°C increase in spring temp. [days/°C]',
                                        ticks=[-12, -9, -6, -3, 0], minor_tick_spacing=1,
                                        left_text='higher sensitivity', right_text='lower sensitivity',
                                        ax=ax, **kw)


def elevation_delay(ax=None, **kw):
    """YlGnBu 0–8 days/100 m: the polar-triplet world map's lapse-rate choropleth ramp."""
    return pu.create_diverging_colorbar(0, 8, cmap='YlGnBu',
                                        label='Delay in runoff onset per 100 meter increase in elevation [days/100m]',
                                        ticks=[0, 2, 4, 6, 8], minor_tick_spacing=1,
                                        left_text='weaker delay', right_text='stronger delay',
                                        ax=ax, **kw)
