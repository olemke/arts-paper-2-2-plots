import os
import sys
from glob import glob

import matplotlib.pyplot as plt
import numpy
import typhon.arts.xml as axml
import typhon.constants


def read_planet_data(planets):
    refdata = {}
    for planet, id in planets:
        refdata[planet] = {}
        refdata[planet]['t'] = axml.load(
            os.path.join(sys.argv[1], 't_field.' + planet + '.xml'))
        refdata[planet]['refindex'] = [axml.load(ifile) for ifile in
                                       sorted(glob(
                                           os.path.join(sys.argv[1],
                                                        "refindexair." + planet + "*.xml")))]
        refdata[planet]['p_grid'] = axml.load(
            os.path.join(sys.argv[1], "p_grid." + planet + ".xml"))

    return refdata


def plot_refindex_p(planets, ax=None):
    if ax is None:
        fig, ax = subplots()

    ax.set_yscale('log')
    ax.set_xscale('log')
    y_min = 1e100
    y_max = -1e100

    for planet, zorder in planets:
        y = refdata[planet]['p_grid']
        x = refdata[planet]['refindex']
        y_max1 = numpy.max(y)
        y_min1 = numpy.min(y)
        if y_max < y_max1: y_max = y_max1
        if y_min > y_min1: y_min = y_min1

        lines, = ax.plot(x, y, linewidth=3, label=planet, zorder=zorder)
        lines.set_solid_capstyle('round')

    ax.set_ylim(y_max, y_min)  # implicitly invert yaxis
    ax.set_xscale('log')
    ax.set_xlabel('Refractivity')
    ax.set_ylabel('Pressure [hPa]')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ylims = ax.get_ylim()
    ax.set_ylim((ylims[0] * 1.5, ylims[1] / 1.5))
    ax.set_yticks([1e0, 1e2, 1e4, 1e6])
    ax.set_xticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2])
    ax.yaxis.set_major_formatter(typhon.plots.HectoPascalFormatter())
    for i in (
            ('Earth', (0.22, 0.84), (-0.14, -0.12)),
            ('Mars', (0.34, 0.75), (0.115, 0.125)),
            ('Venus', (0.82, 0.2), (0.14, 0.12)),
            ('Jupiter', (0.75, 0.23), (-0.17, -0.1)),
    ):
        ax.annotate(i[0], xy=i[1],
                    xytext=(i[1][0] + i[2][0], i[1][1] + i[2][1]),
                    xycoords='axes fraction', textcoords='axes fraction',
                    ha='center', va='center',
                    arrowprops=dict(facecolor='black', shrink=0.1, width=1,
                                    headlength=6, headwidth=6),
                    )
    # ax.legend()


def plot_refindex_n(planets, ax=None):
    if ax is None:
        fig, ax = subplots()

    ax.set_yscale('log')
    ax.set_xscale('log')
    y_min = 1e100
    y_max = -1e100

    nscale = 1
    for planet, zorder in planets:
        y = refdata[planet]['p_grid'] / refdata[planet][
            't'].flatten() / typhon.constants.boltzmann / nscale
        x = refdata[planet]['refindex']
        y_max1 = numpy.max(y)
        y_min1 = numpy.min(y)
        if y_max < y_max1: y_max = y_max1
        if y_min > y_min1: y_min = y_min1

        lines, = ax.plot(x, y, linewidth=3, label=planet, zorder=zorder)
        lines.set_solid_capstyle('round')

    # n for dry air:
    # ax.plot([1e-10, 1e-2], [0.025*1e27, 0.025*1e27], linewidth=3)
    ax.set_ylim(y_max, y_min)  # implicitly invert yaxis
    ax.set_xscale('log')
    ax.set_xlabel('Refractivity')
    ax.set_ylabel(
        'Number density $\\left[\\frac{1}{\\mathrm{m}^3}\\right]$')
    # 'Number density $\\left[\\frac{10^{27}}{\\mathrm{m}^3}\\right]$')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ylims = ax.get_ylim()
    ax.set_ylim((ylims[0] * 1.5, ylims[1] / 1.5))
    ax.set_yticks([1e20 / nscale, 1e22 / nscale, 1e24 / nscale, 1e26 / nscale])
    ax.set_xticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2])
    for i in (
            ('Earth', (0.205, 0.84), (-0.14, -0.12)),
            ('Mars', (0.32, 0.75), (0.115, 0.125)),
            ('Venus', (0.8, 0.2), (0.14, 0.12)),
            ('Jupiter', (0.70, 0.23), (-0.17, -0.1)),
    ):
        ax.annotate(i[0], xy=i[1],
                    xytext=(i[1][0] + i[2][0], i[1][1] + i[2][1]),
                    xycoords='axes fraction', textcoords='axes fraction',
                    ha='center', va='center',
                    arrowprops=dict(facecolor='black', shrink=0.1, width=1,
                                    headlength=6, headwidth=6),
                    )
    # ax.legend()


planets = (('Earth', 2), ('Mars', 3), ('Venus', 0), ('Jupiter', 1))
refdata = read_planet_data(planets)

fig, ax = plt.subplots(1, 1, figsize=(4.9, 3.4))
plot_refindex_n(planets, ax=ax)
fig.tight_layout(pad=1)
plt.show()
fig.savefig(os.path.join(sys.argv[1], 'refractivity_n.pdf'), dpi=300)

fig, ax = plt.subplots(1, 1, figsize=(4.9, 3.4))
plot_refindex_p(planets, ax=ax)
fig.tight_layout(pad=1)
plt.show()
fig.savefig(os.path.join(sys.argv[1], 'refractivity_p.pdf'), dpi=300)
