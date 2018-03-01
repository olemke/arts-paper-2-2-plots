"""Create and plot absorption cross sections for different planets.

Plots the 183 GHz H2O line for different planets at 100 Pa.

Usage: python plot_xsec.py OUTDIR

The plot is created in the given output directory. OUTDIR is created if it
doesn't exist.

The following environment variables must be set:

export ARTS_INCLUDE_PATH=.:/home/oliver/arts/controlfiles/general:/home/oliver/arts-xml-data
export ARTS_BUILD_PATH=/home/oliver/arts/build

Other environment variables:
    PLT_SHOW: Set to 1 to display plot.
    PLT_NOTEX: Set to 1 to disable TeX font rendering.

Author: oliver.lemke@uni-hamburg.de
"""
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy
from matplotlib.ticker import FuncFormatter
from typhon.arts.workspace import Workspace
from typhon.plots import styles

PLANET_SETUP = {
    'Earth': {'species': ['H2O', 'O2', 'CO2', 'N2'],
              'basename': 'planets/Earth/Fascod/tropical/tropical',
              'include': 'planet_earth.arts'},
    'Mars': {'species': ['H2O', 'O2', 'CO2', 'H2', 'N2'],
             'basename': 'planets/Mars/MPS/Mars.Ls0.day.dust-medium/'
                         'Mars.Ls0.day.dust-medium.sol-avg/'
                         'Mars.Ls0.day.dust-medium.sol-avg',
             'include': 'planet_mars.arts'},
    'Venus': {'species': ['H2O', 'O2', 'CO2', 'N2'],
              'basename': 'planets/Venus/MPS/Venus.vira.day/Venus.vira.day',
              'include': 'planet_venus.arts'},
    'Jupiter': {'species': ['H2O', 'CO2', 'H2', 'He'],
                'basename': 'planets/Jupiter/MPS/Jupiter.mean/Jupiter.mean',
                'include': 'planet_jupiter.arts'},
}


def CenteredGigaHertzFormatter(center=0.):
    @FuncFormatter
    def _CenteredGigaHertzFormatter(x, pos):
        if numpy.isclose(x, center):
            return r'$\mathrm{f}_0=' + '{:.1f}$'.format(x / 1e9)
        else:
            return r'$\mathrm{f}_0' + '{:+.2g}$'.format((x - center) / 1e9)

    return _CenteredGigaHertzFormatter


def CenteredGigaHertzFormatter_bak(center=0.):
    @FuncFormatter
    def _CenteredGigaHertzFormatter(x, pos):
        if numpy.isclose(x, center):
            return '{:g}'.format(x / 1e9)
        else:
            return '${:+.2g}$'.format((x - center) / 1e9)

    return _CenteredGigaHertzFormatter


def plot_xsec(lookups, pressure=100, ax=None):
    """Plots absorption cross section from lookup tables."""
    if ax is None:
        ax = plt.gca()

    waterline = 183.31010705357e9

    zorders = {'Jupiter': -3, 'Earth': -2, 'Mars': -1, 'Venus': 0}

    for planet, abs_lookup in lookups.items():
        selected_pressure = numpy.isclose(abs_lookup.pressuregrid, pressure)
        species_h2o = numpy.where(
            [s == ['H2O-*-*-*'] for s in abs_lookup.speciestags])

        xsec_h2o = abs_lookup.absorptioncrosssection[
                   species_h2o, 0, :, selected_pressure].flatten()

        ax.plot(abs_lookup.frequencygrid, xsec_h2o, linewidth=3,
                label=planet, rasterized=True, zorder=zorders[planet])

    # Sort legend entries by their cross section peak values
    ax.legend(*zip(*((h, l) for _, h, l in
                     sorted(
                         zip([numpy.max(h.get_ydata()) for h in
                              ax.get_legend_handles_labels()[0]],
                             *ax.get_legend_handles_labels()),
                         reverse=True))), frameon=False, loc=(0.7, 0.4))

    ax.xaxis.set_major_formatter(CenteredGigaHertzFormatter(waterline))
    ax.set_xlim((waterline - 0.3e9, waterline + 0.3e9))
    ax.xaxis.set_ticks(
        [x * 1e9 + waterline for x in numpy.arange(-0.2, 0.21, 0.1)])
    ax.set_ylim((0., ax.get_ylim()[1]))
    ax.plot([waterline, waterline], [0, ax.get_ylim()[1]],
            zorder=-10, color='k', linewidth=0.8)
    # ax.text(waterline, ax.get_ylim()[1],
    #         f'${waterline/1e9:.2f}' + r'\,\mathrm{GHz}$',
    #         horizontalalignment='center',
    #         fontsize='x-small')
    ax.set_xlabel(r'Frequency $\left[\mathrm{GHz}\right]$')
    ax.set_ylabel(r'Absorption cross section $\left[\mathrm{m}^2\right]$')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def do_arts_common_setup(ws, pressures):
    """Initializes common settings for all planets."""
    # ws.verbositySetScreen(level=0)
    ws.execute_controlfile('general.arts')
    ws.execute_controlfile('continua.arts')
    ws.execute_controlfile('agendas.arts')
    ws.jacobianOff()
    ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
    ws.IndexSet(ws.stokes_dim, 1)
    # Load 183 GHz water vapor line data
    ws.abs_speciesSet(species=['H2O'])
    ws.abs_linesReadFromArts(filename='input/Perrin_H2O_183.xml',
                             fmin=0., fmax=1e12)
    ws.abs_lines_per_speciesCreateFromLines()
    ws.AtmosphereSet1D()
    ws.VectorNLinSpace(ws.f_grid, 10000, 182910107053.57, 183710107053.57)
    ws.VectorSet(ws.p_grid, pressures)
    ws.sensorOff()
    ws.StringSet(ws.iy_unit, '1')


def calc_arts_lookup_table_for_planet(ws, include, species, basename):
    """Planet specific setup and lookup table calculation."""
    ws.execute_controlfile(include)
    # Enable all available broadening species
    ws.abs_speciesSet(species=species)
    ws.AtmRawRead(basename=basename)
    ws.abs_lines_per_speciesCreateFromLines()
    ws.AtmFieldsCalc()
    ws.abs_xsec_agenda_checkedCalc()
    ws.atmfields_checkedCalc()
    ws.abs_lookupSetup()
    ws.abs_lookupCalc()


def str_vmr_table(vmrs, colwidth='12s'):
    all_species = sorted(set(
        [y for pl in PLANET_SETUP.values() for y in pl['species']]))

    strtable = []
    strtable.append(format('Planet', colwidth)
                    + ''.join(format(s, colwidth) for s in all_species))
    for planet in vmrs:
        s = [planet]
        for species in all_species:
            if species in PLANET_SETUP[planet]['species']:
                s.append(format(vmrs[planet][numpy.array(
                    PLANET_SETUP[planet]['species']) == species][0], '4.2e'))
            else:
                s.append('n/a')
        strtable.append(''.join(format(x, colwidth) for x in s))
    return strtable


def main():
    """Main program."""
    outdir = sys.argv[1] if len(sys.argv) == 2 else '.'

    plt.rc('text', usetex='PLT_NOTEX' not in os.environ)
    matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}',
                                                  r'\sansmath']
    plt.style.use(styles('typhon'))

    abs_lookups = {}

    ws = Workspace(verbosity=0)

    wanted_pressure = 700
    # We are only interested in the wanted Pa level, but need more to satisfy
    # the abslookup interpolation order
    pressures = numpy.array((700, 600, 500, 400, 300, 200, 100))

    print('Performing ARTS calculation')
    do_arts_common_setup(ws, pressures)

    vmrs = {}
    for planet, item in PLANET_SETUP.items():
        calc_arts_lookup_table_for_planet(ws, **item)
        abs_lookups[planet] = ws.abs_lookup.to_typhon()
        vmrs[planet] = ws.vmr_field.to_typhon()[:,
                       numpy.isclose(pressures, wanted_pressure), 0,
                       0].flatten()

    print('VMRs:')
    for s in str_vmr_table(vmrs):
        print(s)

    os.makedirs(outdir, exist_ok=True)

    print('Plotting')
    fig, ax = plt.subplots()
    plot_xsec(abs_lookups, wanted_pressure, ax=ax)
    fig.savefig(os.path.join(outdir, f'xsec_{wanted_pressure:03.0f}.pdf'),
                dpi=300)

    if 'PLT_SHOW' in os.environ:
        plt.show()


if __name__ == '__main__':
    main()
