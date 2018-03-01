#!/usr/bin/env python3
"""Create and plot absorption cross sections for different planets.

Plots the 183 GHz H2O line for different planets at 100 Pa.

Usage: python plot_xsec.py -h

The following environment variables must be set:

export ARTS_INCLUDE_PATH=.:/home/oliver/arts/controlfiles/general:/home/oliver/arts-xml-data
export ARTS_BUILD_PATH=/home/oliver/arts/build

Author: oliver.lemke@uni-hamburg.de
"""
import argparse
import os

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
    ax.xaxis.set_ticks(
        [x * 1e9 + waterline for x in numpy.arange(-0.2, 0.21, 0.1)])
    ax.set_xlim((waterline - 0.3e9, waterline + 0.3e9))
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


def arts_common_setup(ws, pressures):
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
    ws.VectorNLinSpace(ws.f_grid, 1000, 182910107053.57, 183710107053.57)
    ws.VectorSet(ws.p_grid, pressures)
    ws.sensorOff()
    ws.StringSet(ws.iy_unit, '1')


def arts_calc_atmfields(ws, include, species, basename):
    """Planet specific setup of atmospheric fields."""
    ws.execute_controlfile(include)
    # Enable all available broadening species
    ws.abs_speciesSet(species=species)
    ws.AtmRawRead(basename=basename)
    ws.abs_lines_per_speciesCreateFromLines()
    ws.AtmFieldsCalc()


def arts_calc_lookup_table(ws):
    """Lookup table calculation."""
    ws.abs_xsec_agenda_checkedCalc()
    ws.atmfields_checkedCalc()
    ws.abs_lookupSetup()
    ws.abs_lookupCalc()


def str_table(vmrs, temps, colwidth='12s'):
    """Create ASCII table with VMR and temperature values."""
    all_species = sorted(set(
        [y for pl in PLANET_SETUP.values() for y in pl['species']]))

    strtable = []
    strtable.append(format('Planet', colwidth) + format('T', colwidth)
                    + ''.join(format(s, colwidth) for s in all_species))
    for planet in vmrs:
        s = [planet, format(temps[planet][0], '3.1f')]
        for species in all_species:
            if species in PLANET_SETUP[planet]['species']:
                s.append(format(vmrs[planet][numpy.array(
                    PLANET_SETUP[planet]['species']) == species][0], '4.2e'))
            else:
                s.append('n/a')
        strtable.append(''.join(format(x, colwidth) for x in s))
    return strtable


def parse_args():
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description='Create and plot absorption cross sections for '
                    'different planets.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outdir', metavar='OUTDIR', type=str, nargs='?',
                        default='.', help='output directory')
    parser.add_argument('-t', '--notex', action='store_true',
                        help='don\'t use TeX rendering')
    parser.add_argument('-p', '--pressure', type=float, default=700,
                        help='pressure used for plot')
    parser.add_argument('-s', '--show', action='store_true',
                        help='display plots instead of just storing them')
    parser.add_argument('-e', '--tearth', action='store_true',
                        help='use Earth\'s temperature for all planets')
    return parser.parse_args()


def main():
    """Main program."""
    args = parse_args()
    USE_EARTH_TFIELD = args.tearth
    outdir = args.outdir
    wanted_pressure = args.pressure

    plt.rc('text', usetex=not args.notex)
    matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}',
                                                  r'\sansmath']
    plt.style.use(styles('typhon'))

    abs_lookups = {}

    ws = Workspace(verbosity=0)

    # We are only interested in the wanted Pa level, but need more to satisfy
    # the abslookup interpolation order
    pressures = numpy.linspace(wanted_pressure, wanted_pressure - 10, num=10,
                               endpoint=True)

    print('Performing ARTS calculation')
    arts_common_setup(ws, pressures)

    vmrs = {}
    temps = {}
    if USE_EARTH_TFIELD:
        arts_calc_atmfields(ws, **PLANET_SETUP['Earth'])
        t_field_earth = ws.t_field.to_typhon()

    for planet, item in PLANET_SETUP.items():
        arts_calc_atmfields(ws, **item)
        if USE_EARTH_TFIELD:
            ws.t_field = t_field_earth
        arts_calc_lookup_table(ws)
        abs_lookups[planet] = ws.abs_lookup.to_typhon()
        vmrs[planet] = ws.vmr_field.to_typhon()[:,
                       numpy.isclose(pressures, wanted_pressure), 0,
                       0].flatten()
        temps[planet] = ws.t_field.to_typhon()[
            numpy.isclose(pressures, wanted_pressure), 0,
            0].flatten()

    for s in str_table(vmrs, temps):
        print(s)

    os.makedirs(outdir, exist_ok=True)

    print('Plotting')
    fig, ax = plt.subplots()
    plot_xsec(abs_lookups, wanted_pressure, ax=ax)
    tearth = temps['Earth'][0]
    filename = os.path.join(outdir,
                            f'xsec_{wanted_pressure:03.0f}Pa'
                            + (f'_Tearth_{tearth:.0f}K'
                               if USE_EARTH_TFIELD else '')
                            + '.pdf')
    print(f'Saving {filename}')
    fig.savefig(filename, dpi=300)

    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
