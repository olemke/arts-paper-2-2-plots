#!/usr/bin/env python3
"""Create and plot absorption spectra for different planets.

Usage: python plot_spectra.py -h

The following environment variables must be set:

export ARTS_INCLUDE_PATH=.:$HOME/arts/controlfiles/general:$HOME/arts-xml-data
export ARTS_BUILD_PATH=$HOME/arts/build

Author: oliver.lemke@uni-hamburg.de
"""
import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy
import typhon
from matplotlib.ticker import FuncFormatter
from typhon.arts.workspace import Workspace, arts_agenda

NPRESSURES = 100

PLANET_SETUP = {
    'Earth': {
        'species': [
            'C2H2', 'C2H6', 'CH3Cl', 'CH4', 'CO', 'CO2', 'CO2-CIA-CO2-0',
            'COF2', 'ClO', 'H2CO',
            'H2O-SelfContStandardType, H2O-ForeignContStandardType, H2O',
            'H2O2', 'HBr', 'HCN', 'HCl', 'HF', 'HI', 'HNO3', 'HOCl', 'N2',
            'N2O', 'NH3', 'NO', 'NO2', 'O2', 'O3', 'OH', 'PH3', 'SF6', 'SO2'],
        'basename': 'planets/Earth/Fascod/tropical/tropical',
        'include': 'planet_earth.arts',
        'p_grid': (NPRESSURES * 2, 1100e2, 0.1),
        'reflectivity': 0.,
    },
    'Mars': {
        'species': [
            'CH4', 'CO', 'CO2', 'CO2-CIA-CO2-0', 'H2', 'H2O', 'H2O2', 'H2S',
            'HCl', 'N2', 'NO2', 'O', 'O2', 'O3', 'OCS', 'OH', 'SO2'],
        'basename': 'planets/Mars/MPS/Mars.Ls0.day.dust-medium/'
                    'Mars.Ls0.day.dust-medium.sol-avg/'
                    'Mars.Ls0.day.dust-medium.sol-avg',
        'include': 'planet_mars.arts',
        'p_grid': (NPRESSURES, 766., 0.1),
        'reflectivity': 0.13,
    },
    'Venus': {
        'species': [
            'CO', 'CO2', 'CO2-CIA-CO2-0', 'H2O', 'H2SO4', 'HCl', 'HF', 'N2',
            'NO', 'NO2', 'O', 'O2', 'O3', 'OCS', 'SO', 'SO2'],
        'basename': 'planets/Venus/MPS/Venus.vira.day/Venus.vira.day',
        'include': 'planet_venus.arts',
        'p_grid': (NPRESSURES * 6, 9.2e6, 0.1),
        'reflectivity': 0.15,
    },
    'Jupiter': {
        'species': [
            'C2H2', 'C2H6', 'C3H8', 'CH4', 'CO', 'CO2', 'H2', 'H2O',
            'H2S', 'HCN', 'He', 'NH3', 'PH3'],
        'basename': 'planets/Jupiter/MPS/Jupiter.mean/Jupiter.mean',
        'include': 'planet_jupiter.arts',
        'p_grid': (NPRESSURES * 6, 1.0e6, 0.1),
        'reflectivity': 0.,
    },
}


def parse_args():
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description='Create and plot absorption spectra for different planets.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outdir', metavar='OUTDIR', type=str, nargs='?',
                        default='.', help='output directory')
    parser.add_argument('--recalc', '-r', action='store_true',
                        help='recalculate lookup tables')
    parser.add_argument('-t', '--notex', action='store_true',
                        help='don\'t use TeX rendering')
    parser.add_argument('-s', '--show', action='store_true',
                        help='display plots instead of just storing them')
    return parser.parse_args()


def GigaHertzFormatter():
    @FuncFormatter
    def _GigaHertzFormatter(x, pos):
        return r'${:.0f}$'.format(x / 1e9)

    return _GigaHertzFormatter


def plot_spectra(y_all, ax=None):
    """Plots absorption cross section from lookup tables."""
    if ax is None:
        ax = plt.gca()

    zorders = {'Jupiter': -3, 'Earth': -2, 'Mars': -1, 'Venus': -4}

    for planet, y_data in y_all.items():
        ax.plot(y_data['f_grid'], y_data['y'],
                label=planet, rasterized=True, zorder=zorders[planet])

    ax.legend(*typhon.plots.sorted_legend_handles_labels(),
              frameon=False)
    ax.xaxis.set_major_formatter(GigaHertzFormatter())
    ax.set_xlabel(r'Frequency $\left[\mathrm{GHz}\right]$')
    ax.set_ylabel(r'Brightness temperature $\left[\mathrm{K}\right]$')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def arts_common_setup(ws):
    """Initializes common settings for all planets."""
    # ws.verbositySetScreen(level=0)
    ws.execute_controlfile('general.arts')
    ws.execute_controlfile('continua.arts')
    ws.execute_controlfile('agendas.arts')
    ws.jacobianOff()
    ws.cloudboxOff()
    # ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
    ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__withCIAextraT)
    ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
    ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
    ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
    ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
    ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
    # ws.Copy(ws.surface_rtprop_agenda,
    # surface_rtprop_agenda__Blackbody_FixT)  # FIXME ok?
    # ws.NumericSet(ws.surface_skin_t, 180.)
    ws.Copy(ws.surface_rtprop_agenda,
            ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface)
    ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
    ws.IndexSet(ws.stokes_dim, 1)
    ws.AtmosphereSet1D()
    ws.VectorNLinSpace(ws.f_grid, 1000, 100e9, 300e9)
    ws.sensorOff()
    # ws.StringSet(ws.iy_unit, '1')
    ws.StringSet(ws.iy_unit, 'PlanckBT')
    ws.sensor_pos = numpy.array([[600e3]])
    ws.sensor_los = numpy.array([[180.]])


def arts_calc_atmfields(ws, include, species, basename, p_grid, reflectivity):
    """Planet specific setup of atmospheric fields."""
    ws.execute_controlfile(include)
    ws.VectorNLinSpace(ws.p_grid, *p_grid)
    # Enable all available broadening species
    ws.abs_speciesSet(species=species)
    ws.AtmRawRead(basename=basename)
    ws.AtmFieldsCalc(vmr_zeropadding=1)
    ws.Extract(ws.z_surface, ws.z_field, 0)
    ws.Extract(ws.t_surface, ws.t_field, 0)

    ws.VectorSetConstant(ws.surface_scalar_reflectivity, 1, reflectivity)

    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()
    ws.atmgeom_checkedCalc()


def arts_calc_lookup_table(ws):
    """Lookup table calculation."""
    ws.abs_xsec_agenda_checkedCalc()
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.abs_lookupSetup()
    ws.abs_lookupCalc()


@arts_agenda
def propmat_clearsky_agenda_lookup_faraday(ws):
    ws.Ignore(ws.rtp_mag)
    ws.Ignore(ws.rtp_los)
    ws.Ignore(ws.rtp_nlte)
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddFromLookup()
    # ws.propmat_clearskyAddFaraday()


@arts_agenda
def surface_rtprop_agenda__Blackbody_FixT(ws):
    ws.Touch(ws.surface_skin_t)
    ws.surfaceBlackbody()


def annotate_lines(ax=None):
    if ax is None:
        ax = plt.gca()

    LINES = (
        ('O$_2$', 119., 208.952),
        ('H$_2$O', 183., 243.323),
        ('CO', 231., 231.678),
        ('PH$_3$', 267., 144.114),
    )

    for label, x, y in LINES:
        ax.text(x * 1e9, y - 9, label, horizontalalignment='center',
                fontsize='xx-small')


def main():
    """Main program."""
    args = parse_args()
    outdir = args.outdir
    recalc_lookups = args.recalc
    os.makedirs(outdir, exist_ok=True)

    plt.rc('text', usetex=not args.notex)
    matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}',
                                                  r'\sansmath']
    plt.style.use(typhon.plots.styles('typhon'))

    ws = Workspace(verbosity=2)
    ws.verbosityInit()

    print('Performing ARTS calculation')
    arts_common_setup(ws)

    y_all = {}

    for planet, item in PLANET_SETUP.items():
        arts_calc_atmfields(ws, **item)

        lookup_file = os.path.join(outdir, planet + '_lookup.xml')
        if os.path.isfile(lookup_file) and not recalc_lookups:
            ws.ReadXML(ws.abs_lookup, filename=lookup_file)
            ws.abs_lookupAdapt()
        else:
            ws.ReadXML(ws.abs_cia_data,
                       filename='spectroscopy/cia/hitran2011/'
                                'hitran_cia2012_adapted.xml.gz')
            ws.abs_linesReadFromSplitArtscat(
                basename='spectroscopy/Perrin/', fmin=0., fmax=1e12)
            ws.abs_lines_per_speciesCreateFromLines()
            arts_calc_lookup_table(ws)
            ws.WriteXML('binary', ws.abs_lookup, lookup_file)

        ws.propmat_clearsky_agenda_checkedCalc()
        ws.yCalc()
        ws.WriteXML("ascii", ws.y,
                    filename=os.path.join(outdir, planet + '.y.xml'))
        y_all[planet] = {
            'f_grid': ws.f_grid.to_typhon(),
            'y': ws.y.to_typhon(),
        }

    print('Plotting')
    fig, ax = plt.subplots()
    plot_spectra(y_all)
    ax.set_xlim(100e9, 300e9)
    annotate_lines()
    filename = os.path.join(outdir, 'planet_spectra.pdf')
    print(f'Saving {filename}')
    fig.savefig(filename, dpi=300)

    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
