#!/usr/bin/env python3
"""Visualize an absorption lookup table.

Usage: python plot_lookup.py -h

Author: oliver.lemke@uni-hamburg.de
"""
import argparse
import os
import re

import matplotlib
import matplotlib.pyplot as plt
import numpy
import typhon
from cycler import cycler
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import interp1d

PLANETS_G = {
    'Earth': typhon.constants.g,
    'Mars': 3.711,
    'Venus': 8.87,
    'Jupiter': 24.79,
}

PLANETS_R = {
    'Earth': typhon.constants.gas_constant / (
            0.78 * 28.01 + 0.21 * 32 + 0.01 * 39.95) * 1000,
    'Mars': typhon.constants.gas_constant / (
            0.95 * 44.01 + 0.027 * 28.01) * 1000,
    'Venus': typhon.constants.gas_constant / (
            0.97 * 44.01 + 0.044 * 28.01) * 1000,
    'Jupiter': typhon.constants.gas_constant / (
            0.86 * 2.016 + 0.14 * 4.003) * 1000
}


def parse_args():
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description='Visualize an absorption lookup table.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('lookup', metavar='LOOKUPTABLEFILE', type=str,
                        help='ARTS lookup table XML file')
    parser.add_argument('-o', '--outdir', metavar='DIRECTORY', type=str,
                        default='.', help='output directory')
    parser.add_argument('-e', '--temppert', metavar='TEMPPERT', type=int,
                        default=0, help='index of temperature pertubation')
    parser.add_argument('-t', '--notex', action='store_true',
                        help='don\'t use TeX rendering')
    parser.add_argument('--opacity', action='store_true',
                        default=False, help='plot opacity')
    parser.add_argument('-s', '--show', action='store_true',
                        help='display plots instead of just storing them')
    parser.add_argument('-T', '--title', metavar='TITLE', type=str,
                        help='plot title')
    parser.add_argument('-P', '--planet', metavar='PLANET', type=str,
                        default='Earth', help='planet, determines gravity and '
                                              'gas constant')
    parser.add_argument('-p', '--pressures', metavar='PYTHONLIST', type=str,
                        help='list of pressures in Pa to plot, '
                             'e.g. [1,10,100,1000,10000].\n'
                             'Only used for cross-section plots.')
    parser.add_argument('-V', '--vmrpert', metavar='VMRPERT', type=int,
                        default=0, help='index of nonlinear species vmr'
                                        'pertubation')
    parser.add_argument('--altitudes', '-z', metavar='ALTITUDEPROFILE',
                        type=str, help='ARTS altitude profile XML file')
    return parser.parse_args()


def GigaHertzFormatter():
    @FuncFormatter
    def _GigaHertzFormatter(x, pos):
        return r'${:.0f}$'.format(x / 1e9)

    return _GigaHertzFormatter


def _calc_lookup_species_count(lookup):
    nlsspecies = lookup.nonlinearspecies
    speciescount = numpy.ones_like(lookup.speciestags, dtype=int)
    if nlsspecies is not None:
        speciescount[nlsspecies] = lookup.nonlinearspeciesvmrpertubations.size
    return speciescount


def _get_lookup_species_index(lookup, species, vmrpert):
    ret = 0
    spindex = lookup.speciestags.index(species)
    nlsspecies = lookup.nonlinearspecies
    speciescount = _calc_lookup_species_count(lookup)
    if nlsspecies is not None and spindex in nlsspecies:
        if vmrpert >= speciescount[spindex]:
            raise RuntimeError(
                'Nonlinear species VMR pertubation index too large')
        ret = vmrpert

    return ret + (numpy.sum(speciescount[0:spindex]) if spindex > 0 else 0)


def plot_lookup_xsec(lookup, ipressures=None, species=None, ax=None, tpert=0,
                     vmrpert=0):
    if ax is None:
        ax = plt.gca()

    ax.set_yscale('log')
    if species is None:
        species = lookup.speciestags

    for tag in species:
        ax.set_prop_cycle(
            cycler('color', [plt.cm.viridis(i) for i in
                             numpy.linspace(0, 1, len(ipressures))]))
        for pi in ipressures:
            xsec = lookup.absorptioncrosssection[
                   tpert,
                   _get_lookup_species_index(lookup, tag, vmrpert), :, pi]
            ax.plot(lookup.frequencygrid, xsec)

    if len(species) > 1:
        ax.legend(fontsize='xx-small', frameon=False)
    else:
        ax.set_title(
            ',\n'.join(re.sub('(-\*)+$', '', s) for s in species[0]),
            y=1. - len(species[0]) * 0.05,
            fontsize='xx-small')

    ax.xaxis.set_major_formatter(GigaHertzFormatter())
    ax.tick_params(axis='both', which='major', labelsize='xx-small')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def plot_lookup_opacity(lookup, opacity, species=None, vmrpert=0, ax=None,
                        oneline=False, total=False):
    if ax is None:
        ax = plt.gca()

    ax.set_yscale('log')
    if species is None:
        species = lookup.speciestags

    for tag in species:
        ax.plot(lookup.frequencygrid,
                opacity[_get_lookup_species_index(lookup, tag, vmrpert), :],
                label=',\n'.join(tag))
    if oneline:
        ax.plot(lookup.frequencygrid, numpy.ones_like(lookup.frequencygrid),
                linewidth=1, linestyle='--', color='k')
    if total:
        if lookup.nonlinearspecies is not None:
            speciescount = _calc_lookup_species_count(lookup)
            spindex = numpy.cumsum(speciescount)
            spindex[1:] = spindex[0:-1]
            spindex[0] = 0
            spindex[lookup.nonlinearspecies] += vmrpert
            o = opacity[spindex]
        else:
            o = opacity
        ax.plot(lookup.frequencygrid, numpy.sum(o, axis=0),
                linewidth=1, color='k')

    if len(species) > 1:
        ax.legend(fontsize='xx-small', frameon=False)
    else:
        ax.set_title(',\n'.join(re.sub('(-\*)+$', '', s) for s in species[0]),
                     y=1. - len(species[0]) * 0.05,
                     fontsize='xx-small')

    ax.xaxis.set_major_formatter(GigaHertzFormatter())
    ax.tick_params(axis='both', which='major', labelsize='xx-small')
    ax.tick_params(axis='both', which='minor', labelsize='xx-small')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def calc_opacity_from_lookup(lookup, z=None, g=typhon.constants.g,
                             r=typhon.constants.gas_constant_dry_air, tpert=0):
    speciescount = _calc_lookup_species_count(lookup)
    vmrs = (numpy.repeat(lookup.referencevmrprofiles, speciescount, axis=0)
            if lookup.nonlinearspecies is not None
            else lookup.referencevmrprofiles)

    ni = (lookup.pressuregrid * vmrs
          / lookup.referencetemperatureprofile / typhon.constants.boltzmann
          ).reshape(sum(speciescount), 1, len(lookup.pressuregrid))

    alpha = ni * lookup.absorptioncrosssection[tpert, :, :, :]

    if z is not None:
        z = interp1d(z.grids[0], z.data[:, 0, 0])(lookup.pressuregrid)
    else:
        # Calculate z from hypsometric formula
        pgrid = lookup.pressuregrid
        z = [r * t / g * numpy.log(p1 / p2)
             for p1, p2, t in zip(pgrid[:-1], pgrid[1:], (
                    lookup.referencetemperatureprofile[
                    1:] + lookup.referencetemperatureprofile[:-1]) / 2.)]
        z = numpy.cumsum(z)
        p = (pgrid[1:] + pgrid[:-1]) / 2.
        z = interp1d(p, z, fill_value='extrapolate')(lookup.pressuregrid)

    return numpy.vstack(numpy.trapz(ialpha, z, axis=1) for ialpha in alpha)


def add_opacity_legend(ax=None):
    if ax is None:
        ax = plt.gca()

    blue_line = matplotlib.lines.Line2D([], [], label='species opacity')
    black_line = matplotlib.lines.Line2D([], [], color='k', linewidth=1.,
                                         label='total opacity')
    dashed_line = matplotlib.lines.Line2D([], [], color='k', linestyle='--',
                                          linewidth=1., label='opacity=1')

    handles = [blue_line, black_line, dashed_line]
    labels = [h.get_label() for h in handles]

    ax.legend(handles=handles, labels=labels, fontsize='xx-small',
              loc='lower left')


def add_xsec_legend(lookup, ipressures, ax=None):
    if ax is None:
        ax = plt.gca()

    pgrid = lookup.pressuregrid
    colors = [plt.cm.viridis(i) for i in numpy.linspace(0, 1, len(ipressures))]
    handles = [matplotlib.lines.Line2D([], [],
                                       color=colors[i],
                                       label=f'{pgrid[ip]/100.:8.3f} hPa')
               for i, ip in enumerate(ipressures)]

    labels = [h.get_label() for h in handles]

    ax.legend(handles=handles, labels=labels, fontsize='xx-small',
              loc='lower left')


def main():
    args = parse_args()

    outdir = args.outdir
    do_opacity = args.opacity
    os.makedirs(outdir, exist_ok=True)

    plt.rc('text', usetex=not args.notex)
    matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}',
                                                  r'\sansmath']
    plt.style.use(typhon.plots.styles('typhon'))

    lookup = typhon.arts.xml.load(args.lookup)
    if do_opacity:
        opacity = calc_opacity_from_lookup(lookup,
                                           typhon.arts.xml.load(args.altitudes)
                                           if args.altitudes is not None
                                           else None,
                                           g=PLANETS_G[args.planet],
                                           r=PLANETS_R[args.planet],
                                           tpert=args.temppert)
    cols = 3
    rows = int(numpy.ceil(len(lookup.speciestags) / cols))
    fig, ax = plt.subplots(rows + 1, cols, figsize=(10, (rows + 1) * 2))

    fig.tight_layout()
    # plot_abs_lookup(lookup, ['PH3-*-*-*'])
    for cax, species in zip(ax[1:, :].flatten(), lookup.speciestags):
        if do_opacity:
            plot_lookup_opacity(lookup, opacity, vmrpert=args.vmrpert,
                                oneline=True, total=True, species=[species],
                                ax=cax)
        else:
            psize = lookup.pressuregrid.size
            if args.pressures is not None:
                ipressures = [numpy.abs(lookup.pressuregrid - p).argmin() for p
                              in eval(args.pressures)]
            else:
                ipressures = (lookup.pressuregrid.size - 1 - (
                    range(psize) if psize <= 5
                    else numpy.linspace(0, lookup.pressuregrid.size,
                                        num=6,
                                        endpoint=False,
                                        dtype=int)))
            plot_lookup_xsec(lookup, ipressures, species=[species], ax=cax,
                             tpert=args.temppert, vmrpert=args.vmrpert)

    if rows * cols > len(lookup.speciestags):
        for cax in ax[1:, :].flatten()[len(lookup.speciestags):]:
            cax.axis('off')

    for cax in ax[0, :]:
        cax.axis('off')

    if do_opacity:
        add_opacity_legend(ax[0, 0])
    else:
        add_xsec_legend(lookup, ipressures, ax[0, 0])

    for cax in ax[-1, :]:
        cax.set_xlabel('Frequency [GHz]', fontsize='xx-small')

    filename = os.path.join(outdir, os.path.splitext(
        os.path.basename(args.lookup))[0])
    if do_opacity:
        filename += '_opacity.pdf'
    else:
        filename += '_xsecs.pdf'

    if args.title is not None:
        fig.suptitle(args.title)
    else:
        fig.suptitle(args.planet)

    print(f'Saving {filename}')
    fig.savefig(filename, bbox='tight', dpi=300)

    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
