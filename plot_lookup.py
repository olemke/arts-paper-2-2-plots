#!/usr/bin/env python3
"""Visualize an absorption lookup table.

Usage: python plot_lookup.py -h

Author: oliver.lemke@uni-hamburg.de
"""
import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy
import typhon
from matplotlib.ticker import FuncFormatter

PLANETS_Z = {
    'Earth:': 'input/Earth.tropical.z.xml',
    'Mars:': 'input/Mars.Ls0.day.dust-medium.sol-avg.z.xml',
    'Venus:': 'input/Venus.vira.day.z.xml',
    'Jupiter:': 'input/Jupiter.mean.z.xml',
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
    parser.add_argument('-t', '--notex', action='store_true',
                        help='don\'t use TeX rendering')
    parser.add_argument('-s', '--show', action='store_true',
                        help='display plots instead of just storing them')
    parser.add_argument('--altitudes', '-z', metavar='ALTITUDEPROFILE',
                        type=str, help='ARTS altitude profile XML file',
                        required=True)
    return parser.parse_args()


def GigaHertzFormatter():
    @FuncFormatter
    def _GigaHertzFormatter(x, pos):
        return r'${:.0f}$'.format(x / 1e9)

    return _GigaHertzFormatter


def plot_abs_lookup(lookup, opacity, z, species=None, ax=None):
    if ax is None:
        ax = plt.gca()

    if species is None:
        species = lookup.speciestags

    for tag in species:
        xsec = lookup.absorptioncrosssection[0,
               lookup.speciestags.index(tag), :, :]

        alpha = opacity[lookup.speciestags.index(tag), :, :]
        z=numpy.interp(lookup.pressuregrid[::-1], z.grids[0][::-1], z.data[:, 0, 0][::-1])[::-1]

        ax.plot(lookup.frequencygrid,
                #numpy.abs(numpy.trapz(alpha, lookup.pressuregrid, axis=1)),
                numpy.abs(numpy.trapz(alpha, z, axis=1)),
                label=',\n'.join(tag),
                rasterized=True)
    ax.legend(fontsize='x-small', frameon=False)
    ax.xaxis.set_major_formatter(GigaHertzFormatter())
    ax.xaxis.set_ticks([100e9, 200e9, 231e9, 300e9])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def calc_opacity_from_lookup(lookup):
    ni = (lookup.pressuregrid * lookup.referencevmrprofiles
          / lookup.referencetemperatureprofile / typhon.constants.boltzmann
          ).reshape(len(lookup.speciestags), 1, len(lookup.pressuregrid))

    return ni * lookup.absorptioncrosssection[0, :, :, :]


def main():
    args = parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    plt.rc('text', usetex=not args.notex)
    matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}',
                                                  r'\sansmath']
    plt.style.use(typhon.plots.styles('typhon'))

    lookup = typhon.arts.xml.load(args.lookup)
    opacity = calc_opacity_from_lookup(lookup)
    cols = 3
    rows = int(numpy.ceil(len(lookup.speciestags) / cols))
    fig, ax = plt.subplots(rows, cols, figsize=(10, rows * 2))

    fig.tight_layout()
    # plot_abs_lookup(lookup, ['PH3-*-*-*'])
    for cax, species in zip(ax.flatten(), lookup.speciestags):
        plot_abs_lookup(lookup, opacity,
                        typhon.arts.xml.load(args.altitudes),
                        species=[species], ax=cax)

    if rows * cols > len(lookup.speciestags):
        for cax in ax.flatten()[len(lookup.speciestags):]:
            cax.axis('off')

    filename = os.path.join(outdir, os.path.splitext(
        os.path.basename(args.lookup))[0] + '_opacity.pdf')
    print(f'Saving {filename}')
    fig.savefig(filename, dpi=300)

    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
