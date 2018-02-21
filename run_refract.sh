#!/bin/bash

ARTS=$HOME/Hacking/arts-build/clang-debug/src/arts

OUTDIR=out-refractivity
mkdir -p $OUTDIR
$ARTS -I$HOME/Dropbox/Hacking/sat/arts/controlfiles \
    -D$HOME/Dropbox/Hacking/sat/arts-xml-data \
    -r022 -o $OUTDIR \
    CalcRefraction.arts
