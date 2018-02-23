#!/bin/bash

ARTS=$ARTS_BUILD_PATH/src/arts

OUTDIR=out-refractivity

mkdir -p $OUTDIR
$ARTS -I$ARTS_INCLUDE_PATH \
    -D$ARTS_XML_DATA_PATH \
    -r022 -o $OUTDIR \
    CalcRefraction.arts

echo -n "Creating figures... "
python plot_refractivity.py $OUTDIR
echo "Done"

