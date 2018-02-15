#!/bin/bash

ARTS=$HOME/Hacking/arts-build/clang-reldebug/src/arts

mkdir -p out-refract
$ARTS -I$HOME/Dropbox/Hacking/sat/arts/controlfiles \
    -D$HOME/Dropbox/Hacking/sat/arts-xml-data \
    -r022 -o out-refract \
    TestRefractPlanets.arts
