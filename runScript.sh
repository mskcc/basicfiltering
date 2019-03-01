#!/bin/bash

SCRIPTDIR=`dirname $0`
case $1 in
    pindel) shift; exec python $SCRIPTDIR/filter_pindel.py "$@" ;;
    mutect) shift; exec python $SCRIPTDIR/filter_mutect.py "$@" ;;
    vardict) shift; exec python $SCRIPTDIR/filter_vardict.py "$@" ;;
    complex) shift; exec python $SCRIPTDIR/filter_complex.py "$@" ;;
    *) echo "pindel, mutect, vardict, or complex?"; exit 1 ;;
esac
