#!/bin/bash

SCRIPTDIR=`dirname $0`
case $1 in
    mutect) shift; exec python $SCRIPTDIR/filter_mutect.py "$@" ;;
    vardict) shift; exec python $SCRIPTDIR/filter_vardict.py "$@" ;;
    complex) shift; exec python $SCRIPTDIR/filter_complex.py "$@" ;;
    *) echo "mutect, vardict, or complex?"; exit 1 ;;
esac
