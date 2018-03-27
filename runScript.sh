#!/bin/bash

SCRIPTDIR=`dirname $0`
case $1 in
    pindel) shift; exec python $SCRIPTDIR/filter_pindel.py "$@" ;;
    mutect) shift; exec python $SCRIPTDIR/filter_mutect.py "$@" ;;
    vardict) shift; exec python $SCRIPTDIR/filter_vardict.py "$@" ;;
    sid) shift; exec python $SCRIPTDIR/filter_sid.py "$@" ;;
    *) echo "pindel, mutect, vardict, or sid?"; exit 1 ;;
esac
