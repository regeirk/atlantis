#!/bin/bash
nproc=$1
fpath="."
fname="dump_$nproc.xyz"
#region="-R20.125/380.125/-89.875/90.125"
region="-R20.125/379.875/-89.875/89.875"
resolution="-I0.25"


# Updates the path to run all GMT scripts
export NETCDFHOME=/usr/lib
export GMTHOME=/usr/lib/gmt
export PATH=$PATH:$GMTHOME/bin

blockmedian $region $resolution $fpath/$fname > block_$nproc.xyz
surface block_$nproc.xyz $region $resolution -Goutput_$nproc.grd -T0.25 -C0.1 -Vl
