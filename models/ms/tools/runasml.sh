#!/bin/sh
# This script is an example showing how to run a sequence of assimilation with
# varying parameters.
# In this example we run the SEIK filter with ensemble size 30. We perform a
# set of experiments with varying values for the forgetting factor.
#
# 2010-02, Lars Nerger, AWI
# $Id$

# Name of executable
EXE="./pdaf_mitsch"

# General settings for all experiments
DEFAULTS="-total_steps 50000 -step_null 0 -dim_ens 30"

FILTER=7

# Run experiments
for FORGET in 1
do
    $EXE $DEFAULTS -filtertype $FILTER -forget $FORGET \
	-file_asml t${FILTER}_N30_f${FORGET}.nc -forcing 0.01
done
