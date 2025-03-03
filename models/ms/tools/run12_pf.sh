#!/bin/sh
# This script is an example showing how to run a sequence of assimilation with
# varying parameters.
# In this example we run the SEIK filter with ensemble size 30. We perform a
# set of experiments with varying values for the forgetting factor.
#
# 2010-02, Lars Nerger, AWI
# $Id: runasml.sh 350 2020-01-30 18:35:17Z lnerger $

# Name of executable
EXE="./pdaf_mitsch"

# General settings for all experiments
DEFAULTS="-total_steps 10 -step_null 1000 -dim_ens 30 -type_ensinit 'eof' -pf_noise_amp 1.5 -pf_noise_type 1"

FILTER=12

# Run experiments
for FORGET in 1.0
do
    $EXE $DEFAULTS -filtertype $FILTER -forget $FORGET -cradius 10 \
	-file_asml t${FILTER}_N30_f${FORGET}.nc
done
