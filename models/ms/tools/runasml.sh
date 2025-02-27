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
DEFAULTS="-total_steps 100000 -step_null 0 -dim_ens 30 -model_err_amp 0.01"

# Run experiments
for FORGET in 1; do
	for FILTER in 6; do #4 5 6 7 9 10 11; do
		for FORCING in 0.01; do
			$EXE $DEFAULTS -filtertype $FILTER -forget $FORGET -forcing ${FORCING} \
				-file_asml FILTER${FILTER}_FORGET${FORGET}_FORCING${FORCING}_ME_0.01.nc
		done
	done
done
