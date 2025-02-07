import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import argparse as ap

if __name__ == "__main__":
	#load state.nc
	file = nc.Dataset('../t1_N30_f1.nc', 'r')

	# print details of state
	print(file)

	# print the variables
	print(file.variables)

	# print the forecast ensemble
	print(file['state_for'][:].data[:])
