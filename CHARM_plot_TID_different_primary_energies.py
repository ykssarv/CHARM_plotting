#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Program for managing FLUKA's USRTRACK scores of the CHARM facility."""

import matplotlib.pyplot as plt
from numpy import asarray

colors = ["navy", "orange", "gold", "blue", "red", "darkblue", "darkred", "darkgreen", "darkorange", "darkviolet", "magenta", "cyan", "crimson"]
markers = ["o", "P", "*", "<", ">", "v", "^", "1", "2", "3", "4", "x", "s"]


vol = 20.0 * 20.0 * 20.0  # Volume of the scoring regions
norm_energy = 1.602176462E-7  # Normalization from GeV/g to Gy


def read_data_frame(directory, input_filename):
	filename = f'{directory}/{input_filename}_usrbin_66'
	print(f'Using file: {filename} for TID values')
	with open(filename) as file:
		data = file.readlines()

	data_values = data[10].split("\s+")
	data_values_numeric = asarray([float(value[:-1]) for value in data_values])
	data_values_numeric *= norm_energy / vol
	# print(data_values_numeric)
	
	data_errors = data[14].split("\s+")
	data_errors_numeric_relative = asarray([float(value[:-1]) for value in data_errors])
	data_errors_numeric_absolute = data_values_numeric * data_errors_numeric_relative / 100 
	# print(data_errors_numeric)

	return data_values_numeric, data_errors_numeric_absolute


REGION = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10', 'R11', 'R12', 'R13',
		  'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6',
		  'G0', 'G0old',
		  'PC0', 'PCint',
		  'T0u', 'T0', 'T0d',
		  # 'F0u', 'F0', 'F0d'
          ]


def main(directory, output_ed_name, input_filename, energies):
	# Activate Plotting of the spectra for a given position
	saveplot = True
	fileEXT = 'png'  # pdf or jpeg or png

	# Create the plot
	wplot = 9.0 / 2.54 
	width = wplot * 3 / 2
	height = wplot
	fsize = 9.5
	msize1 = 2
	msize2 = 3
	mewidth = 0.2
	llwidth = 0.5
	legsize = 14
	lwidth = 0.15

	# Set the test location for plotting.
	for i in range(len(REGION)):
		fig1 = plt.figure(num=1, figsize=(width, height))
		sub1 = fig1.add_subplot(111, aspect='auto')

		for energy in energies:
			directory_analysis = f"{directory}/{output_ed_name}{energy}/analysis"
			filename = f"{input_filename}{energy}"
			# plot_lethargy_spectra_legacy(directory, target, config, rack, saveplot, fileEXT)
			data_locations, err_locations = read_data_frame(directory_analysis, filename)

			rack = REGION[i]

			value_at_region = data_locations[i]
			error_at_region = err_locations[i]

			sub1.errorbar(energy, value_at_region, yerr=error_at_region,
			              linestyle='', linewidth=llwidth, color=colors[i],
			              marker=markers[1], markeredgecolor='black', markersize=msize2, markeredgewidth=mewidth,
			              drawstyle='steps-mid', label=f'{energy} GeV')

		# Create the labels
		# sub1.axis([1e-14,1e2, 1e-2, 1e-13], fontsize=fsize)
		axis_font = {'fontname': 'Times New Roman', 'size': 18, 'weight': 'normal'}
		# sub1.axis([1e-13, 1e2, 1e-12, 1e-2])  # fontsize=fsize
		#plt.yscale("log")
		sub1.tick_params(labelsize=fsize)
		plt.xlabel('Primary Beam Energy [GeV]', fontsize=fsize)
		plt.ylabel('TID [Gy/POT]', fontsize=fsize)
		# plt.title('CHARM - Configuration {0:<8s} - Test Location {1:<5s}'.format(config, rack), fontstyle='italic',fontsize=fsize)

		sub1.set_xlim(xmax=24)

		# Create Grid
		sub1.grid()
		sub1.grid(b=True, axis="both", which='both', linestyle='-', linewidth=lwidth)
		# sub1.grid(b=True, which='minor', linestyle='-', linewidth=lwidth)

		# Create Legend
		from matplotlib import rcParams
		rcParams['font.family'] = 'sans-serif'
		rcParams['font.sans-serif'] = ['Times New Roman']
		# rcParams['font.size'] = 20
		# rcParams['axes.linewidth'] = 0.1 #set the value globally
		for axis in ['top', 'bottom', 'left', 'right']:
			sub1.spines[axis].set_linewidth(0.5)
		sub1.xaxis.set_tick_params(width=0.5)
		sub1.yaxis.set_tick_params(width=0.5)
		# leg1=sub1.legend(loc='lower left', numpoints=1, prop={'size':12,'family':'sans-serif'}, shadow=True)
		# leg1.get_frame().set_linewidth(0.0)
		#plt.legend(fontsize=fsize, loc=7)

		# sub1.annotate('(b)', xy=(0.90,0.90), xycoords='axes fraction', fontweight='bold', fontsize=9)

		# Make/Save the plot
		plt.tight_layout(pad=3, w_pad=1.0, h_pad=2)

		if saveplot:
			outfile = f"plots/TID/CHARM_{rack}_TID.{fileEXT}"
			plt.savefig(outfile, dpi=250, bbox_inches='tight')
		# plt.savefig(outfile)
		else:
			plt.show()


if __name__ == "__main__":
	# Absolute path to the data files

	directory = '/scratch/charmflk/2022-summer-student-project/simulations'

	output_ed_name = "output_cu_"
	input_filename = 'charm_new_cu'

	energies = [1, 9, 12, 15]

	main(directory, output_ed_name, input_filename, energies)
