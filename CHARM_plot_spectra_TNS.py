#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Program for managing FLUKA's USRTRACK scores of the CHARM facility."""

import matplotlib.pyplot as plt
import numpy as np


def extract_value_at_rack(
		rack, energy_average_geometric_neutron, energy_difference_neutron,
		energy_average_geometric_neutron2, energy_difference_neutron2):
	"""Sets the variable for the reading of the test location data."""
	if rack == 'R1':
		loc = 0
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R2':
		loc = 1
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R3':
		loc = 2
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R4':
		loc = 3
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R5':
		loc = 4
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R6':
		loc = 5
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R7':
		loc = 6
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R8':
		loc = 7
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R9':
		loc = 8
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R10':
		loc = 9
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R11':
		loc = 10
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R12':
		loc = 11
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'R13':
		loc = 12
		energy_average_geometric = energy_average_geometric_neutron
		energy_difference = energy_difference_neutron
	elif rack == 'M0':
		loc = 13
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'M1':
		loc = 14
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
		# energy_average_geometric = list(energy_average_geometric_neutron) + list(energy_average_geometric_neutron2)
		# energy_difference = list(energy_difference_neutron) + list(energy_difference_neutron2)
	elif rack == 'M2':
		loc = 15
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'M3':
		loc = 16
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'M4':
		loc = 17
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'M5':
		loc = 18
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'M6':
		loc = 19
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'G0':
		loc = 20
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'G0old':
		loc = 21
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'PC0':
		loc = 22
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'PCint':
		loc = 23
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'T0u':
		loc = 24
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'T0':
		loc = 25
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'T0d':
		loc = 26
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'F0u':
		loc = 27
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'F0':
		loc = 28
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	elif rack == 'F0d':
		loc = 29
		energy_average_geometric = energy_average_geometric_neutron2
		energy_difference = energy_difference_neutron2
	else:
		raise Exception('ERROR: Unkown Test Location')

	return loc, energy_average_geometric, energy_difference


def compute_lethargy_values(value_raw, error_raw, energy_geo, energy_threshold = 1e-4):
	"""Computes lethargy values."""
	for index in range(len(energy_geo)):
		energy_level = energy_geo[index]
		if energy_level < energy_threshold:
			continue
		else:
			break
	value = value_raw[index:] * energy_geo[index:] / vol
	error = value * error_raw[index:] / 100.0
	return value, error, energy_geo[index:]

colors = ["navy", "orange", "gold", "blue", "red", "darkblue", "darkred", "darkgreen", "darkorange", "darkviolet", "magenta", "cyan", "crimson"]
markers = ["o", "P", "*", "<", ">", "v", "^", "1", "2", "3", "4", "x", "s"]


vol = 20.0 * 20.0 * 20.0  # Volume of the scoring regions
norm_energy = 1.602176462E-7  # Normalization from GeV/g to Gy

particle_bins = {
	"proton": 80,
	"neutron": 81,
	"photon": 82,
	"electron": 83,
	"positron": 84,
	"pion_positive": 85,
	"pion_negative": 86,
	"muon_positive": 87,
	"muon_negative": 88,
	"kaon_positive": 89,
	"kaon_negative": 90,
	"kaon_neutral": 91,
	"alpha": 92,
	"hehad": 93,
	"heheq": 94,
	"thneq": 95,
}

energy_thresholds = {
	"proton": 1e-4,
	"neutron": 1e-14,
	"photon": 1e-4,
	"electron": 1e-4,
	"positron": 1e-4,
	"pion_positive": 1e-4,
	"pion_negative": 1e-4,
	"muon_positive": 1e-4,
	"muon_negative": 1e-4,
	"kaon_positive": 1e-4,
	"kaon_negative": 1e-4,
	"kaon_neutral": 1e-4,
	"alpha": 1e-4,
	"hehad": 1e-14,
	"heheq": 1e-14,
	"thneq": 1e-14,
}


def read_data_frame(directory, input_filename, particle):
	filename = f'{directory}/{input_filename}_usrtrk_{particle_bins[particle]}_tab.lis'
	print(f'Using file: {filename} for {particle} spectra')
	df = np.genfromtxt(filename, names=['binmin', 'binmax', 'diffF', 'relerr'])
	return df


def read_energy(df, binno):
	"""Reads and the computes the energy values.
	# Note [0:binno] is equivalent to [0:80] where the last number is not inclusive
	# -> from line 0 to line 79
	"""
	energy_min = df[0:binno]['binmin']
	energy_max = df[0:binno]['binmax']
	energy_delta = (energy_max-energy_min)
	energy_average = (energy_min+energy_max)/2.0
	energy_average_geometric = np.sqrt(energy_min*energy_max)
	return energy_min, energy_max, energy_delta, energy_average, energy_average_geometric


def read_energy_different_binning(df, binno, block, lastblock):
	"""Reads and the computes the energy values for different binning.
	# Note [0:binno] is equivalent to [0:80] where the last number is not inclusive
	# -> from line 0 to line 79
	"""
	energy_min = df[lastblock:lastblock+(binno*block)]['binmin']
	energy_max = df[lastblock:lastblock+(binno*block)]['binmax']
	energy_delta = (energy_max-energy_min)
	energy_average = (energy_min+energy_max)/2.0
	energy_average_geometric = np.sqrt(energy_min*energy_max)
	return energy_min, energy_max, energy_delta, energy_average, energy_average_geometric


def read_block(df, binno, block):
	"""Reads block."""
	fluka_value = df[binno*(block-1):binno*block]['diffF']
	fluka_err = df[binno*(block-1):binno*block]['relerr']

	return fluka_value, fluka_err

REGION = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10', 'R11', 'R12', 'R13',
		  'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6',
		  'G0', 'G0old',
		  'PC0', 'PCint',
		  'T0u', 'T0', 'T0d',
		  # 'F0u', 'F0', 'F0d'
          ]

def load_particle(particle, directory, input_filename):
	df = read_data_frame(directory, input_filename, particle)

	binno = 80  # Number of Energy bins

	### Energy Array
	Emin, Emax, DE, Eave, Egeo = read_energy(df, binno)

	### Fluence Array
	data_locations = list()
	err_locations = list()
	for i in range(len(REGION)):
		fluka, err = read_block(df, binno, i + 1)

		data_locations.append(fluka)
		err_locations.append(err)

	print(f'\nFinished loading {particle} spectra file! \n')
	return DE, Emin, Egeo, data_locations, err_locations


def load_neutron(directory, config):
	df = read_data_frame(directory, config, 'neutron')

	binno = 304  # Number of Energy bins

	# Energy Array
	Emin, Emax, DE, Eave, Egeo = read_energy(df, binno)

	# Fluence Array
	data_locations = list()
	err_locations = list()
	number_of_R_racks = 13
	# for i in range(number_of_R_racks):
	for i in range(len(REGION)):
		fluka, err = read_block(df, binno, i + 1)

		data_locations.append(fluka)
		err_locations.append(err)

	# G0
	# Note: The binning is changing from this point
	# block = 1
	# lastblock = 304*number_of_R_racks
	# binno = 279
	#
	# Emin_2, Emax_2, DE_2, Eave_2, Egeo_2 = read_energy_different_binning(df, binno, block, lastblock)
	#
	# for i in range(len(REGION)-number_of_R_racks):
	# 	block = i + 1
	# 	fluka, err = read_block_different_binning(df, binno, block, lastblock)
	#
	# 	data_locations.append(fluka)
	# 	err_locations.append(err)

	print('\nFinished loading NEUTRON spectra file! \n')
	# return DE, Emin, Egeo, data_locations, err_locations, DE_2, Emin_2, Egeo_2
	return DE, Emin, Egeo, data_locations, err_locations


def plot_lethargy_spectra(directory, input_filename, rack, saveplot, fileEXT):
	# Load RAW data
	DE_proton, Emin_proton, Egeo_proton, proton_RAW, proton_err = load_particle("proton", directory, input_filename)
	# DE_neutron, Emin_neutron, Egeo_neutron, neutron_RAW, neutron_err, DE_neutron2, Emin_neutron2, Egeo_neutron2 = load_neutron(directory, config)
	DE_neutron, Emin_neutron, Egeo_neutron, neutron_RAW, neutron_err = load_neutron(directory, input_filename)
	DE_photon, Emin_photon, Egeo_photon, photon_RAW, photon_err = load_particle("photon", directory, input_filename)
	DE_electron, Emin_electron, Egeo_electron, electron_RAW, electron_err = load_particle("electron", directory, input_filename)
	DE_positron, Emin_positron, Egeo_positron, positron_RAW, positron_err = load_particle("positron", directory, input_filename)
	DE_pion_positive, Emin_pion_positive, Egeo_pion_positive, pion_positive_RAW, pion_positive_err = load_particle("pion_positive", directory, input_filename)
	DE_pion_negative, Emin_pion_negative, Egeo_pion_negative, pion_negative_RAW, pion_negative_err = load_particle("pion_negative", directory, input_filename)
	DE_muon_positive, Emin_muon_positive, Egeo_muon_positive, muon_positive_RAW, muon_positive_err = load_particle("muon_positive", directory, input_filename)
	DE_muon_negative, Emin_muon_negative, Egeo_muon_negative, muon_negative_RAW, muon_negative_err = load_particle("muon_negative", directory, input_filename)
	DE_kaon_positive, Emin_kaon_positive, Egeo_kaon_positive, kaon_positive_RAW, kaon_positive_err = load_particle("kaon_positive", directory, input_filename)
	DE_kaon_negative, Emin_kaon_negative, Egeo_kaon_negative, kaon_negative_RAW, kaon_negative_err = load_particle("kaon_negative", directory, input_filename)
	DE_kaon_neutral, Emin_kaon_neutral, Egeo_kaon_neutral, kaon_neutral_RAW, kaon_neutral_err = load_particle("kaon_neutral", directory, input_filename)
	DE_alpha, Emin_alpha, Egeo_alpha, alpha_RAW, alpha_err = load_particle("alpha", directory, input_filename)
	DE_heh, Emin_heh, Egeo_heh, heh_RAW, heh_err = load_particle("heheq", directory, input_filename)


	# Computation of the data for the selected test location
	# loc, Egeo_neutron, DE_neutron = extract_value_at_rack(rack, Egeo_neutron, DE_neutron, Egeo_neutron2, DE_neutron2)
	loc = REGION.index(rack)

	# Plot of the selected test location
	print('Plot for test location:', rack)

	# Create the Lethargy arrays
	pr, pr_err, Egeo_proton = compute_lethargy_values(proton_RAW[loc], proton_err[loc], Egeo_proton, energy_thresholds["proton"])
	ne, ne_err, Egeo_neutron = compute_lethargy_values(neutron_RAW[loc], neutron_err[loc], Egeo_neutron, energy_thresholds["neutron"])
	ph, ph_err, Egeo_photon = compute_lethargy_values(photon_RAW[loc], photon_err[loc], Egeo_photon, energy_thresholds["photon"])
	el, el_err, Egeo_electron = compute_lethargy_values(electron_RAW[loc], electron_err[loc], Egeo_electron, energy_thresholds["electron"])
	positron, positron_err, Egeo_positron = compute_lethargy_values(positron_RAW[loc], positron_err[loc], Egeo_positron, energy_thresholds["positron"])
	pion_positive, pion_positive_err, Egeo_pion_positive = compute_lethargy_values(pion_positive_RAW[loc], pion_positive_err[loc], Egeo_pion_positive, energy_thresholds["pion_positive"])
	pion_negative, pion_negative_err, Egeo_pion_negative = compute_lethargy_values(pion_negative_RAW[loc], pion_negative_err[loc], Egeo_pion_negative, energy_thresholds["pion_negative"])
	muon_positive, muon_positive_err, Egeo_muon_positive = compute_lethargy_values(muon_positive_RAW[loc], muon_positive_err[loc], Egeo_muon_positive, energy_thresholds["muon_positive"])
	muon_negative, muon_negative_err, Egeo_muon_negative = compute_lethargy_values(muon_negative_RAW[loc], muon_negative_err[loc], Egeo_muon_negative, energy_thresholds["muon_negative"])
	kaon_positive, kaon_positive_err, Egeo_kaon_positive = compute_lethargy_values(kaon_positive_RAW[loc], kaon_positive_err[loc], Egeo_kaon_positive, energy_thresholds["kaon_positive"])
	kaon_negative, kaon_negative_err, Egeo_kaon_negative = compute_lethargy_values(kaon_negative_RAW[loc], kaon_negative_err[loc], Egeo_kaon_negative, energy_thresholds["kaon_negative"])
	kaon_neutral, kaon_neutral_err, Egeo_kaon_neutral = compute_lethargy_values(kaon_neutral_RAW[loc], kaon_neutral_err[loc], Egeo_kaon_neutral, energy_thresholds["kaon_neutral"])
	alpha, alpha_err, Egeo_alpha = compute_lethargy_values(alpha_RAW[loc], alpha_err[loc], Egeo_alpha, energy_thresholds["alpha"])
	heh, he_err, Egeo_heh = compute_lethargy_values(heh_RAW[loc], heh_err[loc], Egeo_heh, energy_thresholds["heheq"])

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

	fig1 = plt.figure(num=1, figsize=(width, height))
	sub1 = fig1.add_subplot(111, aspect='auto')
	# PAPER STYLE
	sub1.errorbar(Egeo_neutron, ne, yerr=(ne_err, ne_err), linestyle='', linewidth=llwidth, color=colors[0],
	                    marker=markers[0], markeredgecolor='black', markersize=msize1, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='n')
	sub1.errorbar(Egeo_proton, pr, yerr=(pr_err, pr_err), linestyle='', linewidth=llwidth, color=colors[1],
	                    marker=markers[1], markeredgecolor='black', markersize=msize2, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='p')
	sub1.errorbar(Egeo_photon, ph, yerr=(ph_err, ph_err), linestyle='', linewidth=llwidth, color=colors[2],
	                    marker=markers[2], markeredgecolor='black', markersize=msize1, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='$\mathregular{\gamma}$')
	sub1.errorbar(Egeo_electron, el, yerr=(el_err, el_err), linestyle='', linewidth=llwidth, color=colors[3],
	                    marker=markers[3], markeredgecolor='black', markersize=msize1, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='$\mathregular{e^{-}}$')
	sub1.errorbar(Egeo_positron, positron, yerr=(positron_err, positron_err), linestyle='', linewidth=llwidth, color=colors[4],
	                    marker=markers[4], markeredgecolor='black', markersize=msize1, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='$\mathregular{e^{+}}$')
	sub1.errorbar(Egeo_pion_negative, pion_negative, yerr=(pion_negative_err, pion_negative_err), linestyle='', linewidth=llwidth, color=colors[5],
	                    marker=markers[5], markeredgecolor='black', markersize=msize2, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='$\mathregular{\pi^{-}}$')
	sub1.errorbar(Egeo_pion_positive, pion_positive, yerr=(pion_positive_err, pion_positive_err), linestyle='', linewidth=llwidth, color=colors[6],
	                    marker=markers[6], markeredgecolor='black', markersize=msize2, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='$\mathregular{\pi^{+}}$')
	# sub1.errorbar(Egeo_muon_negative, muon_negative, yerr=(muon_negative_err, muon_negative_err), linestyle='', linewidth=llwidth, color=colors[7],
	#                     marker=markers[7], markeredgecolor='black', markersize=msize1, markeredgewidth=mewidth, drawstyle='steps-mid',
	#                     label='$\mathregular{\mu^{-}}$')
	# sub1.errorbar(Egeo_muon_positive, muon_positive, yerr=(muon_positive_err, muon_positive_err), linestyle='', linewidth=llwidth, color=colors[8],
	#                     marker=markers[8], markeredgecolor='black', markersize=msize1, markeredgewidth=mewidth, drawstyle='steps-mid',
	#                     label='$\mathregular{\mu^{+}}$')
	sub1.errorbar(Egeo_kaon_negative, kaon_negative, yerr=(kaon_negative_err, kaon_negative_err), linestyle='', linewidth=llwidth, color=colors[9],
	                    marker=markers[9], markeredgecolor='black', markersize=msize1, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='$\mathregular{k^{-}}$')
	sub1.errorbar(Egeo_kaon_positive, kaon_positive, yerr=(kaon_positive_err, kaon_positive_err), linestyle='', linewidth=llwidth, color=colors[10],
	                    marker=markers[10], markeredgecolor='black', markersize=msize1, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='$\mathregular{k^{+}}$')
	sub1.errorbar(Egeo_kaon_neutral, kaon_neutral, yerr=(kaon_neutral_err, kaon_neutral_err), linestyle='', linewidth=llwidth, color=colors[11],
	                    marker=markers[11], markeredgecolor='black', markersize=msize1, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='$\mathregular{k^{0}}$')
	sub1.errorbar(Egeo_alpha, alpha, yerr=(alpha_err, alpha_err), linestyle='', linewidth=llwidth, color=colors[12],
	                    marker=markers[12], markeredgecolor='black', markersize=msize1, markeredgewidth=mewidth,
	                    drawstyle='steps-mid', label='$\mathregular{\\alpha}$')

	# Create the labels
	# sub1.axis([1e-14,1e2, 1e-2, 1e-13], fontsize=fsize)
	axis_font = {'fontname': 'Times New Roman', 'size': 18, 'weight': 'normal'}
	# sub1.axis([1e-13, 1e2, 1e-12, 1e-2])  # fontsize=fsize
	sub1.loglog()
	sub1.tick_params(labelsize=fsize)
	plt.xlabel('Energy [GeV]', fontsize=fsize)
	plt.ylabel('Fluence per unit lethargy $\mathregular{[cm^{-2} \, / POT]}$', fontsize=fsize)
	# plt.title('CHARM - Configuration {0:<8s} - Test Location {1:<5s}'.format(config, rack), fontstyle='italic',fontsize=fsize)

	sub1.set_xlim(xmax=1e+5)

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
	plt.legend(fontsize=fsize, loc=7)

	# sub1.annotate('(b)', xy=(0.90,0.90), xycoords='axes fraction', fontweight='bold', fontsize=9)

	# Make/Save the plot
	plt.tight_layout(pad=3, w_pad=1.0, h_pad=2)

	if saveplot:
		outfile = f"plots/locations/CHARM_{rack}_TNS.{fileEXT}"
		plt.savefig(outfile, dpi=100, bbox_inches='tight')
		# plt.savefig(outfile)
		plt.close()
	else:
		plt.show()


def main():
	# Absolute path to the data files
	# directory = '/scratch/hlillepa/2022-summer-student-project/charm_test/output_cu/analysis'
	# input_filename = "charm_v2022_cu"

	directory = directory = '/scratch/hlillepa/2022-summer-student-project/charm_test/output15_usrtrack/analysis'
	input_filename = 'charm_new15'

	# Activate Plotting of the spectra for a given position
	saveplot = True
	fileEXT = 'png'  # pdf or jpeg or png

	# Set the test location for plotting.
	for rack in REGION:
		# plot_lethargy_spectra_legacy(directory, target, config, rack, saveplot, fileEXT)
		plot_lethargy_spectra(directory, input_filename, rack, saveplot, fileEXT)


if __name__ == "__main__":
	main()
