#!/usr/bin/env python3
#
# Deck.py
# A script that generates the Serpent input deck for our SMTF-MSBR
	
import math #HOW DARE SOMEONE USE MATH AND NOT NUMPY. 
import lattice, surfs, cells, materials


def write_deck(fsf = 0.07, relba = 0.08,\
    pitch = 11.500, \
    slit = 0.2, temp=700, r2 = 3.3, rs = 0.9, \
    rfuel = 150, rcore = 215, zcore = 400, refl_ht = 100, \
    name = 'Test deck'):
	'''Write the actual Serpent deck
	Inputs:
* channel_pitch:  hexagonal pitch of fuel cells [cm]
* channel_r2: inner radius of the outer fuel channel [cm]
* channel_top_r: radius of the channel holder knob [cm]                          
* channel_top_h: height the channel holder knob [cm]                          
* blanket_slit:   thickness of blanket salt slit [cm]
* salt_fraction:  salt fraction in the central channel [fraction]
* crod_positions: tuple of [x,y] describing CR channel position(s), indexed in channel numbers. [0,0] is the center.
* r_core:  radius of the inner core (fuel only) [cm]
* r_tank:  radius of the core + hastelloy container [cm]
* h_plena: height of each of the lower plena [cm]
* z_core: height of the graphite channels [cm]
* z_tank: height of the hastelloy tank [cm]
* tank_thick: hastelloy thickness [cm]
* holder_thick: holder plate thickness [cm]
* name:		string containing name of this model

# OLD		fsf:		fuel salt fraction
#		pitch:		hexagonal pitch (cm)
#		slit:		blanket salt slit width (cm)
#		rfuel:		radius of the inner core of fuel cells (cm)	
#		rcore:		radius of the entire core, fuel + blanket (cm)
#		zcore:		height of the active core, fuel (cm)
#		refl_ht:	height of the axial reflector (cm)
#		name:		string containing name of this model
	Outputs:
		output:		String containing the entire input deck
	'''	
# TODO 
#    channel_r1 = 
#    channel_r3 = 

	# Read the initial values from some external source.
	# Right now, I'm just plugging them in to test the script.
	FSF = 	fsf
	PITCH = pitch # cm
	SLIT = 	slit  # cm
	RELBA = relba
	TEMP = temp

	fuel_cells = int(rfuel/PITCH)
	blan_cells = 1

	# radius (outer): auxiliary fuel channel radius
    #ro = ri / math.sqrt(6)
	ro = 1.1
	#--------------------------------
	# Begin writing the input deck
	
	output = '''\
set title "{name}"
/*This is a  model of the:
Molten Salt iso-Breeder Reactor
Core Design Team: Dallas Moser, Igor Gussev
Reprocessing: Devon Drey
Advisor: Dr. Ondrej Chvala
*/\n
	'''
	
		
	LATS = range(33,33+7)
	
	# define the relevant universe numbers
	ub, uf, uc, uup, uuc, ulc = range(1,7)
	ul1, ul2, ul3, ul4 = range(25,29)
	ulp = 10
	uh  = 11
	ubb = 7 # blank blanket universe
	uhsu = 8 # holding shafts upper
	uhsl = 9 # holding shafts lower
	uhp  = 12 # holding plate
	
	
	# Tuple of all the universe numbers
	UNIVERSES = (
	ub, 	# blanket cell
	uf,		# fuel cell 
	uc, 	# control rod
	uup,	# upper channel fuel
	ul1,	# lower channel 1
	ul2, 	# lower channel 2
	ul3, 	# inlet plenum penetration
	ul4, 	# outlet plenum penetration
	ulp, 	# lower plenum (pure fuel)
	uuc,	# upper control
	ulc,	# lower	control
	uh)		# pure hastelloy hex
	
	rcore_inner = rcore - 2*2.54
	rcore_outer = rcore
	rfuel = 150
	plenum_vol = 37*28316.8 	# 37 ft^3 to cm^2
	# Height of each plenum: inlet and outlet
	plenum_ht = plenum_vol / (2*math.pi*rfuel**2)
	gt = 6*2.54 # thickness of graphite: cm
	ht = 3*2.54 # thickness of hastelloy, placeholder
	fuel_cells = int(rfuel/PITCH)
	blan_cells = 1
	rgref = rcore + gt
	rhast = rgref + ht
	
	surface_cards = surfs.write_surfs(FSF, RELBA, PITCH, SLIT, TEMP, ro, r2, rs, \
									  rfuel, rcore_inner, rcore_outer, \
									  zcore, plenum_ht, refl_ht)
	output += surface_cards
	
	cell_cards = cells.write_cells(UNIVERSES, LATS, 30, 31, 32, 33)
	output += cell_cards
	
	# Create the middle/active core
	lattice_cards = lattice.write_lattice(rfuel, PITCH, rcore_inner, LATS[0], ub, uf, uc, fuel_cells, blan_cells)
	# Create the upper plenum
	lattice_cards += lattice.write_lattice(rfuel, PITCH, rcore_inner,LATS[1], ub, uup, uuc, fuel_cells, blan_cells)
	# Create the lower level -1
	lattice_cards += lattice.write_lattice(rfuel, PITCH, rcore_inner,LATS[2], ubb, ul1, ulc, fuel_cells, blan_cells)
	# Create the lower level -2
	lattice_cards += lattice.write_lattice(rfuel, PITCH, rcore_inner,LATS[3], ubb, ul2, ulc, fuel_cells, blan_cells)
	# Create the lower level -3: penetration to inlet plenum
	lattice_cards += lattice.write_lattice(rfuel, PITCH, rcore_inner,LATS[4], ubb, ul3, uh,  fuel_cells, blan_cells)
	# Create the lower level -4: penetration to the outlet plenum
	lattice_cards += lattice.write_lattice(rfuel, PITCH, rcore_inner,LATS[5], uh, ul4, uh,  fuel_cells, blan_cells)
	# Create the lower fuel plena (identical for both)
	lattice_cards += lattice.write_lattice(rfuel, PITCH, rcore_inner,LATS[6], uh, ulp, ulp, fuel_cells, blan_cells)
	# Create the upper holding shafts
	lattice_cards += lattice.write_lattice(rfuel, PITCH, rcore_inner,40, ubb, uhsu, uhsu, fuel_cells, blan_cells)
	# Create lower holding shafts
	lattice_cards += lattice.write_lattice(rfuel, PITCH, rcore_inner,41, ubb, uhsl, uhsl, fuel_cells, blan_cells)
	# Create holding plate
	lattice_cards += lattice.write_lattice(rfuel, PITCH, rcore_inner,42, 11, uhp, uhp, fuel_cells, blan_cells)
	output += lattice_cards
	
	mat_cards = materials.write_materials(TEMP)
	output += mat_cards
	
	
	data_cards = '''
set gcu  0
%set sym  12
set nfg  2  0.625E-6

% --- Neutron population and criticality cycles:

set pop 2000 500 20


%% --- Data Libraries
set acelib "sss_endfb7u.xsdata"
set nfylib "sss_endfb7.nfy"
set declib "sss_endfb7.dec"
	'''
	output += data_cards
	
	
	plot_cards = '''
% PLOTS
plot 1 3000 3000 0  -300 300  -80 560
plot 2 3000 3000 0  -300 300  -80 560
plot 3 3000 3000 29 %[250 -100 100 -100 100]
	'''
	output += plot_cards
	
	output = output.format(**locals())

	return output


if __name__ == '__main__':
	print("This is the Serpent deck writing function for the MSR project.")
	raw_input("Press Ctrl+C to exit, or Enter to test it. ")
	print(write_deck())
