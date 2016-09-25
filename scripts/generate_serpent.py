#! /usr/bin/python
#
# Generate Serpent
# A script that generates the Serpent input deck for our SMTF-MSBR
#
# Calls the deck_writing function for each case we want to run

import deck
import os

# Parameters from the MCNP optimization
FSF = 0.07
PITCH = 11.500
R2 = 3.3
#SLIT = 0.323
RFUEL = 150
RCORE = 215
ZCORE = 400
ZREFL = 100

FILENAME = "msbr.inp"
QSUBLOC = "/home/tlabossi/MSR/qsub.sh"

# Before iterating, do this quick test
fails = os.system('mkdir s0.108')
if fails:
	print "ERROR: Cannot create file/directory for writing. Aborting."
	exit()
# Otherwise, clean up and continue with the script
os.system('rmdir s0.108')


# From .108 cm to .327 cm
slits = range(108, 327, 3)
#slits = range(108, 111, 3) #just to test
# Generate for each slit width we want to try
for s in slits:
	# Scale by a factor of 1000 and convert to float
	s *= .001

	title = "SMTFMSBR: slit = " + str(s)


	# Make the deck
	serp_deck = deck.write_deck(fsf = FSF, pitch = PITCH, slit = s, \
				rfuel = RFUEL, rcore = RCORE, r2 = R2, zcore = ZCORE, refl_ht = ZREFL,
				name = title)
	

	dirname = 's' + str(s) + '/'
	fails = os.system('mkdir ' + dirname)
	if fails:
		print "Could not create directory", dirname, " - trying the next one"
	else:
		try:
			fname = dirname + FILENAME
			f = open(fname, 'w')
			f.write(serp_deck)
			f.close()
		except IOError as e:
			print "Unable to write to file", fname
			print e
		else:
			# Do the qsub bit
			os.system('cd ' + dirname + '; cp ' + QSUBLOC + ' .; qsub qsub.sh')

print "All done."
