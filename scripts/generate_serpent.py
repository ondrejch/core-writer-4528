#! /usr/bin/python
#
# Generate Serpent
# A script that generates the Serpent input deck for our SMTF-MSBR
#
# Calls the deck_writing function for each case we want to run

import deck
import os
import numpy as np

# Parameters from the MCNP optimization
FSF = 0.07
PITCH = 11.500
R2 = 3.3
#SLIT = 0.323
RFUEL = 150
RCORE = 215
ZCORE = 400
ZREFL = 100


#Job submission settings:
FILENAME = "msbr.inp"
qsubtemplate = 'qsubtemplate.txt'
nperjob = 2 #num nodes per job
ncpu = 8 #cpu per job
queue = 'gen5'

#now the qsub text is constructed.
qtext=[]
with open(qsubtemplate, 'r') as temp:
    text=temp.read()
    for line in text:
       qtext.append(line.format(**locals())) 

# From .108 cm to .327 cm
# 73 different slit widths are attempted.
slits = np.linspace(0.108, 0.327, 73)

for s in slits:

    title = "SMTFMSBR: slit = " + str(s)


    # Make the deck
    serp_deck = deck.write_deck(fsf = FSF, pitch = PITCH, slit = s,
                            rfuel = RFUEL, rcore = RCORE, r2 = R2, zcore = ZCORE, refl_ht = ZREFL,
                            name = title)
    
    # create the directory
    dirname = 's' + str(s) + '/'
    os.mkdir(dirname)

    # go into the directory to run this lattice in
    os.chdir(dirname)

    # now write out the input file.
    with open(FILENAME, 'w') as f:
        f.write(serp_deck)

    # place a qsub script with a descriptive name here.
    with open('{}.sh'.format(dirname),'w') as f:
        f.write(qtext)

    #and finally, submit the job.
    #os.system('qsub {}'.format(FILENAME+'.sh'))

    #and now return to the original directory.
    os.chdir('..')
        
                        
                        

print "All done."
