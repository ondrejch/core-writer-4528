#! /usr/bin/python
#
# Lattice: Script containing a function to write a hexagonal SERPENT lattice for the MSR project.

import math
from math import ceil,sqrt,pi
def write_lattice(rmax = 100, pitch = 11.5, core = 200, nlat = '1', \
    ub = '1', uf = '2', uc = '3', \
    nf = 6, nb = 1):
	'''Function to write the hexagonal lattice of the given type
	Accepts as input:
		nlat:		number identifying the lattice
		ub, uf, uc:	the universe numbers of the cells for the MSR core:
					blanket, fuel, central (int);
		nf, nb:		the number of fuel/blanket cells in radius (integer);
		p:			the pitch (cm) between cells (float)
	Outputs:
		lattice:	string containing the generated lattice card'''
	
	# Old lattice calculation
	'''	# The user will normally specify a lone integer.
	# Correct the universes to strings plus spaces
	ub = str(ub) + ' '
	uf = str(uf) + ' '
	uc = str(uc) + ' '

	# Total width: [nb]-[nf]-center-[nf]-[nb]
	nc = 2*(nf + nb) + 1
	# Adjust nx to account for blanket buffers
	nx = int(1.5*nc)
	ny = nc + 2 

	# Scaling factor: relates circumradius to hexagon radius
	# This is necessary because HEXYC are taller than they are wide.
	# Thus, we need fewer of them in on dimension.
	c = math.sqrt(3.0)/2.0 	# about 0.866


	# Structure of a SERPENT lattice card
	\'''lat <u0> <type> <x0> <y0> <nx> <ny> <p>
	#
	# <u0> is the universe number of the lattice
	# <type> is the lattice type: 2 is our hexagon (hexyc)
	# <x0>, <y0> are the x/y coords of the origin--ours are (0,0)
	# <nx>, <ny> are the number of lattice elements in each direction
	# <p> is the lattice pitch\'''

	# The card begins
	lattice = \'''
% This lattice was generated by a script.
lat {nlat}  2  0.0 0.0  {nx} {ny}  {p} \n\'''

	# Make a buffer row of blanket salt
	lattice += nx*ub + '\n'

	# Here is the interesting part.
	# The fuel cells need to follow the equation of a circle
	for j in range(nc):
		# Move down the lattice core in the y
		
		# Create a blank row to fill
		#row = ''			# non-indented version
		row = '' + (j+1)*' ' 	# indent it a bit to make the hex

		# Create a buffer of whitespace and blanket cells
		buf = int(round(nc/2.0))
		bl = int(math.floor((nc - j)/2))*ub 		# left buffer
		br = int(math.floor(j/2.0))*ub + (nc - j)*' '	# right buffer
		
		# Ones that didn't work
		#bl = int(round((nc - j)/2.0))*ub			# left buffer
		#br = int(math.floor((nc - j)/2))*' '	# right buffer

		
		#if False: #j == nf + nb:
		if j == nf + nb:
			# Then this is the special central row
			row += bl + (nb+1)*ub + (nf-1)*uf + uc + (nf-1)*uf + (nb+1)*ub + br
		else:
			# Cylindrical coordinates
			# Trace out the boundary of the circle
			r = nc/2.0 - 0.5
			y = abs(r-float(j)) #- 0.5
			x = math.sqrt(r**2 - y**2)*c

			# x2 is the integer number of fuel cells in row
			if x <= 0:
				x2 = 0
			else:
				#x2 = int(round(2.0*x)) - nb
				x2 = int(math.ceil(2.0*x)) - nb
			
			# The width of the blanket salt on each side	
			w1 = int(math.floor( (nc - x2)/2.0) + (nc-j)/2.0)
			w2 = int(math.ceil( (nc - x2)/2.0) + (j)/2.0)
			#w1 = 0
			#w2 = 0

			# Put it all together to make the new row
			row += (w1)*ub + x2*uf + (w2)*ub
		

		# Then, finally, concatenate the new row to the lattice
		lattice += row + '\n'


	# Make a final buffer row of blanket salt
	lattice += (j+2)*' ' + nx*ub + '\n'
	
	# Plug in the numbers and submit.
	lattice = lattice.format(**locals())
	return(lattice)'''

	# New lattice calculation !!!WORK IN PROGRESS!!!

	''' lattice dependencaies:
	pitch
	rmax
	core
	'''	
	DEBUG=False
#	core = 2
#	pitch = 10
#	rmax = 100
# adjust universe numbers to include a space
	ub = str(ub) + ' '
	uf = str(uf) + ' '
	uc = str(uc) + ' '
	# and these go in a x type hex lattice
	n = 2*int(ceil(core/pitch))+10  #seems to do the trick
	n += 1 if n%2==0 else 0 #make sure it is an odd pin num to center it
	
	# needs to be changed for appropriate lattice structures
	lattice = '''
% This lattice was generated by a script.
lat {nlat} 2  0.0 0.0  {n} {n}  {pitch} \n'''

# iterate through lattice positions, only place pins where they lie fully in-mod
        # lattice starts at bottom left
	x0 = 0.0
	y0 = 0.0

        # yes, yes it did take well over 7 hours to get x0, y0 right
	x0 -= (n+0.5*n)*pitch/2.0 - .75*pitch
	y0 -= n*sqrt(3.0)/4.0*pitch - sqrt(3.0)/4.0*pitch
	for i in range(n): #y positions
		thisline='' # lattice line in input file
		thisline += ' '*i #indent to show hex better
		for j in range(n): # x positions
			y = y0 + i*sqrt(3.0)/2.0*pitch #y pos
			x = x0 + (j+0.5*i)*pitch #x pos
			lat_r = sqrt(x**2+y**2)  #radius from lattice center
			if DEBUG:
				print x, y, lat_r, lat_r < rmax
			if lat_r < rmax:             
                   # if lattice radius is in the center, write the unique central channel
				if int(x) == 0 and int(y) == 0:
					#thisline += '0 '
					thisline += uc
				else:
					#thisline += '5 ' #add a channel to geometry
					thisline += uf
                  # self.channels.append(CircularChannel(x,y,r)) #add channel to core object
			else:
				#thisline += '6 ' #add a blank
				thisline += ub
		thisline += '\n'
		lattice += thisline
		
	lattice = lattice.format(**locals())
	return(lattice)

# This executes if someone tries to run the module
if __name__ == '__main__':
	print "This is a module which generates lattice cards for Serpent."
	a = raw_input("Press 'enter' to exit, or type anything to run a test. ")
	if a:
		print '_____________________________\n'
		print write_lattice()
