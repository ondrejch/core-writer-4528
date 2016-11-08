#! /usr/bin/python
#
# Surfs
# Calculate dimensions for each of the surfaces based off the MSBR cell dimensions
# and output the surface cards for the Serpent input deck.


import math


def write_surfs(pitch, slit, d, ro, rox, roy, r1, r2, c, rfuel, rcore, rgref, rhast, zcore, pht, zrefl):
	'''Function to write the surfaces for our MSBR Serpent model
	Inputs:
		pitch:  hexagonal pitch of fuel cells
		slit:   thickness of blanket salt slit
		d:      circumradius of hexagon
		ri:     central fuel channel radius
		r1:     central fuel channel radius
		r2:     radius of concentric graphite ring
		r3:     radius of outer fuel channel ring
		ro:     auxiliary fuel channel radius
		rox:    x coordinate of aux center
		roy:    y coordinate of aux center
		c:      ratio of (rox,roy) to d
		rfuel:  radius of the inner core (fuel only)  (rfuel)
		rcore:  radius of the outer core (fuel+blanket)
		rgref:  radius of the core + graphite reflector (rcore_inner)
		rhast:  radius of the core + graphite + hastelloy container (rcore_outer)
		pht:	height of each of the lower plena
		zcore:	height of the core
		zrefl:	height of the axial reflector
	Output:
		surfaces:   string containing the surface cards for the MSR'''
	
'''	# Radial reflector scaling term
	rs = 0.9
		
	# r2 is the outer fuel radius; thast = hastelloy thickness (1/8 in)
	thast = 1.0/8 * 2.54							# hastelloy thickness
	r2 = math.sqrt(2*r1**2 + 2*r1*thast + thast**2) # outer fuel cylinder
	# Radius of outer fuel ring with equal volume to inner fuel channel
	r3 = math.sqrt(r1**2 + r2**2)
	rdiff = (r3 - r2)
	# Establish a few additional dimensions
	hexs = pitch/2.0    # radius of cell, outside slit
	#ry = c*d
	hexg = hexs - slit   # radius of graphite, inside slit
	ry = c*d			# y coord of vertical channel

	# hexf: top channel hexagon
	hexf = math.sqrt(3)/2.0*(d*c) + rdiff
	zhexf1 = zcore + 2.0*rdiff	# top of fuel channel (salt)
	# Make the "roof" as thick as the "wall" of the graphite hexagon
	zhexf2 = zhexf1 + (hexg - hexf) # top of fuel channel (graphite)
	# Top of the core: Axial reflector
	zrefltop = zhexf2 + zrefl
	zgreftop = zrefltop + (rgref - rcore)
	zhasttop = zgreftop + (rhast - rgref)


	# Bottom (floor) with the transition later
	# We want this to converge to another channel with concentric
	# cylinders at the lower plenum
	# Graphite thickness between central channel and aux channel
	gt1 = d*c - rdiff - r1
	# Next layer down concentric circles with half the graphite thickness
	# thast = hastelloy thickness (1/8 in)
	thast = 1.0/8 * 2.54							# hastelloy thickness
	rh = r1 + thast									# radius of hastelloy cyl
	rg = (gt1/2.0 + r2)    				   			# radius of graphite hex
	# And cut off the bottom 10 inches below that
	ztrans1 = 0 - 10*2.54
	# Next layer down, the concentric circles w/ hastelloy
	rh2 = r2 + thast		# Outer hastelloy pipe
	# Cut off this layer 3 inches down
	ztrans2 = ztrans1 - 3*2.54

	# Then, the lower plena
	zitop = ztrans2 - pht	# z of the top of the inlet plenum
	zibot = zitop - pht		# z of the bottom of the inlet plenum
	zotop = zibot - pht		# z of top of outlet plenum
	zobot = zotop - pht		# z of bottom of outlet plenum
	# Then, the very bottom of the entire core
	zbot = zobot - (rhast - rgref)

	# Reflector apothem 
	rr = rs*hexg
	
	# Top reflector thickness
	botrt = zhexf2
	toprt = botrt + 10
'''

'''	surfaces = \'''
%------define the hexagon and fuel channel cells----
surf 10 hexxc 0   0   {hexg}	     % HEX FOR GRAPHITE
surf 11 hexxc 0   0   {hexs}	      % HEX FOR SLIT
surf 12 hexxc 0   0   {hexf}          % HEX FOR FUEL used in channel cap
surf 20 cyl   0   0   {ri}	    % CENTER HOLE
surf 21 cyl   0   0   {r2}     % INTERMEDIATE GRAPHITE RING
surf 22 cyl   0   0   {r3}     % OUTER FUEL RING
surf 27 pz    0			    % BOTTOM OF CORE
surf 28 pz    {zcore}		     % HEIGHT OF CORE
surf 30 cyl   0   0   {rfuel}   % CYLINDRICAL BOUNDS FOR FUEL LATTICE
surf 31 cyl   0   0   {rcore}   % CYLINDRICAL BOUNDS FOR FUEL AND BLANKET
surf 32 cyl   0   0   {rgref}   % CYLINDRICAL BOUNDS FOR GRAPHITE + PREVIOUS 2
surf 33 cyl   0   0   {rhast}   % ENTIRE CORE CYLIDNRICAL BOUNDS
surf 41 pz    {zhexf1}          % BOTTOM OF GRAPHITE CAP
surf 42 pz    {zhexf2}          % TOP OF GRAPHITE CAP
surf 51 cone  0 0 0 {hexg}  -{hexg}  % NOT USED
surf 52 pz    {ztrans1}       % BOTTOM OF A LATTICE
surf 53 cyl   0   0   {rh}   % hastelloy tube radius
surf 54 cyl   0   0   {r2}   % outer fuel cyl radius
surf 55 hexxc 0   0   {rg}	% Hex for fuel transition
surf 56 pz    {ztrans2}
surf 57 cyl   0   0   {rh2}
surf 60 pz    {zbot}
surf 61 pz    {zitop}
surf 62 pz    {zibot}
surf 63 pz    {zotop}
surf 64 pz    {zobot}
surf 80 pz    {zrefltop}
surf 81 pz    {zgreftop}
surf 82 pz    {zhasttop}
surf 90 cyl   0   {ry}   {ro}	    % OUTER HOLES x 6 CONTROL RODS
surf 91 cyl   0  -{ry}   {ro}
surf 92 cyl   {rox}  {roy}  {ro}
surf 93 cyl  -{rox}  {roy}  {ro}
surf 94 cyl  -{rox} -{roy}  {ro}
surf 95 cyl   {rox} -{roy}  {ro}
surf 101 hexxc 0   0   {rr}	     % HEX FOR RADIAL GRAPHITE REFLECTOR
surf 102 cylz 0  0 {rfuel} {botrt} {toprt}
\'''
'''

	surfaces = '''
%------ main universe ------
surf 1 cyl   0   0   {rfuel}         % CYLINDRICAL BOUNDS FOR FUEL LATTICE
surf 2 cyl   0   0   {rcore_inner}   % CYLINDRICAL BOUNDS FOR FUEL AND REFLECTOR
surf 3 cyl   0   0   {rcore_outer}   % CYLINDRICAL BOUNDS FOR ENTIRE CORE

%------ graphite hexagon and fuel channel cells ------	
surf 101 hexxc 0   0   {hexg}	     % HEX FOR GRAPHITE
surf 102 hexxc 0   0   {hexs}	     % HEX FOR SLIT
surf 103 cyl   0   0   {r1}	     % CENTER HOLE
surf 104 cyl   0   0   {r2}          % INTERMEDIATE GRAPHITE RING
surf 105 cyl   0   0   {r3}          % OUTER FUEL RING
surf 198 pz    0		     % BOTTOM OF CORE
surf 199 pz    {zcore}		     % HEIGHT OF CORE

%------ graphite hexagon and control rod channels ------
surf 201 hexxc 0   0   {hexg}	     % HEX FOR GRAPHITE
surf 202 hexxc 0   0   {hexs}	     % HEX FOR SLIT
surf 203 cyl   0   0   {r1}  	     % CENTER HOLE
surf 204 cyl   0   {ry}   {ro}	     % OUTER HOLES x 6 CONTROL RODS
surf 205 cyl   0  -{ry}   {ro}       %           ||
surf 206 cyl   {rox}  {roy}  {ro}    %           ||
surf 207 cyl  -{rox}  {roy}  {ro}    %          _||_
surf 208 cyl  -{rox} -{roy}  {ro}    %          \  /
surf 209 cyl   {rox} -{roy}  {ro}    % __________\/____________
surf 298 pz    0		     % BOTTOM OF CORE
surf 299 pz    {zcore}		     % HEIGHT OF CORE

%------ graphite hexagon for radial graphite reflector ------
surf 301 hexxc 0   0   {rs*hex}      % HEX FOR REFLECTOR GRAPHITE
surf 398 pz    0	             % BOTTOM OF CORE
surf 399 pz    {zcore}		     % HEIGHT OF CORE

%------ fuel channel cap and holding shaft ------
surf 401 cyl   0   0   {r3}          % FUEL CAP ABOVE FUEL CHANNELS
surf 402 hexxc 0   0   {hexg}        % GRAPHITE HEX AROUND FUEL CAP
surf 402 cyl   0   0   {r_shaft}     % CYLINDRICAL HOLDING SHAFT
surf 496 pz    {zcore}		     % TOP OF FUEL IN CORE
surf 497 pz    {z_cap_fuel}          % TOP OF INNER CAP
surf 498 pz    {z_cap_graphite}      % TOP OF GRAPHITE IN CAP
surf 499 pz    {z_shaft}             % TOP OF HOLDING SHAFT

%------ top graphite reflector ------
surf 501 cyl   0   0   {topr_r}      % TOP AXIAL GRAPHITE REFLECTOR
surf 598 pz    {z_topr_bot}          % BOTTOM OF REFELECTOR
surf 599 pz    {z_topr_top}          % TOP OF REFLECTOR

%------ holding plate ------
surf 601 cyl   0   0   {r_plate}     % TOP HOLDING PLATE
surf 698 pz    {z_cap_graphite}      % BOTTOM OF HOLDING PLATE
surf 699 pz    {z_plate}             % TOP OF HOLDING PLATE

'''




	surfaces = surfaces.format(**locals())
	
	return surfaces


if __name__ == '__main__':
	print "This is a module to write surfaces for the MSR core."
	raw_input("Press Ctrl+C to quit, or enter else to test it. ")
	print write_surfs(pitch = 11.5, slit = 0.323, d=6.267, ri = 2.2, ro = 0.9, \
					rox = 3.54, roy = 2.05, c = 0.65, \
					rfuel = 400, rcore = 427, rgref = 450, rhast = 470, \
					zcore = 100, pht = 10, zrefl = 35)
