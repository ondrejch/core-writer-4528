#! /usr/bin/python
#
# Surfs
# Calculate dimensions for each of the surfaces based off the MSBR cell dimensions
# and output the surface cards for the Serpent input deck.


import math


def write_surfs(fsf, pitch, slit, ro, r2, rs, rfuel, rcore_inner, rcore_outer, zcore, pht, zrefl):
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
	dgr = 15	
	plenum_vol = 37*28316.8 	# 37 ft^3 to cm^2
	# Height of each plenum: inlet and outlet
	plenum_ht = plenum_vol / (2*math.pi*rfuel**2)
	gt = 6*2.54 # thickness of graphite: cm
	ht = 3*2.54 # thickness of hastelloy, placeholder
	rgref = rcore_inner + gt
	rhast = rgref + ht
	# Calculate the actual dimensions based on our parameters
	# Have benchmarked this against existing perl script--should be right
	#
	# inradius: half the pitch
	hpitch = pitch/2.0
	# circumradius: distance from center to corner
	d = (hpitch - slit) * 2.0/math.sqrt(3)  
	# radius (inner): central fuel channel radius
	hexarea = 2.0 * math.sqrt(3.0) * 10**2
	#r1 = math.sqrt(hexarea*fsf/(2.0*math.pi) )
	r1 = 2.06
	# radius (outer): auxiliary fuel channel radius
#	ro = ri / math.sqrt(6)
	ro = 1.1
	# c: a constant that determines how far along the circumradius the channels appear
	c = (r1 + d*math.sqrt(3)/2.0) / ( d*(1.0 + math.sqrt(3)/2.0) )
	# X and Y coordinates of ro
	rox = c*d*math.sqrt(3)/2.0
	roy = c*d*0.5
	
	# Radial reflector scaling term
	rs = 0.9
 # r2 is the outer fuel radius; thast = hastelloy thickness (1/8 in)
	thast = 1.0/8 * 2.54							# hastelloy thickness
	#r2 = math.sqrt(2*r1**2 + 2*r1*thast + thast**2) # outer fuel cylinder
	r2 = r1 + 1.127
	# Radius of outer fuel ring with equal volume to inner fuel channel
	r3 = math.sqrt(r1**2 + r2**2)
	rdiff = (r3 - r2)
	# Establish a few additional dimensions
	hexs = pitch/2.0    # radius of cell, outside slit
	#ry = c*d
	#hexg = hexs - slit   # radius of graphite, inside slit
	hexg = 6.83
	ry = c*d			# y coord of vertical channel

	# hexf: top channel hexagon
	hexf = math.sqrt(3)/2.0*(d*c) + rdiff
	zhexf1 = zcore + 2.0*rdiff	# top of fuel channel (salt)
	# Make the "roof" as thick as the "wall" of the graphite hexagon
	zhexf2 = zhexf1 + (hexg - hexf) # top of fuel channel (graphite)
	# Top of the core: Axial reflector
	zrefltop = zhexf2 + zrefl
	zgreftop = zrefltop + (rgref - rcore_inner)
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

    # radial reflector 
	rs = 0.9 
	rr = rs*hexg
	axial_top = zhexf2+dgr
	zplate = axial_top + 15
	zshaft = zplate + 30

	surfaces = '''
%------define the hexagon and fuel channel cells----
surf 10 hexxc 0   0   {hexg}	     % HEX FOR GRAPHITE 6.82625cm for def
surf 11 hexxc 0   0   {hexs}	      % HEX FOR SLIT
surf 12 hexxc 0   0   {hexf}          % HEX FOR FUEL
surf 20 cyl   0   0   {r1}	    % CENTER HOLE
surf 21 cyl   0   0   {r2}     % INTERMEDIATE GRAPHITE RING
surf 22 cyl   0   0   {r3}     % OUTER FUEL RING
surf 27 pz    0			    % BOTTOM OF CORE
surf 28 pz    {zcore}		     % HEIGHT OF CORE
surf 30 cyl   0   0   {rfuel}
surf 31 cyl   0   0   {rcore_inner}
surf 32 cyl   0   0   {rgref}
surf 33 cyl   0   0   {rhast}
surf 41 pz    {zhexf1}
surf 42 pz    {zhexf2}
surf 51 cone  0 0 0 {hexg}  -{hexg}
surf 52 pz    {ztrans1}
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
surf 101 hexxc 0      0   {rr}      % HEX FOR RADIAL REFLECTOR
surf 102 cylz 0 0 {rgref} {zhexf2} {axial_top} % top axial reflector from top of the upper channel
surf 200 hexxc 0 0 {hexg}
surf 201 cyl 0 0 {rgref} {axial_top} {zplate} % top holding plate
surf 202 cyl 0 0 {r1}     % holding shaft
surf 203 pz {zplate}   % top of holding plate
surf 204 pz {zshaft}   % top of holding shafts
surf 205 cyl 0 0 {rgref} % blanket above core
'''


# New calculations for surface definitions commented sections are a work in progress
	''' Variables that need to be given to the function:
	ro                    # Radius of the center hole for the centeral cell
	r2                    # Radius of the concentric graphtie ring
	rs                    # scaling factor for radial reflector hexagons
	rfuel                 # Radius of the fuel lattice
	rcore_inner           # Radius of the inner core barrel
	rcore_outer           # Radius of the outer core barrel
	zcore                 # height of core
	slit                  # Width of slit 
	pitch                 # Lattice pitch
	fsf                   # fuel salt fraction
	pht                   # height of each of the lower plena
	zrefl                 # height of the axial reflector
	'''

	'''plenum_vol = 37*28316.8 	# 37 ft^3 to cm^2

# Height of each plenum: inlet and outlet
plenum_ht = plenum_vol / (2*math.pi*rfuel**2)

# thast = hastelloy thickness (1/8 in)
thast = 1.0/8 * 2.54							# hastelloy thickness

# radius (inner): central fuel channel radius
hexarea = 2.0 * math.sqrt(3.0) * 10**2
r1 = math.sqrt(hexarea*fsf/(2.0*math.pi) )

# Radius of outer fuel ring with equal volume to inner fuel channel
r3 = math.sqrt(r1**2 + r2**2)
rdiff = (r3 - r2)

# Establish a few additional dimensions
hexs = pitch/2.0     # radius of cell, outside slit

# circumradius: distance from center to corner
d = (pitch/2.0 - slit) * 2.0/math.sqrt(3)

# c: a constant that determines how far along the circumradius the channels appear
c = (ri + d*math.sqrt(3)/2.0) / ( d*(1.0 + math.sqrt(3)/2.0) )

hexg = hexs - slit   # radius of graphite, inside slit
ry = c*d			 # y coord of vertical channel

# X and Y coordinates of control rods
ro_x = c*d*math.sqrt(3)/2.0
ro_y = c*d*0.5

# hexf: top channel hexagon
hexf = math.sqrt(3)/2.0*(d*c) + rdiff
z_cap_fuel = zcore + 2.0*rdiff	# top of fuel channel (salt)

# Make the "roof" as thick as the "wall" of the graphite hexagon
z_cap_graphite = z_cap_fuel + (hexg - hexf) # top of fuel channel (graphite)

# Top of the core: Axial reflector
z_topr_top = z_cap_graphite + zrefl
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

# New surface definitions
surfaces = \'''
%------ main universe ------
surf 1 cyl   0   0   {rfuel}         % CYLINDRICAL BOUNDS FOR FUEL LATTICE
surf 2 cyl   0   0   {rcore_inner}   % CYLINDRICAL BOUNDS FOR FUEL AND REFLECTOR
surf 3 cyl   0   0   {rcore_outer}   % CYLINDRICAL BOUNDS FOR ENTIRE CORE
surf 4 pz    {zbot}                  % LOWER PLENUM BOUNDARY PARAMETERS (4-11)
surf 5 pz    {zitop}
surf 6 pz    {zibot}
surf 7 pz    {zotop}
surf 8 pz    {zobot}
surf 9 pz    {zrefltop}
surf 10 pz    {zgreftop}
surf 11 pz    {zhasttop}

%------ graphite hexagon and fuel channel cells ------	
surf 101 hexxc 0   0   {hexg}	     % HEX FOR GRAPHITE
surf 102 hexxc 0   0   {hexs}	     % HEX FOR SLIT
surf 103 cyl   0   0   {r1}	         % CENTER HOLE
surf 104 cyl   0   0   {r2}          % INTERMEDIATE GRAPHITE RING
surf 105 cyl   0   0   {r3}          % OUTER FUEL RING
surf 198 pz    0		             % BOTTOM OF CORE
surf 199 pz    {zcore}		         % HEIGHT OF CORE

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
surf 298 pz    0		             % BOTTOM OF CORE
surf 299 pz    {zcore}		         % HEIGHT OF CORE

%------ graphite hexagon for radial graphite reflector ------
surf 301 hexxc 0   0   {rs*hex}      % HEX FOR REFLECTOR GRAPHITE
surf 398 pz    0	                 % BOTTOM OF CORE
surf 399 pz    {zcore}		         % HEIGHT OF CORE

%------ fuel channel cap and holding shaft ------
surf 401 cyl   0   0   {r3}          % FUEL CAP ABOVE FUEL CHANNELS
surf 402 hexxc 0   0   {hexg}        % GRAPHITE HEX AROUND FUEL CAP
surf 402 cyl   0   0   {r_shaft}     % CYLINDRICAL HOLDING SHAFT
surf 403 hexxc 0   0   {hexs}	     % HEX FOR SLIT
surf 496 pz    {zcore}		         % TOP OF FUEL IN CORE
surf 497 pz    {z_cap_fuel}          % TOP OF INNER CAP
surf 498 pz    {z_cap_graphite}      % TOP OF GRAPHITE IN CAP
surf 499 pz    {z_shaft}             % TOP OF HOLDING SHAFT

%------ top graphite reflector ------
surf 501 cyl   0   0   {rfuel}      % TOP AXIAL GRAPHITE REFLECTOR
surf 598 pz    {z_cap_graphite}     % BOTTOM OF REFELECTOR
surf 599 pz    {z_topr_top}         % TOP OF REFLECTOR

%------ holding plate ------
surf 601 cyl   0   0   {r_plate}     % TOP HOLDING PLATE
surf 698 pz    {z_cap_graphite}      % BOTTOM OF HOLDING PLATE
surf 699 pz    {z_plate}             % TOP OF HOLDING PLATE

%------ lower channel 1 ------
surf 701 hexxc 0   0   {hexs}      % HEX FOR SLIT
surf 702 hexxc 0   0   {rg}	       % HEX FOR FUEL TRANSITION
surf 703 cyl   0   0   {r1}	       % CENTER HOLE
surf 704 cyl   0   0   {rh}        % HASTELLOY PIPE RADIUS
surf 705 cyl   0   0   {r2}        % OUTER FUEL CYL RADIUS

%------ lower channel 2 ------
surf 801 hexxc 0   0   {hexs}	   % HEX FOR SLIT
surf 802 cyl   0   0   {r1}	       % CENTER HOLE
surf 803 cyl   0   0   {rh}        % HASTELLOY PIPE RADIUS
surf 804 cyl   0   0   {r2}        % OUTER FUEL CYL RADIUS
surf 805 cyl   0   0   {rh2}       % OUTER HASTELLOY PIPE

%------ lower plenum ------
surf 901 cyl   0   0   {r1}	   % CENTER HOLE

%------ penetration to inlet plenum ------
surf 1001 hexxc 0   0   {hexs}	   % HEX FOR SLIT
surf 1002 cyl   0   0   {r2}       % OUTER FUEL CYL RADIUS
surf 1003 cyl   0   0   {r1}	   % CENTER HOLE
surf 1004 cyl   0   0   {rh}       % HASTELLOY PIPE RADIUS

%------ penetration to outlet plenum ------
surf 1101 hexxc 0   0   {hexs}	   % HEX FOR SLIT
surf 1102 cyl   0   0   {r1}	   % CENTER HOLE

	'''




	surfaces = surfaces.format(**locals())
	
	return surfaces


if __name__ == '__main__':
	print "This is a module to write surfaces for the MSR core."
	raw_input("Press Ctrl+C to quit, or enter else to test it. ")
	print write_surfs(pitch = 11.5, slit = 0.323, d=6.267, r1 = 2.2, ro = 0.9, \
					rox = 3.54, roy = 2.05, c = 0.65, \
					rfuel = 400, rcore = 427, rgref = 450, rhast = 470, \
					zcore = 100, pht = 10, zrefl = 35)
