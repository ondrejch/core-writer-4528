#! /usr/bin/python
#
# Surfs
# Calculate dimensions for each of the surfaces based off the MSBR cell dimensions
# and output the surface cards for the Serpent input deck.


import math


def write_surfs(fsf, relba, pitch, slit,temp, ro, r2, rs, rfuel, rcore_inner, rcore_outer, zcore, pht, zrefl):
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
		gr_exp: graphite expansion coefficient m/m K
	Output:
		surfaces:   string containing the surface cards for the MSR'''
	gr_exp = 4.14*10**(-6)
	l=pitch/2.0
	dgr = 15	
	plenum_vol = 37*28316.8 	# 37 ft^3 to cm^2
	# Height of each plenum: inlet and outlet
	plenum_ht = plenum_vol / (2*math.pi*rfuel**2)
	gt = 6*2.54 # thickness of graphite: cm
	ht = 2*2.54 # thickness of hastelloy, placeholder
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
	hexarea = 2.0 * math.sqrt(3.0) * l**2
	r1 = math.sqrt(hexarea*fsf/(2.0*math.pi) )
	#r1 = 2.06
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
	#r2 = r1 + 1.127
	# Radius of outer fuel ring with equal volume to inner fuel channel
	r3 = math.sqrt(r1**2 + r2**2)
	rdiff = (r3 - r2)
	# Establish a few additional dimensions
	hexs = l   # radius of cell, outside slit
	#ry = c*d
	blanketfraction = 1.06923
	blanketA0 = blanketfraction * r1**2 *math.pi
	blanketarea = blanketA0 * relba
	l2 = math.sqrt( l**2 - blanketarea / (2.0 * math.sqrt(3.0)))
	hexg = l2  + (temp - 700)*gr_exp*l2 # radius of graphite, inside slit with thermal expansion 700C nominal temp
	#hexg = 6.83
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
	r3_bot = math.sqrt(r1**2 + rh**2)               # outer fuel radius in lower plenum
	rg = (gt1/2.0 + r3_bot)    				   			# radius of graphite hex
	# And cut off the bottom 10 inches below that
	ztrans1 = 0 - 10*2.54
	# Next layer down, the concentric circles w/ hastelloy
	rh2 = r3_bot + thast		# Outer hastelloy pipe
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
	rr = rs*hexg  # reflector size
	hsl = 3.97  # holding shaft lower radius
	hsu = hsl*0.8     # holding shaft upper radius
	axial_top = zhexf2+30.48
	zplate = axial_top + 12
	zshaft = zplate + 5 
	
	# lower plenum cone calculations: naming scheme goes outward from center
	conez1 = ((0-ztrans1)/(r2-rh))*r2      
	conez2 = ((0-ztrans1)/(r3-r3_bot))*r3
	conez3 = ((0-ztrans1)/(hexg-rh2))*hexs
	
	# Old definitions, currently not in use
	'''surfaces = \'''
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
surf 202 cyl 0 0 {hsu}     % holding shaft
surf 203 pz {zplate}   % top of holding plate
surf 204 pz {zshaft}   % top of holding shafts
surf 205 cyl 0 0 {rgref} % blanket above core
surf 206 cyl 0 0 {hsl}
surf 207 cylz 0 0 {hsu}
surf 208 pz {zhexf2}
surf 209 pz {axial_top}
surf 2501 cyl 0  0 {r3_bot} % outer fuel radius for lower plenum
surf 3501 cone 0 0 0 {r2} -{conez1}
surf 3502 cone 0 0 0 {r3} -{conez2}
surf 3503 cone 0 0 0 {hexg} -{conez3}
'''


# New calculations for surface definitions: comment block below is not up to date
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



# New surface definitions
	surfaces = '''
%------ main universe ------
%% vertical core plane divisions

surf 16 pz    {zhasttop}
surf 15 pz    {zgreftop}
surf 14 pz {zshaft}   % top of holding shafts
surf 13 pz {zplate}   % top of holding plate
surf 12 pz {axial_top}
surf 11 pz {zhexf2}
surf 10 pz    {zhexf2}
surf 9 pz    {zcore}		     % HEIGHT OF CORE
surf 8 pz    0			    % BOTTOM OF CORE
surf 7 pz    {ztrans1}
surf 6 pz    {ztrans2}
surf 5 pz    {zitop}
surf 4 pz    {zibot}
surf 3 pz    {zotop}
surf 2 pz    {zobot}
surf 1 pz    {zbot}

%% radial bounds
surf 17 cyl 0 0 {rgref} % blanket above core
surf 18 cyl   0   0   {rfuel}
surf 19 cyl   0   0   {rgref}
surf 20 cyl   0   0   {rhast}


%------ blanket cell Universe 1 ------
surf 101 hexxc 0   0   {rr}          % HEX FOR REFLECTOR GRAPHITE


%------ graphite hexagon and fuel channel cells Universe 2 ------	
surf 201 hexxc 0   0   {hexg}	     % HEX FOR GRAPHITE
surf 202 hexxc 0   0   {hexs}	     % HEX FOR SLIT
surf 203 cyl   0   0   {r1}	         % CENTER HOLE
surf 204 cyl   0   0   {r2}          % INTERMEDIATE GRAPHITE RING
surf 205 cyl   0   0   {r3}          % OUTER FUEL RING


%------ graphite hexagon and control rod channels Universe 3 ------
surf 301 hexxc 0   0   {hexg}	     % HEX FOR GRAPHITE
surf 302 hexxc 0   0   {hexs}	     % HEX FOR SLIT
surf 303 cyl   0   0   {r1}  	     % CENTER HOLE
surf 304 cyl   0   {ry}   {ro}	     % OUTER HOLES x 6 CONTROL RODS
surf 305 cyl   0  -{ry}   {ro}       %           ||
surf 306 cyl   {rox}  {roy}  {ro}    %           ||
surf 307 cyl  -{rox}  {roy}  {ro}    %          _||_
surf 308 cyl  -{rox} -{roy}  {ro}    %          \  /
surf 309 cyl   {rox} -{roy}  {ro}    % __________\/____________


%------ upper channel Universe 4 ------
surf 401 hexxc 0   0   {hexg}	      % HEX FOR GRAPHITE
surf 402 hexxc 0   0   {hexs}	      % HEX FOR SLIT
surf 403 cyl   0   0   {r3}           % OUTER FUEL RING
surf 404 pz    {zcore}		          % HEIGHT OF CORE
surf 405 pz    {zhexf1}
surf 406 pz    {zhexf2}


%------ upper control Universe 5 ------
surf 501 hexxc 0   0   {hexg}	     % HEX FOR GRAPHITE
surf 502 hexxc 0   0   {hexs}	      % HEX FOR SLIT
surf 503 hexxc 0   0   {hexf}          % HEX FOR FUEL
surf 504 pz    {zcore}		     % HEIGHT OF CORE
surf 505 pz    {zhexf1}
surf 506 pz    {zhexf2}


%------ lower control Universe 6 ------
surf 601 hexxc 0   0   {hexs}	      % HEX FOR SLIT
surf 602 cyl   0   0   {r1}	    % CENTER HOLE
surf 603 cyl   0   0   {rh}   % hastelloy tube radius
surf 604 cyl 0  0 {r3_bot} % outer fuel radius for lower plenum
surf 605 hexxc 0   0   {rg}	% Hex for fuel transition


%------ blank blanket cell Universe 7 ------
surf 701 hexxc 0 0 {hexg}


%------ holding shafts on top of plate Universe 8 ------
surf 801 cyl 0 0 {hsu}     % holding shaft


%------ holding shafts under plate Universe 9 ------
surf 901 cyl 0 0 {hsl}


%------ hastelloy hex Universe 11 ------
surf 1101 hexxc 0   0   {hexs}	      % HEX FOR SLIT


%------ lower plenum bottom Universe 10 ------
surf 1001 hexxc 0   0   {hexs}	      % HEX FOR SLIT


%------ holding plate Universe 12 ------
surf 1201 cylz 0 0 {hsu}


%------ lower channel 1 Universe 25 ------
surf 2501 hexxc 0   0   {hexs}      % HEX FOR SLIT
surf 2502 cyl   0   0   {r1}	       % CENTER HOLE
surf 2503 cyl   0   0   {rh}        % HASTELLOY PIPE RADIUS
surf 2504 cyl   0   0   {r2}        % OUTER FUEL CYL RADIUS
surf 2505 cone 0 0 0 {r2} -{conez1}
surf 2506 cone 0 0 0 {r3} -{conez2}
surf 2507 cone 0 0 0 {hexg} -{conez3}


%------ lower channel 2 Universe 26------
surf 2601 hexxc 0   0   {hexs}	   % HEX FOR SLIT
surf 2602 cyl   0   0   {r1}	   % CENTER HOLE
surf 2603 cyl   0   0   {rh}       % HASTELLOY PIPE RADIUS
surf 2604 cyl   0   0   {r3_bot}   % outer fuel radius for lower plenum
surf 2605 cyl   0   0   {rh2}      % OUTER HASTELLOY PIPE


%------ penetration to inlet plenum Universe 27 ------
surf 2701 hexxc 0   0   {hexs}	      % HEX FOR SLIT
surf 2702 cyl   0   0   {r1}	    % CENTER HOLE
surf 2703 cyl 0  0 {r3_bot} % outer fuel radius for lower plenum
surf 2704 cyl   0   0   {rh}   % hastelloy tube radius


%------ penetration to outlet plenum Universe 28 ------
surf 2801 hexxc 0   0   {hexs}	      % HEX FOR SLIT
surf 2802 cyl   0   0   {r1}	    % CENTER HOLE


%------ lower plenum top Universe 1111 ------
surf 111101 hexxc 0   0   {hexs}	      % HEX FOR SLIT
surf 111102 cyl   0   0   {r1}	    % CENTER HOLE
surf 111103 cyl   0   0   {rh}   % hastelloy tube radius

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
