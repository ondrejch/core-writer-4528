#!/usr/bin/env python3
#
# Cells: Write the cell cards and universes for the Serpent input deck


def write_cells(universes=range(1,1+12), lattices=range(33,33+7), \
				surffuel=30, surfcore=31, surfgref=32, surfhast=29):
	'''Function to write cell cards for Serpent input deck.
	Inputs:
		universes:		Tuple(12) containing the following.
		 ub, uf, uc:	 Universe number of blanket, fuel, central cells
		 uup, ul#:		 Universe number of upper/lower plenum assembly
		 uuc, ulc:		 Universe number of upper/lower control cells
		lattices:		Tuple(7) containing the following:
		 latmid:	 	 Lattice number of the middle fuel cells
		 lattop:		 Lattice number of the top fuel cells
		 latbot#: 		 Lattice number of the bottom fuel cells
		surffuel:		Surface no. of inner core (fuel cells)
		surfcore:		Surface no. of entire core, fuel+blanket
		surfgref:		Surface no. of graphite reflector
	Outputs:
		cells:		String containing cell cards'''
	

	# Unpack the universe tuple
	ub, uf, uc, uup, \
	ul1, ul2, ul3, ul4, \
	ulp, uuc, ulc, uh 		 = universes
	# Unpack the lattice tuple
	latmid, lattop, latbot1, latbot2, latbot3, \
	latbot4, latplen = lattices
	
	'''cells = \'''
%------define cells--------------------

% Universe {uf}: FUEL CELL U2
cell 10 {uf} graphite -10 20 22
cell 13 {uf} graphite 20 -21
cell 11 {uf} blanket  -11 10 
%cell 12 {uf} outside   11 
cell 12 {uf} void  11 
cell  7 {uf} fuel -20 
cell  1 {uf} fuel 21 -22 

% Universe 7: Blank Blanket Cell
cell 701 7 blanket -200
cell 702 7 blanket 200

% Universe 8: Holding Shafts on top of plate
cell 801 8 graphite -202
cell 802 8 blanket  202

% Universe 9: Holding Shafts under plate
cell 901 9 graphite -206
cell 902 9 blanket 206

% Universe 12: Holding Plate
cell 1201 12 graphite -207
cell 1202 12 hastelloy 207

% Universe {ub}: BLANKET CELL
cell 15 {ub} graphite -101 %27 -28
cell 16 {ub} blanket     101 %27 -28

% Universe {uh}: HASTELLOY HEX
cell 18 {uh} hastelloy -11
cell 19 {uh} void       11

% Universe {uc}: CONTROL ROD 
% Similar to fuel cell, but with helium in channels
cell 30 {uc} graphite -10 20 90 91 92 93 94 95  
cell 41 {uc} blanket  -11 10 
cell 401 {uc} blanket -20 
cell 42 {uc} void   11 %27 -28
%cell 38 {uc} void  -27
%cell 39 {uc} void   28
cell  34 {uc} he -90 %27 -28
cell  35 {uc} he -91 
cell  36 {uc} he -92 
cell  37 {uc} he -93 
cell  38 {uc} he -94 
cell  39 {uc} he -95 





% Universe {uup}: UPPER CHANNEL
cell 43 {uup} fuel     -22    	28 -41	% fuel cap
cell 44 {uup} graphite  22 -10 	28 -41	% graphite,level1
cell 45 {uup} graphite -10      41 -42	% graphite cap,level2
cell 46 {uup} blanket  -11 10 	28 -42  % slit all the way through
cell 47 {uup} void      42
cell 48 {uup} void     -28
cell 49 {uup} void      11 		28 -42

% Universe {ul1}: LOWER CHANNEL 1
cell 51 {ul1} blanket -11 3503 
cell 52 {ul1} fuel -20    % -27 56  % Central channel
cell 53 {ul1} hastelloy -53  20 % -27 56  % Hast. pipe
cell 58 {ul1} graphite -3501 53
cell 54 {ul1} fuel -3502  3501%53 % 2501  % outer fuel channel
cell 55 {ul1} graphite -3503 3502 % 3501 was 55 2501 % graphite hex
cell 59 {ul1} void 11

% Universe {ul2}: LOWER CHANNEL 2
cell 251 {ul2} blanket -11 57 
cell 252 {ul2} fuel -20    % -27 56  % Central channel
cell 253 {ul2} hastelloy -53  20 % -27 56  % Hast. pipe
cell 254 {ul2} fuel -2501  53 % -27 56  % outer fuel channel
cell 257 {ul2} hastelloy -57 2501 % -27 56 % outer pipe
cell 259 {ul2} void 11

% Universe 1111: LOWER PLENUM TOP 
cell 2611 1111 fuel -11 53
cell 2621 1111 fuel -20          % Central channel
cell 2631 1111 hastelloy -53  20 % Inner pipe
cell 2651 1111 void 11

% Universe {ulp}: LOWER PLENUM BOTTOM
cell 250 {ulp} fuel -11
cell 260 {ulp} void  11

% Universe {ul3}: PENENETRATION TO INLET PLENUM
cell 261 {ul3} hastelloy -11 2501
cell 262 {ul3} fuel -20          % Central channel
cell 263 {ul3} hastelloy -53  20 % Inner pipe
cell 264 {ul3} fuel -2501  53      % outer fuel channel
cell 265 {ul3} void 11
% Universe {ul4}: PENETRATION TO OUTLET PLENUM
cell 266 {ul4} hastelloy -11 20
cell 267 {ul4} fuel -20
cell 268 {ul4} void 11

% Universe {uuc}: UPPER CONTROL
cell 61 {uuc} he       -12    	28 -41	% helium gap
cell 62 {uuc} graphite  12 -10 	28 -42	% graphite,level1
cell 63 {uuc} graphite -12      41 -42	% graphite cap,level2
cell 64 {uuc} blanket  -11 10 	28 -42  % slit all the way through
cell 65 {uuc} void      42
cell 66 {uuc} void     -28
cell 67 {uuc} void      11 		28 -42

% Universe {ulc}: LOWER CONTROL
cell 71 {ulc} graphite -20               % Central channel - was he
cell 72 {ulc} graphite -53  20    % Hast. pipe - was hastelelloy
cell 73 {ulc} graphite -2501  53           % outer fuel channel - was he
cell 74 {ulc} graphite -55 2501      % graphite hex
cell 75 {ulc} blanket  -11 55      % blanket reflector
cell 79 {ulc} void  11 


% The main universe
%cell 100 0 fill    {latmid} -{surffuel}     27 -28
cell 100 0 fill    {latmid} -{surfgref}     27 -28
%cell 101 0 blanket {surffuel} -{surfcore}   61 -42 
%cell 102 0 graphite {surfcore} -{surfgref}  61 -81
cell 104 0 fill   {lattop} -{surfgref}      28 -42
cell 105 0 fill   {latbot1} -{surfgref}     52 -27
cell 106 0 fill   {latbot2} -{surfgref}     56 -52
cell 107 0 fill   {latbot3} -{surfgref}	  	61 -56
cell 108 0 fill   43 -{surffuel}	  	62 -61
cell 109 0 fill   {latbot4} -{surffuel}	  	63 -62
cell 110 0 fill   {latplen} -{surffuel}	  	64 -63
cell 111 0 hastelloy       -{surfgref}		60 -64
cell 112 0 hastelloy {surffuel} -{surfgref} 64 -61
cell 113 0 hastelloy {surfgref} -{surfhast} 60 -81
%cell 120 0 blanket  	 	    -{surfgref} 42 -81 102 %-80  -{surfcore}
%cell 121 0 blanket 	 	    -{surfcore} 80 -81  % was top reflector
cell 122 0 hastelloy	 	    -{surfhast} 81 -82
cell 123 0 fill 41 -{surfgref} 208 -209% top graphite reflector
cell 124 0 fill 42 -{surfgref} 209 -203 % top holding plate
cell 125 0 fill 40 -{surfgref} 203 -204
cell 126 0 blanket -205 204 -81


cell 999 0 outside {surfhast} 60 -82
cell 998 0 outside -60
cell 997 0 outside  82
	'''



# New cell definitions - WARNING!!! THESE ARE NOT USED IN THE CURRENT VERSION
	cells = '''
%------define cells--------------------
% Universe {ub}: BLANKET CELL U1
cell 15 {ub} graphite   -101
cell 16 {ub} blanket     101 

	
% Universe {uf}: FUEL CELL U2
cell 10 {uf} graphite -201 203 205
cell 13 {uf} graphite 203 -204
cell 11 {uf} blanket  -202 201  
cell 12 {uf} void  202 
cell  7 {uf} fuel -203 
cell  1 {uf} fuel 204 -205 


% Universe {uc}: CONTROL ROD U3 
% Similar to fuel cell, but with helium in channels
cell 30 {uc} graphite -301 303 304 305 306 307 308 309  
cell 41 {uc} blanket  -302 301 
cell 401 {uc} blanket -303 
cell 42 {uc} void   302 
cell  34 {uc} he -304 
cell  35 {uc} he -305 
cell  36 {uc} he -306 
cell  37 {uc} he -307 
cell  38 {uc} he -308 
cell  39 {uc} he -309 


% Universe {uup}: UPPER CHANNEL U4
cell 43 {uup} fuel     -403  404 -405	    % fuel cap
cell 44 {uup} graphite  403 -401  404 -405	% graphite,level1
cell 45 {uup} graphite -401  405 -406	    % graphite cap,level2
cell 46 {uup} blanket  -402  401  404 -406  % slit all the way through
cell 47 {uup} void      406
cell 48 {uup} void     -404
cell 49 {uup} void      402  404 -406


% Universe {uuc}: UPPER CONTROL U5
cell 61 {uuc} he       -503  504 -505	    % helium gap
cell 62 {uuc} graphite  503 -501  504 -506	% graphite,level1
cell 63 {uuc} graphite -503  505 -506	    % graphite cap,level2
cell 64 {uuc} blanket  -502  501  504 -506  % slit all the way through
cell 65 {uuc} void      506
cell 66 {uuc} void     -504
cell 67 {uuc} void      502  504 -506


% Universe {ulc}: LOWER CONTROL U6
cell 71 {ulc} graphite -602       % Central channel - was he
cell 72 {ulc} graphite -603  602  % Hast. pipe - was hastelelloy
cell 73 {ulc} graphite -604  603  % outer fuel channel - was he
cell 74 {ulc} graphite -605  604  % graphite hex
cell 75 {ulc} blanket  -601  605  % blanket reflector
cell 79 {ulc} void      601 


% Universe 7: Blank Blanket Cell U7
cell 701 7 blanket -701
cell 702 7 blanket 701


% Universe 8: Holding Shafts on top of plate U8
cell 801 8 graphite -801
cell 802 8 blanket  801


% Universe 9: Holding Shafts under plate U9
cell 901 9 graphite -901
cell 902 9 blanket 901


% Universe {ulp}: LOWER PLENUM BOTTOM U10
cell 250 {ulp} fuel -1001
cell 260 {ulp} void  1001


% Universe {uh}: HASTELLOY HEX U11
cell 18 {uh} hastelloy -1101
cell 19 {uh} void       1101


% Universe 12: Holding Plate U12
cell 1201 12 graphite -1201
cell 1202 12 hastelloy 1201

% Universe {ul1}: LOWER CHANNEL 1 U25
cell 51 {ul1} blanket   -2501  2507 
cell 52 {ul1} fuel      -2502        % Central channel
cell 53 {ul1} hastelloy -2503  2502  % Hast. pipe
cell 58 {ul1} graphite  -2505  2503
cell 54 {ul1} fuel      -2506  2505  % outer fuel channel
cell 55 {ul1} graphite  -2507  2506  % graphite hex
cell 59 {ul1} void       2501

% Universe {ul2}: LOWER CHANNEL 2  U26
cell 251 {ul2} blanket   -2601  2605 
cell 252 {ul2} fuel      -2602        % Central channel
cell 253 {ul2} hastelloy -2603  2602  % Hast. pipe
cell 254 {ul2} fuel      -2604  2603  % outer fuel channel
cell 257 {ul2} hastelloy -2605  2604  % outer pipe
cell 259 {ul2} void       2601



% Universe {ul3}: PENENETRATION TO INLET PLENUM U27
cell 261 {ul3} hastelloy -2701  2703
cell 262 {ul3} fuel      -2702        % Central channel
cell 263 {ul3} hastelloy -2704  2702  % Inner pipe
cell 264 {ul3} fuel      -2703  2704  % outer fuel channel
cell 265 {ul3} void       2701


% Universe {ul4}: PENETRATION TO OUTLET PLENUM U28
cell 266 {ul4} hastelloy -2801  2802
cell 267 {ul4} fuel      -2802
cell 268 {ul4} void       2801


% Universe 1111: LOWER PLENUM TOP U1111
cell 2611 1111 fuel      -111101  111103
cell 2621 1111 fuel      -111102          % Central channel
cell 2631 1111 hastelloy -111103  111102  % Inner pipe
cell 2651 1111 void       111101


% The main universe
cell 100 0 fill       {latmid}   -{surfgref} 8 -9
cell 104 0 fill       {lattop}   -{surfgref} 9 -10
cell 105 0 fill       {latbot1}  -{surfgref} 7 -8
cell 106 0 fill       {latbot2}  -{surfgref} 6 -7
cell 107 0 fill       {latbot3}  -{surfgref} 5 -6
cell 108 0 fill   43 -{surffuel}  4 -5
cell 109 0 fill       {latbot4}  -{surffuel} 3 -4
cell 110 0 fill       {latplen}  -{surffuel} 2 -3
cell 111 0 hastelloy -{surfgref}		1 -2
cell 112 0 hastelloy  {surffuel} -{surfgref} 2 -5
cell 113 0 hastelloy  {surfgref} -{surfhast} 1 -15
cell 122 0 hastelloy -{surfhast} 15 -16
cell 123 0 fill   41 -{surfgref} 11 -12 % holding shaft below plate
cell 124 0 fill   42 -{surfgref} 12 -13 % top holding plate
cell 125 0 fill   40 -{surfgref} 13 -14
cell 126 0 blanket -17 14 -15


cell 999 0 outside {surfhast} 1 -16
cell 998 0 outside -1
cell 997 0 outside  16
	'''

	cells = cells.format(**locals())
	return cells

if __name__ == '__main__':
	print("This is a module to write cells for the MSR core.")
	input("Press Ctrl+C to quit, or enter else to test it. ")
	print (write_cells())
