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
	
	cells = '''
%------define cells--------------------

% Universe {uf}: FUEL CELL
cell 10 {uf} graphite -10 20 22
cell 13 {uf} graphite 20 -21
cell 11 {uf} blanket  -11 10 
%cell 12 {uf} outside   11 
cell 12 {uf} void  11 
cell  7 {uf} fuel -20 
cell  1 {uf} fuel 21 -22 



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
cell 43 {uup} fuel     -12    	28 -41	% fuel cap
cell 44 {uup} graphite  12 -10 	28 -42	% graphite,level1
cell 45 {uup} graphite -12      41 -42	% graphite cap,level2
cell 46 {uup} blanket  -11 10 	28 -42  % slit all the way through
cell 47 {uup} void      42
cell 48 {uup} void     -28
cell 49 {uup} void      11 		28 -42

% Universe {ul1}: LOWER CHANNEL 1
cell 51 {ul1} blanket -11 55 
cell 52 {ul1} fuel -20    % -27 56  % Central channel
cell 53 {ul1} hastelloy -53  20 % -27 56  % Hast. pipe
cell 54 {ul1} fuel -54  53 % -27 56  % outer fuel channel
cell 55 {ul1} graphite -55 54 % 52 % graphite hex
cell 59 {ul1} void 11

% Universe {ul2}: LOWER CHANNEL 2
cell 251 {ul2} blanket -11 57 
cell 252 {ul2} fuel -20    % -27 56  % Central channel
cell 253 {ul2} hastelloy -53  20 % -27 56  % Hast. pipe
cell 254 {ul2} fuel -54  53 % -27 56  % outer fuel channel
cell 257 {ul2} hastelloy -57 54 % -27 56 % outer pipe
cell 259 {ul2} void 11

% Universe {ulp}: LOWER PLENUM
cell 250 {ulp} fuel -11
cell 260 {ulp} void  11

% Universe {ul3}: PENENETRATION TO INLET PLENUM
cell 261 {ul3} hastelloy -11 54
cell 262 {ul3} fuel -20          % Central channel
cell 263 {ul3} hastelloy -53  20 % Inner pipe
cell 264 {ul3} fuel -54  53      % outer fuel channel
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
cell 71 {ulc} he -20               % Central channel
cell 72 {ulc} hastelloy -53  20    % Hast. pipe
cell 73 {ulc} he -54  53           % outer fuel channel
cell 74 {ulc} graphite -55 54      % graphite hex
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
cell 107 0 fill   {latbot3} -{surffuel}	  	61 -56
cell 108 0 fill   {latplen} -{surffuel}	  	62 -61
cell 109 0 fill   {latbot4} -{surffuel}	  	63 -62
cell 110 0 fill   {latplen} -{surffuel}	  	64 -63
cell 111 0 hastelloy       -{surfgref}		60 -64
cell 112 0 hastelloy {surffuel} -{surfgref} 64 -61
cell 113 0 hastelloy {surfgref} -{surfhast} 60 -81
cell 120 0 blanket  	 	    -{surfgref} 42 -81 %-80  -{surfcore}
%cell 121 0 blanket 	 	    -{surfcore} 80 -81  % was top reflector
cell 122 0 hastelloy	 	    -{surfhast} 81 -82



cell 999 0 outside {surfhast} 60 -82
cell 998 0 outside -60
cell 997 0 outside  82
	'''



# New cell definitions
	'''cells = \'''
% Universe {uf}: FUEL CELL
cell 10 {uf} graphite -101  103 105
cell 13 {uf} graphite  103 -104
cell 11 {uf} blanket  -102  101 
cell 12 {uf} void      102 
cell  7 {uf} fuel     -103 
cell  1 {uf} fuel      104 -105 

% Universe {ub}: BLANKET CELL
cell 15 {ub} graphite -301 
cell 16 {ub} blanket   301 

% Universe {uh}: HASTELLOY HEX
cell 18 {uh} hastelloy -102
cell 19 {uh} void       102

% Universe {uc}: CONTROL ROD 
cell 30 {uc} graphite -201 203 204 205 206 207 208 209
cell 41 {uc} blanket  -202 201 
cell 401 {uc} blanket -203
cell 42 {uc} void      202
%cell 38 {uc} void    -298
%cell 39 {uc} void     299
cell  34 {uc} he      -204 
cell  35 {uc} he      -205 
cell  36 {uc} he      -206
cell  37 {uc} he      -207 
cell  38 {uc} he      -208 
cell  39 {uc} he      -209 

% Universe {uup}: UPPER CHANNEL !! holding shaft not included !!
cell 43 {uup} fuel     -401  496 -497       % fuel cap
cell 44 {uup} graphite  401 -402  496 -497	% graphite,level1
cell 45 {uup} graphite -402  497 -498       % graphite cap,level2
cell 46 {uup} blanket  -403  402  496  498  % slit all the way through
cell 47 {uup} void      498
cell 48 {uup} void     -496
cell 49 {uup} void      403  496 -498

% Universe {ul1}: LOWER CHANNEL 1
cell 51 {ul1} blanket   -701 702
cell 52 {ul1} fuel      -703     % Central channel
cell 53 {ul1} hastelloy -704 703 % Hast. pipe
cell 54 {ul1} fuel      -705 704 % outer fuel channel
cell 55 {ul1} graphite  -702 705 % graphite hex
cell 59 {ul1} void       701

% Universe {ul2}: LOWER CHANNEL 2
cell 251 {ul2} blanket   -801 805
cell 252 {ul2} fuel      -802     % Central channel
cell 253 {ul2} hastelloy -803 802 % Hast. pipe
cell 254 {ul2} fuel      -804 803 % outer fuel channel
cell 257 {ul2} hastelloy -805 804 % outer pipe
cell 259 {ul2} void       801

% Universe {ulp}: LOWER PLENUM
cell 250 {ulp} fuel -901
cell 260 {ulp} void  901

% Universe {ul3}: PENENETRATION TO INLET PLENUM
cell 261 {ul3} hastelloy -1001 1002
cell 262 {ul3} fuel      -1003        % Central channel
cell 263 {ul3} hastelloy -1004 1003 % Inner pipe
cell 264 {ul3} fuel      -1004 1002  % outer fuel channel
cell 265 {ul3} void       1001

% Universe {ul4}: PENETRATION TO OUTLET PLENUM
cell 266 {ul4} hastelloy -1101 1102
cell 267 {ul4} fuel      -1102
cell 268 {ul4} void       1101
	'''

	cells = cells.format(**locals())
	return cells

if __name__ == '__main__':
	print("This is a module to write cells for the MSR core.")
	input("Press Ctrl+C to quit, or enter else to test it. ")
	print (write_cells())
