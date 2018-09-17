#!/usr/bin/python3
#
# Class that generates a ORNL-4528 core.
#
# Ondrej Chvala, ochvala@utk.edu


import math
import numpy as np
import sys
import os

my_debug:bool = False
min_graphite_width:float = 0.05     # Minimum graphite width
min_fuelsalt_width:float = 0.04     # Minimum fuel channel width
min_blanket_width:float  = 0.01     # Minimum blanket slit width
std_bf:float             = 0.928002 # Nominal blanket fraction

def validate_geometry(l, r1, r2, r3, l4) ->bool:
    'Geometric check on the lattice parameters'
    if r1 < r2 < r3 < l4 < l:
        pass
    else:
        print("ERR: Lattice dimensions are inconsistent!")
        return False
    if abs((r1**2 - r3**2+r2**2)/(r1**2 - r3**2+r2**2)) > 1e3:
        print("ERR: Fuel salt rings have different area!")
        return False
    if r1<min_fuelsalt_width or (r3-r2)<min_fuelsalt_width:
        print("ERR: Fuel channel too small!")
        return False
    if (r2-r1)<min_graphite_width:
        print("ERR: Graphite ring too small!")
        return False
    if (l4-r3)<min_graphite_width:
        print("ERR: Outer graphite thickness too small!")
        return False
    if (l-l4)<min_blanket_width:
        print("ERR: Blanket slit too small!")
        return False
    return True

class CoreGen(object):
    'ORNL-4528 core writer class'
    def __init__(self,  l:float=6.045, r1:float=2.550, r2:float=3.758, r3:float=4.534, l4:float=5.669, \
        tempC:float=700, rgr_scale:float=0.9, rfuel:float=160.0, rcore:float=220.0, zcore:float=400.0, \
        refl_ht:float = 100.0, deckname = 'ORNL-4528 deck'):

        # Do basic geometry checks first
        if not validate_geometry(l, r1, r2, r3, l4):
            raise ValueError("Lattice geometry failed checks!")
        if rgr_scale > 1.0:
            raise ValueError("Blanket graphite scaling error!")
        if rfuel > rcore:
            raise ValueError("Core geometry error!")

        # Core lattice geometry
        self.l:float        = l         # Hex lattice apothem size / half pitch [cm]
        self.r1:float       = r1        # Radius of the inner fuel channel log [cm]
        self.r2:float       = r2        # Radius of the first graphite ring [cm]
        self.r3:float       = r3        # Radius of the outer fuel channel [cm]
        self.l4:float       = l4        # Apothem of the graphite hex [cm]
        # Core parameters
        self.tempC:float    = tempC     # Salt temperature [C]
        self.rgs:float      = rgr_scale # Radial reflector scaling fraction
        self.rfuel:float    = rfuel     # Radius of fuel elements [cm]
        self.rcore:float    = rcore     # Radius of entire core (with blanket) [cm]
        self.zcore:float    = zcore     # Axial core height [cm]
        self.reflht:float   = refl_ht   # Axial reflector height [cm]
        self.deckname:str   = deckname  # Title of the Serpent deck

        self.deck_path:str  = '.'       # Where to run the  deck
        self.deck_file_name:str = 'ornl4528.inp'    # Serpent input file name
        self.geomplots:bool = False     # Do Serpent geometry plots
        self.meshplots:bool = False     # Do Serpent mesh plots
        self.queue:str      = 'fill'    # NEcluster torque queue
        self.ompcores:int   = 8         # OMP core count
        self.blanket_drained:bool = False # If True, blanket salt is replaced by helium
        self.absorber:str   = 'enrb4c'  # Set absorber material ['natb4c', 'enrb4c', 'boronmetal']

        # Control rods' lattice positions can be set in the tuple below.
        self.control_rods = [] # example: [(0,0),(4,-4),(4,0),(0,4),(-4,0),(0,-4),(-4, 4)]
        # Control rod can be all out 0, cluster in 1, floating center in 2, all in 3
        # CRod universes are encoded as (crod_state*30 + uc)
        self.crod_state   = [] # example: [  0  ,   0  ,  0  ,  0  ,   0  ,  0   ,   0   ]

        self.fuel_cells = int(self.rfuel/self.pitch())
        self.blan_cells = 1

        self.LATS           = range(50,50+7)
        # define the relevant universe numbers
        ub, uf, uc, uup, uuc, ulc = range(1,7)
        ul1, ul2, ul3, ul4 = range(25,29)
        ulp  = 10
        uh   = 11

        # Tuple of all the universe numbers
        self.UNIVERSES = (
        ub,     # blanket cell
        uf,     # fuel cell
        uc,     # control rod
        uup,    # upper channel fuel
        ul1,    # lower channel 1
        ul2,    # lower channel 2
        ul3,    # inlet plenum penetration
        ul4,    # outlet plenum penetration
        ulp,    # lower plenum (pure fuel)
        uuc,    # upper control
        ulc,    # lower control
        uh)     # pure hastelloy hex

        # Height of each plenum: inlet and outlet
        plenum_vol = 37.0*28316.8       # 37 ft^3 to cm^2
        gt = 6.0*2.54                   # Thickness of graphite reflector: cm
        ht = 2.0*2.54                   # Thickness of hastelloy, placeholder
        self.plenum_ht = plenum_vol / (2.0*math.pi*rfuel)**2
        self.rgref = self.rcore + gt
        self.rhast = self.rgref + ht
        self.rcore_inner = self.rcore - ht
        self.rcore_outer = self.rcore

        self.k:float        = -1.0  # Infinite lattice multiplication factor
        self.kerr:float     = -1.0  # Sigma k
        self.cr:float       = -1.0  # Infinite lattice converison ratio
        self.crerr:float    = -1.0  # Sigma CR

    def lattice_tuple(self):
        'Returns a tuple that can generate core lattice = LatGen(*(core.lattice_tuple())) object'
        return self.l, self.r1, self.r2, self.r3, self.l4

    def pitch(self) ->float:
        'Lattice pitch'
        return 2.0*self.l  # Hex pitch [cm]

    def hexarea(self) ->float:
        'Area of the lattice [cm2]'
        return 2.0 * math.sqrt(3.0) * self.l**2

    def R(self) ->float:
        'Outside radius / side [cm]'
        return 2.0 * self.l / math.sqrt(3.0)

    def slit(self) ->float:
        'Blanket slit width [cm]'
        return self.l-self.l4

    def sf(self) ->float:
        'Fuel salt fraction inside the lattice'
        return math.pi*(self.r1**2 + self.r3**2-self.r2**2) / self.hexarea()

    def bf(self) ->float:
        'Blanket salt fraction inside the lattice'
        return self.l4**2 / self.l**2

    def rel_bf(self) ->float:
        'Relative Blanket salt fraction '
        return self.bf() / std_bf

    def nominal_l4(self) ->float:
        'Return nominal l4 using standard blanket fraction'
        return math.sqrt(std_bf * self.l**2 )

    def save_deck(self):
        'Saves Serpent deck into an input file'
        s2_deck = self.serpent_deck()
        if not os.path.exists(self.deck_path):
            try:                    # Create directory for the job
                os.makedirs(self.deck_path)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise Exception("Unable to create directory: ", self.deck_path, str(e))
                    return False
        try:
            fh = open(self.deck_path + '/' + self.deck_file_name, 'w')
            fh.write(s2_deck)
            fh.close()
        except IOError as e:
            print("[ERROR] Unable to write to file: ", self.deck_file_name)
            print(e)

    def run_deck(self, q=None, ppn=None, qsub_script_name='serp_core.sh'):
        'Saves and runs deck on necluster'
        if q is None:
            q = self.queue
        if ppn is None:
            ppn = self.ompcores
        self.save_deck()
        qsub_script = '''\
#!/bin/bash

#PBS -V
#PBS -q {q}
#PBS -l nodes=1:ppn={ppn}
hostname
'''.format(**locals())
        if q is not 'local':
            qsub_script += '''
module load mpi
module load serpent
cd ${PBS_O_WORKDIR}

'''
        qsub_script += '''rm -f done.out
sss2 -omp {ppn} {self.deck_file_name}'''.format(**locals())
        qsub_script += '''
awk 'BEGIN{ORS="\\t"} /ANA_KEFF/ || /CONVERSION/ {print $7" "$8;}' '''
        qsub_script += '''{self.deck_file_name}_res.m > done.out
grep LFTR {self.deck_file_name} | sed -e s/[a-Z,\\"]//g  >> done.out
'''.format(**locals())

        try:
            fh = open(self.deck_path +"/"+ qsub_script_name, 'w')
            fh.write(qsub_script)
            fh.close()
        except IOError as e:
            print("[ERROR] Unable to write to file: ", qsub_script_name)
            print(e)
            raise
            return False
        if q is 'local':    # Run the deck locally
            os.system('cd ' + self.deck_path + '; sh ./' + qsub_script_name)
        else:               # Submit the job on the cluster
            os.system('cd ' + self.deck_path + ';  qsub ' + qsub_script_name)

    def get_calculated_values(self) -> bool:
        'Fill k and cr for lattice if calculated'
        if os.path.exists(self.deck_path+'/done.out') and os.path.getsize(self.deck_path+'/done.out') > 0:
            pass
        else:                   # Calculation not done yet
            return False
        res_file_name = self.deck_path + '/' + "{self.deck_file_name}_res.m".format(**locals())
        fh = open(res_file_name, 'r')
        for line in fh:
            if "ANA_KEFF" in line:
                self.k     = float(line.split()[6])
                self.kerr  = float(line.split()[7])
            if "CONVERSION_RATIO" in line:
                self.cr    = float(line.split()[6])
                self.crerr = float(line.split()[7])
        if my_debug:
            print("[DEBUG Lat] ---> k = {self.k}, CR = {self.cr}".format(**locals()))
        return True

    def serpent_deck(self) -> str:
        'Returns core deck'
        output     = '''set title "{self.deckname}" '''

        surface_cards = self.write_surfs()
        output += surface_cards

        cell_cards = self.write_cells(18, 31, 17, 20)
        output += cell_cards

        # Unpack the universe tuple
        ub, uf, uc, uup, \
        ul1, ul2, ul3, ul4, \
        ulp, uuc, ulc, uh   = self.UNIVERSES
        ubb  =  7 # blank blanket universe
        uhsu =  8 # holding shafts upper
        uhsl =  9 # holding shafts lower
        uhp  = 12 # holding plate

        # Create the middle/active core
        lattice_cards =  self.write_lattice(self.LATS[0], ub, uf, uc)
        # Create the upper plenum
        lattice_cards += self.write_lattice(self.LATS[1], ub, uup, uc)
        # Create the lower level -1
        lattice_cards += self.write_lattice(self.LATS[2], ub, ul1, uc)
        # Create the lower level -2
        lattice_cards += self.write_lattice(self.LATS[3], ub, ul2, uc)
        # Create the lower level -3: penetration to inlet plenum
        lattice_cards += self.write_lattice(self.LATS[4], ub, ul3, uc)
        # Create the lower level -4: penetration to the outlet plenum
        lattice_cards += self.write_lattice(self.LATS[5], uh, ul4, uc)
        # Create the lower fuel plena (identical for both)
        lattice_cards += self.write_lattice(43, uh, 1111,ulp)
        lattice_cards += self.write_lattice(self.LATS[6], uh, ulp, ulp)
        # Create the upper holding shafts
        lattice_cards += self.write_lattice(40, uhsu, uhsu, uhsu)
        # Create lower holding shafts
        lattice_cards += self.write_lattice(41, uhsl, uhsl, uhsl)
        # Create holding plate - first uhp was 11
        lattice_cards += self.write_lattice(42, uhp, uhp, uhp)
        output += lattice_cards

        mat_cards = self.write_materials()
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
'''
        if self.geomplots:
            plot_cards += '''%plot 1 3000 3000 0  -300 300  -80 560
%plot 2 3000 3000 0  -300 300  -80 560
%plot 3 3000 3000 29 %[250 -100 100 -100 100]
plot 1 6000 6000 0  -300 300  -80 560
plot 2 6000 6000 0  -300 300  -80 560
plot 3 6000 6000 29 %[250 -100 100 -100 100]
'''
        if self.meshplots:
            plot_cards += '''mesh 1 5000 5000
mesh 2 5000 5000
mesh 3 5000 5000
'''
        output += plot_cards

        output = output.format(**locals())

        return output


    def write_lattice(self, nlat = '1', ub = '1', uf = '2', uc = '3'):
        '''Function to write the hexagonal lattice of the given type
    Accepts as input: these are old
        nlat:       number identifying the lattice
        ub, uf, uc: the universe numbers of the cells for the MSR core:
                    blanket, fuel, central (int);
        nf, nb:     the number of fuel/blanket cells in radius (integer);
        p:          the pitch (cm) between cells (float)
    Outputs:
        lattice:    string containing the generated lattice card'''

        # Structure of a SERPENT lattice card
        '''lat <u0> <type> <x0> <y0> <nx> <ny> <p>
        #
        # <u0> is the universe number of the lattice
        # <type> is the lattice type: 2 is our hexagon (hexyc)
        # <x0>, <y0> are the x/y coords of the origin--ours are (0,0)
        # <nx>, <ny> are the number of lattice elements in each direction
        # <p> is the lattice pitch'''

        # Local vars
        rmax    = self.rfuel
        pitch   = self.pitch()
        core    = self.rcore_inner
        nf      = self.fuel_cells
        nb      = self.blan_cells

        my_lat_DEBUG = 2
        if my_lat_DEBUG > 2:
            print (self.control_rods)

        # Lattice calculation
        ''' lattice dependencies:  pitch,  rmax,   core    '''
        # and these go in a x type hex lattice
        n = 2*int(math.ceil(core/pitch))+10
        n += 1 if n%2==0 else 0 # make sure it is an odd pin num to center it

        lattice = '''
% This lattice was generated by a script.
lat {nlat} 2  0.0 0.0  {n} {n}  {pitch} \n'''

        # iterate through lattice positions, only place pins where they lie fully in-mod
        # lattice starts at bottom left
        x0 = 0; y0 = 0
        x0 -= (n+0.5*n)*pitch/2.0 - 0.75*pitch
        y0 -= n*math.sqrt(3.0)/4.0*pitch - math.sqrt(3.0)/4.0*pitch
        for i in range(n):                              # y positions
            y  = y0 + i*math.sqrt(3.0)/2.0*pitch        # y pos [cm]
            iy = math.ceil(i-n/2)                       # lattice position along y from the center point
            thisline=''                                 # lattice line in input file
            thisline += ' '*i                           # indent to show hex better
            for j in range(n):                          # x positions
                x  = x0 + (j+0.5*i)*pitch               # x pos [cm]
                ix = math.ceil(j-n/2)                   # lattice position along x from the center point
                lat_r = math.sqrt(x**2+y**2)            # radius from lattice center
                if my_lat_DEBUG > 10:
                    print( x, y, lat_r, lat_r < rmax)
                if lat_r < rmax:                        # if we are inside active core
                    # if lattice radius is in the center, write the unique central channel
                    if (ix,iy) in self.control_rods:
                        if uc is self.UNIVERSES[2]:     # actual CRod
                            my_cridx = self.control_rods.index((ix,iy)) # find CR state
                            my_uc = 30*self.crod_state[my_cridx] + uc
                            thisline += str(my_uc) + ' '# special control rod element
                        else:
                            thisline += str(uc) + ' '
                    else:
                        thisline += str(uf) + ' '       # add a fuel channel to geometry
                elif lat_r < core:
                    thisline += str(ub) + ' '           # add a blanket
                else:
                    thisline += str(ub) + ' '           # universe 7 is a blank blanket
            thisline += '\n'
            lattice += thisline

        lattice = lattice.format(**locals())
        return(lattice)


    def write_cells(self, surffuel=30, surfcore=31, surfgref=32, surfhast=29):
        '''Function to write cell cards for Serpent input deck.
    Inputs: these are old
        universes:      Tuple(12) containing the following.
         ub, uf, uc:     Universe number of blanket, fuel, central cells
         uup, ul#:       Universe number of upper/lower plenum assembly
         uuc, ulc:       Universe number of upper/lower control cells
        lattices:       Tuple(7) containing the following:
         latmid:         Lattice number of the middle fuel cells
         lattop:         Lattice number of the top fuel cells
         latbot#:        Lattice number of the bottom fuel cells
        surffuel:       Surface no. of inner core (fuel cells)
        surfcore:       Surface no. of entire core, fuel+blanket
        surfgref:       Surface no. of graphite reflector
    Outputs:
        cells:      String containing cell cards'''

        # Unpack the universe tuple
        ub, uf, uc, uup, \
        ul1, ul2, ul3, ul4, \
        ulp, uuc, ulc, uh   = self.UNIVERSES
        # Unpack the lattice tuple
        latmid, lattop, latbot1, latbot2, latbot3, \
        latbot4, latplen    = self.LATS

        # Control rod states
        uc0  = 30*0 + uc    # all out
        uc1  = 30*1 + uc    # cluster in
        uc2  = 30*2 + uc    # floating center in
        uc3  = 30*3 + uc    # all in

        # Blanket is helium?
        b_mat = 'blanket'
        if self.blanket_drained:
            b_mat = 'he'

        # Absorber
        absorber = self.absorber

        # cell definitions
        cells = '''
%------define cells--------------------
% Universe {ub}: BLANKET CELL U1
cell 15 {ub} graphite   -101
cell 16 {ub} {b_mat}     101

% Universe {uf}: FUEL CELL U2
cell 10 {uf} graphite -201 203 205
cell 13 {uf} graphite 203 -204
cell 11 {uf} {b_mat}  -202 201
cell 12 {uf} void  202
cell  7 {uf} fuel -203
cell  1 {uf} fuel 204 -205

% Universe {uc0}: CONTROL ROD U3 - all out
% helium in rods
cell  30 {uc0} graphite -302 303 304 305 306 307 308 309
cell  41 {uc0} {b_mat}  -303
cell  42 {uc0} void      302
cell  34 {uc0} he -304
cell  35 {uc0} he -305
cell  36 {uc0} he -306
cell  37 {uc0} he -307
cell  38 {uc0} he -308
cell  39 {uc0} he -309

% Universe {uc1}: CONTROL ROD U33 - cluster in
cell 330 {uc1} graphite   -302 303 304 305 306 307 308 309
cell 341 {uc1} {b_mat}    -303
cell 342 {uc1} void        302
cell 334 {uc1} {absorber} -304
cell 335 {uc1} {absorber} -305
cell 336 {uc1} {absorber} -306
cell 337 {uc1} {absorber} -307
cell 338 {uc1} {absorber} -308
cell 339 {uc1} {absorber} -309

% Universe {uc2}: CONTROL ROD U63 - central element in
cell 630 {uc2} graphite -302 303 304 305 306 307 308 309
cell 641 {uc2} {absorber} -303
cell 642 {uc2} void        302
cell 634 {uc2} he -304
cell 635 {uc2} he -305
cell 636 {uc2} he -306
cell 637 {uc2} he -307
cell 638 {uc2} he -308
cell 639 {uc2} he -309

% Universe {uc3}: CONTROL ROD U93 - all absorbers in
cell 930 {uc3} graphite -302 303 304 305 306 307 308 309
cell 941 {uc3} {absorber} -303
cell 942 {uc3} void        302
cell 934 {uc3} {absorber} -304
cell 935 {uc3} {absorber} -305
cell 936 {uc3} {absorber} -306
cell 937 {uc3} {absorber} -307
cell 938 {uc3} {absorber} -308
cell 939 {uc3} {absorber} -309


% Universe {uup}: UPPER CHANNEL U4
cell 43 {uup} fuel     -403  404 -405       % fuel cap
cell 44 {uup} graphite  403 -401  404 -405  % graphite,level1
cell 45 {uup} graphite -401  405 -406       % graphite cap,level2
cell 46 {uup} {b_mat}  -402  401  404 -406  % slit all the way through
cell 47 {uup} void      406
cell 48 {uup} void     -404
cell 49 {uup} void      402  404 -406

% Universe {ulc}: LOWER CONTROL U6
cell 71 {ulc} graphite -602       % Central channel - was he
cell 72 {ulc} graphite -603  602  % Hast. pipe - was hastelelloy
cell 73 {ulc} graphite -604  603  % outer fuel channel - was he
cell 74 {ulc} graphite -605  604  % graphite hex
cell 75 {ulc} {b_mat}  -601  605  % {b_mat} reflector
cell 79 {ulc} void      601

% Universe 7: Blank Blanket Cell U7
cell 701 7 {b_mat} -701
cell 702 7 {b_mat} 701

% Universe 8: Holding Shafts on top of plate U8
cell 801 8 graphite -801
cell 802 8 {b_mat}  801

% Universe 9: Holding Shafts under plate U9
cell 901 9 graphite -901
cell 902 9 {b_mat} 901

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
cell 51 {ul1} {b_mat}   -2501  2507
cell 52 {ul1} fuel      -2502        % Central channel
cell 53 {ul1} hastelloy -2503  2502  % Hast. pipe
cell 58 {ul1} graphite  -2505  2503
cell 54 {ul1} fuel      -2506  2505  % outer fuel channel
cell 55 {ul1} graphite  -2507  2506  % graphite hex
cell 59 {ul1} void       2501

% Universe {ul2}: LOWER CHANNEL 2  U26
cell 251 {ul2} {b_mat}   -2601  2605
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
cell 111 0 hastelloy -{surfgref}        1 -2
cell 112 0 hastelloy  {surffuel} -{surfgref} 2 -5
cell 113 0 hastelloy  {surfgref} -{surfhast} 1 -15
cell 122 0 hastelloy -{surfhast} 15 -16
cell 123 0 fill   41 -{surfgref} 11 -12 % holding shaft below plate
cell 124 0 fill   42 -{surfgref} 12 -13 % top holding plate
cell 125 0 fill   40 -{surfgref} 13 -14
cell 126 0 {b_mat} -17 14 -15

cell 999 0 outside {surfhast} 1 -16
cell 998 0 outside -1
cell 997 0 outside  16
    '''

        cells = cells.format(**locals())
        return cells


    def write_surfs(self):
        '''Function to write the surfaces for our Serpent model
    Inputs: these are old
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
        pht:    height of each of the lower plena
        zcore:  height of the core
        zrefl:  height of the axial reflector
        gr_exp: graphite expansion coefficient m/m K
    Output:
        surfaces:   string containing the surface cards for the MSiBR'''

        # Local vars
        gr_exp      = 4.14*10**(-6)
        dgr         = 15
        center_cr   = 2.6
        # Calculate the actual dimensions based on our parameters
        # circumradius: distance from center to corner
        d   = (self.l - self.slit()) * 2.0/math.sqrt(3.0)
        # radius (outer): auxiliary fuel channel radius
        ro = 1.8
        # c: a constant that determines how far along the circumradius the channels appear
        c = (self.r1 + d*math.sqrt(3.0)/2.0) / ( d*(1.0 + math.sqrt(3.0)/2.0) )
        # X and Y coordinates of ro
        rox = c*d*math.sqrt(3.0)/2.0
        roy = c*d*0.5

        # thast = hastelloy thickness (1/8 in)
        thast = 1.0/8 * 2.54                            # hastelloy thickness
        rdiff = (self.r3 - self.r2)
        # Establish a few additional dimensions
        ry = c*d            # y coord of vertical channel

        # hexf: top channel hexagon
        hexf = math.sqrt(3.0)/2.0*(d*c) + rdiff
        zhexf1 = self.zcore + 2.0*rdiff      # top of fuel channel (salt)
        # Make the "roof" as thick as the "wall" of the graphite hexagon
        zhexf2 = zhexf1 + (self.l - hexf) # top of fuel channel (graphite)
        # Top of the core: Axial reflector
        zrefltop = zhexf2 + self.reflht
        zgreftop = zrefltop + (self.rgref - self.rcore_inner)
        zhasttop = zgreftop + (self.rhast - self.rgref)

        # Bottom (floor) with the transition later
        # Graphite thickness between central channel and aux channel
        gt1 = d*c - rdiff - self.r1
        # Next layer down concentric circles with half the graphite thickness
        rh = self.r1 + thast                        # radius of hastelloy cyl
        r3_bot = math.sqrt(self.r1**2 + rh**2)      # outer fuel radius in lower plenum
        rg = (gt1/2.0 + r3_bot)                     # radius of graphite hex
        # Cut off the bottom 10 inches below that
        ztrans1 = 0 - 10*2.54
        # Next layer down, the concentric circles w/ hastelloy
        rh2 = r3_bot + thast        # Outer hastelloy pipe
        # Cut off this layer 3 inches down
        ztrans2 = ztrans1 - 3.0*2.54

        # Then, the lower plena
        zitop = ztrans2 - self.plenum_ht    # z of the top of the inlet plenum
        zibot = zitop - self.plenum_ht      # z of the bottom of the inlet plenum
        zotop = zibot - self.plenum_ht      # z of top of outlet plenum
        zobot = zotop - self.plenum_ht      # z of bottom of outlet plenum
        # Then, the very bottom of the entire core
        zbot = zobot - (self.rhast - self.rgref)

        # Top axial reflector - HARDCODED
        rr          = self.rgs*self.l   # reflector graphite log size
        hsl         = 3.97              # holding shaft lower radius
        hsu         = hsl*0.8           # holding shaft upper radius
        axial_top   = zhexf2+30.48
        zplate      = axial_top + 12.0
        zshaft      = zplate + 5.0
        ro          = 1.8               # radius (outer): auxiliary fuel channel radius

        # lower plenum cone calculations: naming scheme goes outward from center
        conez1 = ((0-ztrans1)/(self.r2-rh))*self.r2
        conez2 = ((0-ztrans1)/(self.r3-r3_bot))*self.r3
        conez3 = ((0-ztrans1)/(self.l4-rh2))*self.l

        # Thermal expansion: TODO
        l:float     = self.l
        r1:float    = self.r1
        r2:float    = self.r2
        r3:float    = self.r3
        l4:float    = self.l4
        zcore:float = self.zcore
        rgref:float = self.rgref
        rfuel:float = self.rfuel
        rhast:float = self.rhast

        # surface definitions
        slit:float     = l - l4
        halfslit:float = slit / 2.0

        # surface definitions
        surfaces = '''
%------ main universe ------
%% vertical core plane divisions

surf 16 pz   {zhasttop}
surf 15 pz   {zgreftop}
surf 14 pz   {zshaft}           % top of holding shafts
surf 13 pz   {zplate}           % top of holding plate
surf 12 pz   {axial_top}
surf 11 pz   {zhexf2}
surf 10 pz   {zhexf2}
surf 9 pz    {zcore}            % HEIGHT OF CORE
surf 8 pz    0                  % BOTTOM OF CORE
surf 7 pz    {ztrans1}
surf 6 pz    {ztrans2}
surf 5 pz    {zitop}
surf 4 pz    {zibot}
surf 3 pz    {zotop}
surf 2 pz    {zobot}
surf 1 pz    {zbot}

%% radial bounds
surf 17 cyl   0   0   {rgref}   % blanket above core
surf 18 cyl   0   0   {rfuel}
surf 20 cyl   0   0   {rhast}


%------ blanket cell Universe 1 ------
surf 101 hexxc 0   0   {rr}     % HEX FOR REFLECTOR GRAPHITE


%------ graphite hexagon and fuel channel cells Universe 2 ------
surf 201 hexxc 0   0   {l4}     % HEX FOR GRAPHITE
surf 202 hexxc 0   0   {l}      % HEX FOR SLIT
surf 203 cyl   0   0   {r1}     % CENTER HOLE
surf 204 cyl   0   0   {r2}     % INTERMEDIATE GRAPHITE RING
surf 205 cyl   0   0   {r3}     % OUTER FUEL RING


%------ graphite hexagon and control rod channels Universe 3 ------
surf 301 hexxc 0   0   {l4}         % HEX FOR GRAPHITE
surf 302 hexxc 0   0   {l}          % HEX FOR SLIT
surf 303 cyl   0   0   {center_cr}  % CENTER HOLE
surf 304 cyl   0   {ry}   {ro}      % OUTER HOLES x 6 CONTROL RODS
surf 305 cyl   0  -{ry}   {ro}      %           ||
surf 306 cyl   {rox}  {roy}  {ro}   %           ||
surf 307 cyl  -{rox}  {roy}  {ro}   %          _||_
surf 308 cyl  -{rox} -{roy}  {ro}   %          \  /
surf 309 cyl   {rox} -{roy}  {ro}   % __________\/____________


%------ upper channel Universe 4 ------
surf 401 hexxc 0   0   {l4}         % HEX FOR GRAPHITE
surf 402 hexxc 0   0   {l}          % HEX FOR SLIT
surf 403 cyl   0   0   {r3}         % OUTER FUEL RING
surf 404 pz    {zcore}              % HEIGHT OF CORE
surf 405 pz    {zhexf1}
surf 406 pz    {zhexf2}


%------ upper control Universe 5 ------
surf 501 hexxc 0   0   {l4}         % HEX FOR GRAPHITE
surf 502 hexxc 0   0   {l}          % HEX FOR SLIT
surf 503 hexxc 0   0   {hexf}       % HEX FOR FUEL
surf 504 pz    {zcore}              % HEIGHT OF CORE
surf 505 pz    {zhexf1}
surf 506 pz    {zhexf2}


%------ lower control Universe 6 ------
surf 601 hexxc 0   0   {l}          % HEX FOR SLIT
surf 602 cyl   0   0   {r1}         % CENTER HOLE
surf 603 cyl   0   0   {rh}         % hastelloy tube radius
surf 604 cyl   0   0   {r3_bot}     % outer fuel radius for lower plenum
surf 605 hexxc 0   0   {rg}         % Hex for fuel transition


%------ blank blanket cell Universe 7 ------
surf 701 hexxc 0   0   {l4}


%------ holding shafts on top of plate Universe 8 ------
surf 801 cyl   0   0   {hsu}        % holding shaft


%------ holding shafts under plate Universe 9 ------
surf 901 cyl   0   0   {hsl}


%------ hastelloy hex Universe 11 ------
surf 1101 hexxc 0   0   {l}         % HEX FOR SLIT


%------ lower plenum bottom Universe 10 ------
surf 1001 hexxc 0   0   {l}         % HEX FOR SLIT


%------ holding plate Universe 12 ------
surf 1201 cylz  0   0   {hsu}


%------ lower channel 1 Universe 25 ------
surf 2501 hexxc 0   0   {l}         % HEX FOR SLIT
surf 2502 cyl   0   0   {r1}        % CENTER HOLE
surf 2503 cyl   0   0   {rh}        % HASTELLOY PIPE RADIUS
surf 2504 cyl   0   0   {r2}        % OUTER FUEL CYL RADIUS
surf 2505 cone  0   0 0 {r2} -{conez1}
surf 2506 cone  0   0 0 {r3} -{conez2}
surf 2507 cone  0   0 0 {l}  -{conez3}


%------ lower channel 2 Universe 26------
surf 2601 hexxc 0   0   {l}         % HEX FOR SLIT
surf 2602 cyl   0   0   {r1}        % CENTER HOLE
surf 2603 cyl   0   0   {rh}        % HASTELLOY PIPE RADIUS
surf 2604 cyl   0   0   {r3_bot}    % outer fuel radius for lower plenum
surf 2605 cyl   0   0   {rh2}       % OUTER HASTELLOY PIPE


%------ penetration to inlet plenum Universe 27 ------
surf 2701 hexxc 0   0   {l}         % HEX FOR SLIT
surf 2702 cyl   0   0   {r1}        % CENTER HOLE
surf 2703 cyl   0   0   {r3_bot}    % outer fuel radius for lower plenum
surf 2704 cyl   0   0   {rh}        % hastelloy tube radius


%------ penetration to outlet plenum Universe 28 ------
surf 2801 hexxc 0   0   {l}         % HEX FOR SLIT
surf 2802 cyl   0   0   {r1}        % CENTER HOLE


%------ lower plenum top Universe 1111 ------
surf 111101 hexxc 0   0   {l}       % HEX FOR SLIT
surf 111102 cyl   0   0   {r1}      % CENTER HOLE
surf 111103 cyl   0   0   {rh}      % hastelloy tube radius

    '''
        surfaces = surfaces.format(**locals())

        return surfaces


    def write_materials(self):
        '''Function to write material cards for Serpent input deck.
Inputs: these are old
    temp: core temperature
    lib:    String containing the neutron cross section library to use.
    scat_lib : thermal scattering library
Outputs:
    mats:    String containing the material cards'''
        tempK = self.tempC + 273.15 # convert C to K to more easily use libraries
        # set neutron cross section library
        if tempK >= 300 and tempK < 600: lib = '03c'
        if tempK >= 600 and tempK < 900: lib = '06c'
        if tempK >= 900 and tempK < 1200: lib = '09c'
        if tempK >= 1200 and tempK < 1500: lib = '12c'

        # set thermal scattering libraries
        if tempK == 400: scat_lib = 'gre7.14t'
        if tempK > 400 and tempK < 600: scat_lib = 'gre7.14t gre7.16t'
        if tempK == 600: scat_lib = 'gre7.16t'
        if tempK > 600 and tempK < 800: scat_lib = 'gre7.16t gre7.18t'
        if tempK == 800: scat_lib = 'gre7.18t'
        if tempK > 800 and tempK < 1000: scat_lib = 'gre7.18t gre7.20t'
        if tempK == 1000: scat_lib = 'gre7.20t'
        if tempK > 1000 and tempK < 1200: scat_lib = 'gre7.20t gre7.22t'
        if tempK == 1200: scat_lib = 'gre7.22t'
        if tempK > 1200 and tempK < 1500: scat_lib = 'gre7.22t gre7.24t'

        # change density of major components to match temperature (reference DMSR project) 973K nominal
        fuel_dens = 2.03434 - (tempK - 973)*1*10**(-3)
        blanket_dens = 4.43711 - (tempK - 973)*1*10**(-3)
        gr_dens = 1.82/((1. + (tempK - 973)*4.14*10**(-6))**3)

        mats = '''
%-------material definition--------------
%NOTE: VOLUMES OR MASS OF EACH MAY NEED TO
%BE CALCULATED FOR BURNUP CALCULATIONS

%  FUEL SALT: 68.5-31.3-0.2 LiF-BeF2-UF4 at 900K
%  DENSITY: 2.03434 g/cc
%  MELT TEMP: 450C or 742.15K
%  MATERIAL INFO FROM ONRL-4528 TALBE 3.1.
%  MAY NEED VOLUME OR MASS FOR BURNUP CALCULATIONS
mat fuel -{fuel_dens}  tmp {tempK}
rgb 130 32 144
3006.{lib}   -0.000725     %  Li6
3007.{lib}  -14.495960     %  Li7
4009.{lib}   -8.508748     %  Be9
9019.{lib}  -75.588074     %  F19
92233.{lib}  -1.265844     %  U233
92234.{lib}  -0.140649     %  U234

%  BLANKET SALT: 71-27-2 LiF-ThF4-BeF2 at 900K
%  DENSITY: 4.43711 g/cc
%  MELT TEMP: 560C or 833.15K
%  MATERIAL INFO FROM ONRL-4528 TALBE 3.1.
%  MAY NEED VOLUME OR MASS FOR BURNUP CALCULATIONS
mat blanket -{blanket_dens} tmp {tempK}
rgb 0 157 254
3006.{lib}   -0.000243     %  Li6
3007.{lib}   -4.855845     %  Li7
4009.{lib}   -0.175712     %  Be9
9019.{lib}  -33.892970     %  F19
90232.{lib} -61.052512     %  Th232
91233.{lib}  -0.022718     %  Pa233

%  NUCLEAR GRAPHITE: Natural concentration of carbon
%  DENSITY: 1.82 G/CC
mat graphite -{gr_dens} moder graph 6000 tmp {tempK}
rgb 130 130 130
6000.{lib} 1
%  THERMAL SCATTERING LIBRARY FOR GRAPHITE
therm graph {tempK} {scat_lib}

%  HELIUM: gas due to alpha particles
%  DENSITY: 54.19 E-6 g/cc
mat he -54.19E-6 tmp {tempK}
rgb 255 245 213
2004.{lib} 1

% CONTROL ROD: NATURAL BORON at 900K
% DENSITY: 2.3 g/cc
% MELT TEMP: 2076C or 2349.15K
% 19.9 B10 and 80.1 B11
mat boronmetal -2.3 tmp {tempK}
rgb 75 75 75
5010.{lib} -0.199
5011.{lib} -0.801

% CONTROL ROD: natural B4C
mat natb4c -2.418 tmp {tempK}
rgb 70 70 70
5010.{lib} -0.14419     %   B10
5011.{lib} -0.63843     %   B11
6000.{lib} -0.21738     %   carbon

% CONTROL ROD: enriched B4C
mat enrb4c -2.296 tmp {tempK}
rgb 65 65 65
5010.{lib} -0.68702     %   B10
5011.{lib} -0.08397     %   B11
6000.{lib} -0.22901     %   carbon

%  Hastelloy
mat hastelloy -8.86 tmp {tempK}
rgb 139 69 19
28058.{lib}  -0.472120   %  Ni
28060.{lib}  -0.181860   %  Ni
28061.{lib}  -0.007905   %  Ni
28062.{lib}  -0.025206   %  Ni
28064.{lib}  -0.006419   %  Ni
42100.{lib}  -0.015408   %  Mo
42092.{lib}  -0.023744   %  Mo
42094.{lib}  -0.014800   %  Mo
42095.{lib}  -0.025472   %  Mo
42096.{lib}  -0.026688   %  Mo
42097.{lib}  -0.015280   %  Mo
42098.{lib}  -0.038608   %  Mo
24050.{lib}  -0.003041   %  Cr
24052.{lib}  -0.058652   %  Cr
24053.{lib}  -0.006651   %  Cr
24054.{lib}  -0.001656   %  Cr
26054.{lib}  -0.002923   %  Fe
26056.{lib}  -0.045877   %  Fe
26057.{lib}  -0.001059   %  Fe
26058.{lib}  -0.000141   %  Fe
14028.{lib}  -0.009223   %  Si
14029.{lib}  -0.000468   %  Si
14030.{lib}  -0.000309   %  Si
25055.{lib}  -0.008000   %  Mn
74182.{lib}  -0.001325   %  W
74183.{lib}  -0.000715   %  W
74184.{lib}  -0.001532   %  W
74186.{lib}  -0.001422   %  W
29063.{lib}  -0.002421   %  Cu
29065.{lib}  -0.001079   %  Cu
'''
        mats = mats.format(**locals())

        return mats

if __name__ == '__main__':
    print("This is a module to write ORNL-4528 core.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    testcore = CoreGen()
    testcore.geomplots = True
    # Insert control rods
    icr = 8  # How far, in lattice locations, are the 6 additional control rods
    # Coltrol rod positions: center + 6 positions around the center
    testcore.control_rods=[(0,0), (icr,-icr), (icr,0), (0,icr), (-icr,0), (0,-icr), (-icr,icr)]
    testcore.crod_state  =[ 0, 0, 0, 0, 0, 0, 0 ]   # All control rods out
    #testcore.crod_state  =[ 0, 1, 2, 3, 0, 2, 2 ]   # Test various control rod states
    testcore.save_deck()            # Save LFTR deck to file
    testcore.run_deck("xeon",32)    # Cluster run, 32 cores
#    testcore.run_deck("local",8)    # Local run, 8 cores

