#!/usr/bin/python3
#
# Class that generates an ORNL-4528 type lattice.
#
# Ondrej Chvala, ochvala@utk.edu

import math
import numpy as np
import sys
import os

my_debug:bool = True
min_graphite_width:float = 0.05     # minimum graphite width
min_fuelsalt_width:float = 0.04     # minimum fuel channel width
min_blanket_width:float  = 0.01     # minimum blanket slit width
std_bf:float             = 0.928002 # nominal blanket fraction

def nominal_l4(l:float) ->float:
    'Return nominal l4 using standard blanket fraction'
    return math.sqrt(std_bf * l**2 )

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

class LatGen(object):
    'ORNL-4518 lattice generator class'
    def __init__(self, l:float=7.09, r1:float=1.94, r2:float=2.86, r3:float=3.45, \
                 l4:float=6.83):
        'Geometry based constructor with sane defaults'
        if my_debug:
             print("[DEBUG lat] Creating: ", l, r1, r2, r3, l4)
        # Do basic geometry checks first
        if not validate_geometry(l, r1, r2, r3, l4):
            raise ValueError("Lattice geometry failed checks!")

        self.l:float        = l     # Hex lattice apothem size / half pitch [cm]
        self.r1:float       = r1    # Radius of the inner fuel channel log [cm]
        self.r2:float       = r2    # Radius of the first graphite ring [cm]
        self.r3:float       = r3    # Radius of the outer fuel channel [cm]
        self.l4:float       = l4    # Apothem of the graphite hex [cm]
        self.lib:str        = '09c' # Default Serpent CE x-section library

        self.deck_path:str  = '.'   # Where to run the lattice deck
        self.deck_file_name:str = 'msbr.inp'    # Serpent input file name
        self.queue:str      = 'fill'    # NEcluster torque queue
        self.ompcores:int   = 8         # OMP core count

        self.k:float        = -1.0  # Infinite lattice multiplication factor
        self.kerr:float     = -1.0  # Sigma k
        self.cr:float       = -1.0  # Infinite lattice converison ratio
        self.crerr:float    = -1.0  # Sigma CR

    def slit(self) ->float:
        'Blanket slit width [cm]'
        return self.l - self.l4

    def hexarea(self) ->float:
        'Area of the lattice [cm2]'
        return 2.0 * math.sqrt(3.0) * self.l**2

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

    def R(self) ->float:
        'Outside radius / side [cm]'
        return 2.0 * self.l / math.sqrt(3.0)

    def save_deck(self):
        'Saves Serpent deck into an input file'
        s2_deck = self.serpent_deck()
        try:
            fh = open(self.deck_path + '/' + self.deck_file_name, 'w')
            fh.write(s2_deck)
            fh.close()
        except IOError as e:
            print("[ERROR] Unable to write to file: ", self.deck_file_name)
            print(e)

    def run_deck(self, q=None, ppn=None, qsub_script_name='serp_lat.sh'):
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
grep ORNL-4528 {self.deck_file_name} | sed -e s/[a-Z,\\"]//g  >> done.out
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
        'Builds and returns the Serpent2 deck file'
        my_sf  = self.sf()          # Salt fraction
        my_bf  = self.bf()          # Blanket fraction
        my_rbf = self.rel_bf()      # Relative blanket fraction
        my_slit= self.slit()
        # Deck header
        output = '''\
set title "ORNL-4528 unit cell, l {self.l}, sf {my_sf}, bf {my_bf}, relBA {my_rbf}, radii {self.r1} {self.r2} {self.r3} "
'''
        output += self.write_surfaces()     # Surfaces
        output += self.write_cells()        # Cells
        output += self.write_materials()    # Materials
        # Data cards
        data_cards = '''
%______________data cards___________________________________________
% Boundary condition
set bc 3

% Neutron population and criticality cycles
set pop 10000 200 40 % 10000 neutrons, 50 cycles, 20 of them inactive

% Data Libraries
set acelib "sss_endfb7u.sssdir"
set declib "sss_endfb7.dec"
set nfylib "sss_endfb7.nfy"

% Analog reaction rate
set arr 2
'''
        output += data_cards
        # Plots
        plot_cards = '''
% Plots
plot 3 1500 1500
%mesh 3 1500 1500
'''
        output += plot_cards
        output = output.format(**locals())
        return output

    def write_surfaces(self) -> str:
        'Returns ORNL-4528 surface cards'
        surfaces = '''
%______________surface definitions__________________________________
surf 1   hexxc  0.0 0.0 {self.l}  % reflective unit cell boundary
surf 2   hexxc  0.0 0.0 {self.l4} % graphite boundary
surf 11  cyl    0.0 0.0 {self.r1} % central fuel channel radius
surf 12  cyl    0.0 0.0 {self.r2} % inner channel inner radius
surf 13  cyl    0.0 0.0 {self.r3} % inner channel outer radius

'''
        surfaces = surfaces.format(**locals())
        return surfaces

    def write_cells(self) -> str:
        'Returns ORNL-4528 cell cards'
        cells = '''
%______________cell definitions_____________________________________
cell 11  0  FLiBeU    -11      % inner fuel channel
cell 31  0  graphite   11 -12  % recursive sleeve
cell 12  0  FLiBeU     12 -13  % outer fuel channel
cell 30  0  graphite   13 -2   % moderator block section
cell 20  0  FLiBeTh     2 -1   % blanket region
cell 99  0  outside     1      % graveyard
'''
        cells = cells.format(**locals())
        return cells

    def write_materials(self) -> str:
        'Returns material cards'
        mats = '''
%______________material definitions_________________________________
% LiF-BeF2-UF4 (68.5-31.3-0.2 mole%)
mat FLiBeU -1.93600 tmp 973.0
rgb 130 32 144
 3006.{self.lib}  -0.000725 % Li-6
 3007.{self.lib} -14.495960 % Li-7
 4009.{self.lib}  -8.508748 % Be-9
 9019.{self.lib} -75.588074 % F-19
92233.{self.lib}  -1.265844 % U-233
92234.{self.lib}  -0.140649 % U-234

% LiF-BeF2-ThF4-PaF4 (71.0-2.0-27.0-0.0 mole%)
mat FLiBeTh -4.39764 tmp 973.
rgb 0 157 254
 3006.{self.lib}  -0.000243 % Li-6
 3007.{self.lib}  -4.855845 % Li-7
 4009.{self.lib}  -0.175712 % Be-9
 9019.{self.lib} -33.892970 % F-19
90232.{self.lib} -61.052512 % Th-232
91233.{self.lib}  -0.022718 % Pa-233

%  NUCLEAR GRAPHITE: Natural concentration of carbon
%  DENSITY: 1.82 g/cm^3
mat graphite -1.82 moder graph 6000 tmp 973
rgb 130 130 130
6000.{self.lib} 1
%  THERMAL SCATTERING LIBRARY FOR GRAPHITE
therm graph gre7.08t
'''
        mats = mats.format(**locals())
        return mats


if __name__ == '__main__':
    print("This is a module to write MSBR-4528 lattice.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    testlat = LatGen()
    testlat.save_deck()
    print("bf: ",testlat.bf(), 1.0/ testlat.bf())
    #testlat.run_deck()

