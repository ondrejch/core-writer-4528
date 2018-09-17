#!/usr/bin/python3
#
# ORNL-4528 parameter search along 2 dimensions
#
# Ondrej Chvala, ochvala@utk.edu

import numpy as np
import matplotlib.pyplot as plt
import os

from latgen4528 import LatGen
from slitcritsearch4528 import SlitCritSearch
from pitchsearch4528 import PitchSearch
from r1search4528 import R1Search
from r2search4528 import R2Search

my_debug = True

class Search2d(object):
    'Find k(CR) relation for two searches'
    def __init__(self, search_path:str='./.runs/'):
        'Constructor'
        self.search_path:str = search_path
        # 2D map
        self.map   = []         # 2D map of crit searches
#        self.nmap  = []         # helping map  of integets, for debugging
        self.crs   = np.array([]) # Array of fixed CR
        self.ks    = np.array([]) # Array k(x1,x2) at fixed CR, one matrix for each CR

        # Map strings to class names
        self.search_classes = { "pitch": PitchSearch, "r1": R1Search, "r2": R2Search }
        self.ax1name:str = "r1"     # 1st scan axis - can be "pitch", "r1, or "r2"
        self.ax2name:str = "pitch"  # 2nd scan axis - can be "pitch", "r1, or "r2"
        self.search1 = None
        self.search2 = None

    def create_searches(self, stretch_geometry_ax1:int=1, stretch_geometry_ax2:int=1):
        'Associates search classes - make sure the names are set first'
        self.search1 = self.search_classes[self.ax1name](stretch_geometry_ax1, self.search_path)
        self.search2 = self.search_classes[self.ax2name](stretch_geometry_ax2, self.search_path)

    def prepare_searches(self):
        'Prepares 2D matrix of searches - make sure the parametrs are set first'
        self.search1.prepare_searches()     # Independent 1D search, axis 1
        self.search2.prepare_searches()     # Independent 1D search, axis 2
        # Interpolate the lattice parameters
        i1 = 0      # Running indeces
        for s1 in self.search1.critsearches:
            self.map.append([])                 # Add list for the new row
#            self.nmap.append([])
            i2 = 0  # Running indeces
            for s2 in self.search2.critsearches:
                s   = SlitCritSearch()
                if my_debug:
                    print("[DBG2D] ", i1, i2)
#                    print("[DBG2D] ", self.nmap)
                # Linear interpolation using [x1, x2]
                s.l  = (s1.l  + s2.l ) /2.0
                s.r1 = (s1.r1 + s2.r1) /2.0
                s.r2 = (s1.r2 + s2.r2) /2.0
                s.r3 = (s1.r3 + s2.r3) /2.0
                s.l4 = (s1.l4 + s2.l4) /2.0
                # Setup latices for the slit crit search
                s.set_path_from_geometry()
                s.build_lattices()
#                self.nmap[i1].append(100+100*i1+i2)    # Add search to the map
                self.map[i1].append(s)    # Add search to the map
                i2 += 1
            i1 += 1
        self.ks = np.empty((0, i1, i2), dtype=np.float64)   # Empty array with the correct shape

    def run_searches(self, sleep_sec:int = -1):
        'Executes the searches'
        # Independent 1D searches first
        self.search1.run_searches()
        self.search2.run_searches()
        # Submit the 2D search map
        for row in self.map:
            for s in row:
                s.qsub_lattices()
        for row in self.map:
            for s in row:
                if sleep_sec>=0:
                    s.sleep_sec = sleep_sec
                s.get_calculated_values()
                s.fit()

    def get_ks_at_cr(self, cr:float=0.97):
        'Evaluate searches'
        if len(self.map) == 0:
            print("[ERROR] Results not ready, prepare/run searches first!")
            return False
        if cr in self.crs:  # Check if we already have calculated this
            pass
        else:
            self.crs = np.append(self.crs, cr)  # Add CR to the list
            my_ks = []
            my_i  = 0
            for row in self.map:
                my_ks.append([])                # Add list for the new row
                for s in row:
                    my_ks[my_i].append( s.get_k_at_cr(cr) )
                my_i += 1
            my_ks = np.flipud(np.array(my_ks).transpose()) # Fixes axis orientations
            if len(self.ks) == 0:
                self.ks = [my_ks]
            else:
                self.ks = np.concatenate((self.ks, [my_ks]), axis=0)
            if my_debug:
                print("[DEBUG S2d] crs:", self.crs)
#                print("[DEBUG] ks:",  self.ks)
#                print("[DEBUG] :",  np.where(self.crs == cr))
        return self.ks[np.where(self.crs == cr)[0][0] ]

    def plot_ks_at_cr(self, cr:float=0.97, plot_file:str='k_vs_r1.pdf', do_1D_plots:bool=True):
        '''Plot k|cr vs ax1/ax2
        #       *** MAJOR WARNING ****
        # The axes IGNORE geometries were REJECTED for geometrical consistency reasons!!!
        # !!! Thus the value on the axes are NOT reliable !!!'''
        if len(self.map) == 0:
            print("[ERROR] Results not ready, prepare/run searches first!")
            return False
        if my_debug:
            print("[DEBUG S2d] plotting:", cr, plot_file)
        # Plot data map
        plt.imshow(self.get_ks_at_cr(cr),extent=(self.search1.min, self.search1.max, self.search2.min, self.search2.max))
        plt.title(r"k$_{inf}$ at CR=%5.2f" % cr)
        plt.xlabel("%s [cm]" % self.ax1name)
        plt.ylabel("%s [cm]" % self.ax2name)
        cbar = plt.colorbar()
        cbar.set_label(r"k$_{inf}$")
        if plot_file == None:
            plt.show()
        else:
            plt.savefig(plot_file, bbox_inches='tight')
        plt.close()
        if do_1D_plots:        # 1D only plots
            self.search1.plot_ks_at_cr(cr, os.path.dirname(plot_file)+"/ax1_"+os.path.basename(plot_file))
            self.search2.plot_ks_at_cr(cr, os.path.dirname(plot_file)+"/ax2_"+os.path.basename(plot_file))

    def get_highest_k_at_cr(self, cr:float=0.97):
        'Find latttice with highest k'
        max_k = -1.0
        max_s = None
        for row in self.map:
            for s in row:
                k = s.get_k_at_cr(cr)
                if k > max_k:
                    max_k = k
                    max_s = s
        if my_debug:
            print("[DEBUG S2d max_k] :", max_k)
        my_l4 = max_s.l - max_s.get_slit_at_cr(cr)
        return (max_k, max_s.l, max_s.r1, max_s.r2, max_s.r3, my_l4)


if __name__ == '__main__':
    print("This is a module to run 2 combined parameter searches on ORNL-4528 lattice.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    s2d = Search2d()        # Call constructor
    s2d.ax1name = "r1"      # 1st scan axis - can be "pitch", "r1, or "r3"
    s2d.ax2name = "pitch"   # 2nd scan axis - can be "pitch", "r1, or "r3"
    s2d.create_searches()   # Associate search types with x/y axes
#    s2d.search1.min =  0.25 # Minimum r1 [cm]   - Change default parameters for the searches
#    s2d.search1.max =  9.0  # Maximum r1 [cm]
#    s2d.search1.N   =   36  # Number of steps in r1
    s2d.search1.stretch_geometry = 1
#    s2d.search2.min =  7.0  # Minimum half-pitch [cm]
#    s2d.search2.max = 17.0  # Maximum half-pitch [cm]
#    s2d.search2.N   =   21  # Number of steps in half-pitch
    s2d.prepare_searches()  # Build all corresponsing slit searches
    s2d.run_searches()      # Execute all slit searches
    print(s2d.get_ks_at_cr())       # Print results
    print(s2d.get_ks_at_cr(1.0))
    print("Highest k at cr = 0.93: ", s2d.get_highest_k_at_cr(0.93))
    print("Highest k at cr = 0.95: ", s2d.get_highest_k_at_cr(0.95))
    print("Highest k at cr = 0.97: ", s2d.get_highest_k_at_cr(0.97))
    s2d.plot_ks_at_cr(0.93,"./plots/lr1_k_at_cr097pct.png")
    s2d.plot_ks_at_cr(0.95,"./plots/lr1_k_at_cr095pct.png")
    s2d.plot_ks_at_cr(0.97,"./plots/lr1_k_at_cr097pct.png")
    s2d.plot_ks_at_cr(0.98,"./plots/lr1_k_at_cr098pct.png")
    s2d.plot_ks_at_cr(1.00,"./plots/lr1_k_at_cr100pct.png")


