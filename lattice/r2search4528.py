#!/usr/bin/python3
#
# Fitness search on graohite sleeve radius (r2)
#
# Ondrej Chvala, ochvala@utk.edu

import math
import numpy as np
import matplotlib.pyplot as plt

from latgen4528 import LatGen, nominal_l4, validate_geometry
from slitcritsearch4528 import SlitCritSearch

my_debug = True

class R2Search(object):
    'Find k(CR) relation for different diameters of graphite sleeve'
    def __init__(self, stretch_geometry:int=1, search_path:str='./.runs/'):
        'Constructor'
        # Starting geometry
        self.l:float    = 7.09  # Hex lattice apothem size / half pitch [cm]
        self.r1:float   = 1.94  # Radius of the inner fuel channel log [cm]
        self.r2:float   = 2.86  # Radius of the first graphite ring [cm]
        self.r3:float   = 3.45  # Radius of the outer fuel channel [cm]
        self.l4:float   = 6.83  # Apothem of the graphite hex [cm]

        self.min:float  =  2.0  # Minimum r1 [cm]
        self.max:float  =  9.0  # Maximum r1 [cm]
        self.N:int      =   29  # Number of steps in r1
        self.critsearches = []  # List of slit-width crticality searches
        if not stretch_geometry in [0,1,2]:
            raise ValueError("Parameter error: stretch geometry has to be 0, 1, or 2!")
        self.stretch_geometry = stretch_geometry # Shoudl we scale the rest of the lattice?
                # 0 - keep l, r1, scale r3
                # 1 - keep l, scale rest
                # 2 - scale everything
        self.search_path:str  = search_path # TODO: THIS IS IGNORED

        self.steps = np.array([])   # Array with r1 steps [cm]
        self.crs   = np.array([])   # Array of fixed CR
        self.ks    = np.array([])   # Array k(l) at fixed CR, one row for each CR

    def prepare_searches(self):
        'Setup critsearches and lattices'
        self.steps = np.linspace(self.min, self.max, self.N)
        for my_r2 in self.steps:
            s   = SlitCritSearch()
            s.r2 = my_r2          # Define geometry for slit search
            if self.stretch_geometry == 2:  # Scale all
                s.l  = self.l  * my_r2/self.r2
                s.r1 = self.r1 * my_r2/self.r2
                s.r3 = self.r3 * my_r2/self.r2
                s.l4 = self.l4 * my_r2/self.r2
            if self.stretch_geometry == 1:  # Keep l, scale rest
                s.l  = self.l
                s.r1 = self.r1 * my_r2/self.r2
                s.r3 = self.r3 * my_r2/self.r2
                s.l4 = nominal_l4(s.l)
#                if my_debug:
#                    print("[DBGr1] 1 ", s.l, s.r1, s.r2, s.r3, s.l4)
            if self.stretch_geometry == 0:  # Keep l and r2, scale r3 and l4
                s.l  = self.l
                s.r1 = self.r1
                s.r3 = math.sqrt(s.r1**2 + s.r2**2) # Adjust channel sleeve area
                s.l4 = nominal_l4(s.l)
            if not validate_geometry(s.l, s.r1, s.r2, s.r3, s.l4): # Make sure the lattice dimensions are consistent
                if my_debug:
                    print("[DBGr2] geometry problem, skipping ", s.l, s.r1, s.r2, s.r3, s.l4)
                r2_to_remove  = np.array([my_r2])   # Remove my_r1 from the steps
                self.steps = np.setdiff1d(self.steps, r2_to_remove)
                continue
            s.set_path_from_geometry()  # TODO: so far we ignore pitchsearch_path
            s.build_lattices()
            self.critsearches.append(s) # Add slit crit search to the list
            if my_debug:
                print("[DBGr2] done ---> ", my_r2)
        self.ks = np.empty((0,self.steps.size), dtype=np.float64)   # Empty array with the correct shape

    def run_searches(self, sleep_sec = -1):
        'Submits jobs for all rsearches and gathers results'
        if self.steps.size == 0:
            print("[ERROR] Search not ready, run prepare_searches()")
            return False
        for s in self.critsearches:    # Submit jobs
            s.qsub_lattices()
        for s in self.critsearches:    # Get results
            if sleep_sec>=0:
                s.sleep_sec = sleep_sec
            s.get_calculated_values()
            s.fit()

    def plot_searches(self):
        'Plotting is separate since it takes time'
        if self.steps.size == 0:
            print("[ERROR] Search not ready, run prepare_searches()")
            return False
        for s in self.critsearches:    # Make plots
            s.plot("fit_r2_%d_%08.5f.png" % (self.stretch_geometry, s.r2))

    def get_ks_at_cr(self, cr:float=0.97):
        'Evaluate searches'
        if self.steps.size == 0:
            print("[ERROR] Search not ready, run prepare_searches()")
            return False
        if cr in self.crs:  # Check if we already have calculated this
            pass
        else:
            self.crs = np.append(self.crs, cr)  # Add CR to the list
            my_ks = []
            for s in self.critsearches:
                my_ks.append( s.get_k_at_cr(cr) )
#            print("[#####]", cr, self.ks, my_ks)
#            print("[#####]", self.ks.shape, len(my_ks))
            self.ks = np.vstack([self.ks, my_ks])
            if my_debug:
                print("[DEBUG] crs:", self.crs)
                print("[DEBUG] ks:",  self.ks)
        return self.ks[np.where(self.crs == cr)].flatten()

    def plot_ks_at_cr(self, cr:float=0.97, plot_file:str='k_vs_r2.pdf'):
        'Plot k|cr vs r2'
        if self.steps.size == 0:
            print("[ERROR] Search not ready, run prepare_searches()")
            return False
        # Plot data points
        plt.plot(self.steps, self.get_ks_at_cr(cr), color="purple",  marker="+", ls="none")
        plt.title("k(r2) at CR=%5.2f" % cr)
        plt.xlabel("r2 [cm]")
        plt.ylabel(r"k$_{inf}$")
        plt.grid(True)
        if plot_file == None:
            plt.show()
        else:
            plt.savefig(plot_file, bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    print("This is a module to run r2 parameter search on ORNL-4528 lattice.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    stretch_geometry = 2            # 2, 1, 0
    r = R2Search(stretch_geometry)  # Call constructor
    r.prepare_searches()            # Build all corresponsing slit searches
    print(r.steps.size, r.steps)    # Show slit search setting
    r.run_searches()                # Execute all slit searches
    r.plot_searches()               # Plot results of all slit searches
    print(r.get_ks_at_cr())         # Print results
    print(r.get_ks_at_cr(1.0))
    r.plot_ks_at_cr(0.97,"./plots/r2_%d_k_at_cr097pct.png" % stretch_geometry)
    r.plot_ks_at_cr(0.95,"./plots/r2_%d_k_at_cr095pct.png" % stretch_geometry)
    r.plot_ks_at_cr(0.98,"./plots/r2_%d_k_at_cr098pct.png" % stretch_geometry)
    r.plot_ks_at_cr(1.00,"./plots/r2_%d_k_at_cr100pct.png" % stretch_geometry)
