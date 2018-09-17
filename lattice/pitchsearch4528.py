#!/usr/bin/python3
#
# LFTR fitness search on lattce pitch.
#
# Ondrej Chvala, ochvala@utk.edu

import numpy as np
import matplotlib.pyplot as plt

from latgen4528 import LatGen, validate_geometry
from slitcritsearch4528 import SlitCritSearch

my_debug:bool = True

class PitchSearch(object):
    'Find k(CR) relation for different pitches, keeping the rest of the ORNL-4528 lattice the scaled or same'
    def __init__(self, stretch_geometry:int=1, pitchsearch_path:str='./.runs'):
        'Constructor'
        # Starting geometry
        self.l:float    = 7.09  # Hex lattice apothem size / half pitch [cm]
        self.r1:float   = 1.94  # Radius of the inner fuel channel log [cm]
        self.r2:float   = 2.86  # Radius of the first graphite ring [cm]
        self.r3:float   = 3.45  # Radius of the outer fuel channel [cm]
        self.l4:float   = 6.83  # Apothem of the graphite hex [cm]

        self.min:float  =   4.0 # Minimum half-pitch [cm]
        self.max:float  =  14.0 # Maximum half-pitch [cm]
        self.N:int      =    11 # Number of steps in half-pitch
        self.critsearches = []  # List of slit-width crticality searches
        self.stretch_geometry:int = stretch_geometry # Shoudl we scale the rest of the lattice?
        self.pitchsearch_path:str = pitchsearch_path # TODO: THIS IS IGNORED

        self.steps = np.array([]) # Array with half-pitch steps [cm]
        self.crs   = np.array([]) # Array of fixed CR
        self.ks    = np.array([]) # Array k(l) at fixed CR, one row for each CR

    def prepare_searches(self):
        'Setup critsearches and lattices'
        self.steps = np.linspace(self.min, self.max, self.N)
        for my_l in self.steps:
            s   = SlitCritSearch()
            s.l = my_l          # Define geometry for slit search
            if self.stretch_geometry:
                s.r1 = self.r1 * my_l/self.l
                s.r2 = self.r2 * my_l/self.l
                s.r3 = self.r3 * my_l/self.l
                s.l4 = self.l4 * my_l/self.l
            else:
                s.r1 = self.r1
                s.r2 = self.r2
                s.r3 = self.r3
                s.l4 = self.l4
            if not validate_geometry(s.l, s.r1, s.r2, s.r3, s.l4): # Make sure the lattice dimensions are consistent
                if my_debug:
                    print("[DBG pitch] geometry problem, skipping ", s.l, s.r1, s.r2, s.r3, s.l4)
                l_to_remove  = np.array([my_l])   # Remove my_r1 from the steps
                self.steps = np.setdiff1d(self.steps, l_to_remove)
                continue
            s.set_path_from_geometry()  # TODO: so far we ignore pitchsearch_path
            s.build_lattices()
            self.critsearches.append(s) # Add slit crit search to the list
        self.ks = np.empty((0,self.steps.size), dtype=np.float64)   # Empty array with the correct shape

    def run_searches(self, sleep_sec:int = -1):
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
            s.plot("fit_l_%08.5f.png" % s.l)

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
            self.ks = np.vstack([self.ks, my_ks])
            if my_debug:
                print("[DEBUG] crs:", self.crs)
                print("[DEBUG] ks:",  self.ks)
        return self.ks[np.where(self.crs == cr)].flatten()

    def plot_ks_at_cr(self, cr=0.97, plot_file='k_vs_l.pdf'):
        'Plot k|cr vs l'
        if self.steps.size == 0:
            print("[ERROR] Search not ready, run prepare_searches()")
            return False
        # Plot data points
        plt.plot(self.steps, self.get_ks_at_cr(cr), color="green",  marker="+", ls="none")
        plt.title("k(l) at CR=%5.2f" % cr)
        plt.xlabel("Hex apothem [cm]")
        plt.ylabel(r"k$_{inf}$")
        plt.grid(True)
#        plt.legend(loc="best", fontsize="medium", title="Parameters of the polynomial fits")
        if plot_file == None:
            plt.show()
        else:
            plt.savefig(plot_file, bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    print("This is a module to run pitch search on ORNL-4528 lattice.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    stretch_geometry = 1            # 1, 0
    p = PitchSearch(stretch_geometry) # Call class constructor
    p.min =  4.0                    # Change default search parameters
    p.max = 25.0                    # Change default search parameters
    p.N   = 43                      # Change default search parameters
    p.prepare_searches()            # Build all corresponsing slit searches
    print(p.steps.size, p.steps)    # Show slit search setting
    p.run_searches()                # Execute all slit searches
    print(p.get_ks_at_cr())         # Print results
    print(p.get_ks_at_cr(1.0))
    p.plot_ks_at_cr(0.97,"./plots/l_%d_k_at_cr097pct.png" % stretch_geometry )
    p.plot_ks_at_cr(0.95,"./plots/l_%d_k_at_cr095pct.png" % stretch_geometry )
    p.plot_ks_at_cr(0.98,"./plots/l_%d_k_at_cr098pct.png" % stretch_geometry )
    p.plot_ks_at_cr(1.00,"./plots/l_%d_k_at_cr100pct.png" % stretch_geometry )

