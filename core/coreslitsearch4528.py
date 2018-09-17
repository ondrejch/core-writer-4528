#!/usr/bin/python3
#
# Class to find relations of reactor physics parameters to
# slit with, given the rest of ORNL-4528 core geometry.
#
# Ondrej Chvala, ochvala@utk.edu

import math
from array import array
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
import os
import errno
import time

from ornl4528core import CoreGen, min_blanket_width

my_debug = True

class CoreSlitSearch(object):
    'Establish k_eff(slit_width) and CR(slit_width)'
    def __init__(self, my_core0:CoreGen=None, my_path:str='./.cruns/'):
        'Constructor with default slit criticality search values'
        self.relslit_min:float = 0.80   # Minimum relative slit width
        self.relslit_max:float = 1.20   # Maximum relative slit width
        self.relslit_N:int     = 11     # Number of steps in slit width
        self.c:CoreGen      = my_core0  # Nominal core
        self.corelist          = []     # List of cores

        self.main_path:str  = my_path   # Main directory to run cases
        self.deck_file_name = 'ornl4528.inp'# Serpent lattice deck name
        self.sleep_sec:int  = 60        # Sleep timer between results read attempts [s]

        self.slits      = array('d')    # Slit sizes
        self.k          = array('d')    # ANALYTICAL_KEFF
        self.kerr       = array('d')    # ANALYTICAL_KEFF Error
        self.cr         = array('d')    # CONVERSION_RATIO
        self.crerr      = array('d')    # CONVERSION_RATIO Error
        self.fit_done   = False         # Fit flag
        self.k_fit_p    = [1.0, -1.0, 0.1]  # k fit parametrs
        self.cr_fit_p   = [1.0,  1.0]       # CR fit parametrs
        self.k_chi2:float  = -0.1        # Reduced \Chi^2 of the k(slit) fit
        self.cr_chi2:float = -0.1        # Reduced \Chi^2 of the CR(slit) fit

    def set_nominal_core(self, hpitch:float = 6.045, r1:float = 2.550, r2:float = 3.758, r3:float = 4.534,
        l4:float=5.669, tempC:float = 700., rgr_scale:float = 0.90,  rfuel:float = 200.0, \
        rcore:float = 225.0, zcore:float = 400, refl_ht:float = 100, deckname:str = 'ORNL-4528 deck'):
        'Sets nominal core parameters'
        self.c = CoreGen(hpitch, r1, r2, r3, l4, tempC, rgr_scale, rfuel, rcore, zcore, refl_ht, deckname)

    def set_path_from_geometry(self):
        'Sets path to directory to run cases based on slit geometry'
        if self.c is None:
            print("[ERROR] Nominal core not assigned, run set_nominal_core()!")
            return False
        self.main_path = self.main_path + "/" + "%08.5f"%self.c.l + "/%08.5f"%self.c.r1 + \
            "/%08.5f"%self.c.r2 + "/%08.5f"%self.c.l4
            # r3 omitted since fixed by flow speed up = down

    def build_cores(self):
        'Builds cores for the slit study case'
        if self.c is None:
            print("[ERROR] Nominal core not assigned, run set_nominal_core()!")
            return False
        rel_slit_array = np.linspace(self.relslit_min, self.relslit_max, self.relslit_N)
        if my_debug:
            print("[DEBUG] Nominal l4 ", self.c.l4 )
        nom_slit:float = self.c.slit()  # nominal slit size
        for rslit in rel_slit_array:
            slit:float = rslit * nom_slit
            if slit < min_blanket_width:
                continue            # skip if slit is too small
            my_l4:float = self.c.l - slit
            if my_debug:
                print("[DEBUG] Building core, rslit {rslit}, slit {slit} cm".format(**locals()))
            my_deckname =  self.c.deckname + ', rslit= ' + str(rslit)     # New core deck name
            my_core = CoreGen(self.c.l, self.c.r1, self.c.r2, self.c.r3, my_l4, \
                       self.c.tempC, self.c.rgs, self.c.rfuel, self.c.rcore, \
                       self.c.zcore, self.c.reflht, my_deckname)
            my_core.deck_path = self.main_path +'/'+ str(my_l4)
            my_core.geomplots = False   # This just wastes time and diskspace, not needed in a search
            my_core.meshplots = False
            my_core.control_rods = self.c.control_rods
            my_core.crod_state   = self.c.crod_state
            self.corelist.append( my_core )

    def qsub_cores(self, rerun=False):
        'Run lattice simulations, rerun if rerun=True'
        if self.c is None:
            print("[ERROR] Nominal core not assigned, run set_nominal_core()!")
            return False

        for my_core in self.corelist:
            try:                    # Create directory for each job
                os.makedirs(my_core.deck_path)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise Exception("Unable to create directory: ", my_core.deck_path, str(e))
                    return False
            if (rerun == False) and my_core.get_calculated_values():
                pass # Avoid recalculating existing data
            else:
                if my_debug:
                    print("[DEBUG] Submitting deck", my_core, " in ", my_core.deck_path)
                my_core.run_deck()

    def get_calculated_values(self, reread=True):
        'Fill K and CR values from Serpent results file'
        cores_done:int = 0
        if reread:
            for my_core in self.corelist:        # Flag all lattices for re-read
                my_core.k = -0.1
        while cores_done < len(self.corelist):
            cores_done = 0
            for my_core in self.corelist:
                if(my_core.k < 0.0):
                    if my_core.get_calculated_values():
                        cores_done += 1
                        if my_debug:
                            print("[DEBUG] Got results for ", my_core.deck_path+'/'+my_core.deck_file_name)
                else: # Core already read in
                    cores_done += 1
            if cores_done < len(self.corelist):
                if my_debug:
                    print("[DEBUG] ", cores_done, " done, sleeping ...")
                time.sleep(self.sleep_sec)  # Wait a minute for Serpent ...
        for my_core in self.corelist:    # Fill arrays for fit
            self.slits.append(my_core.slit())
            self.k.append(my_core.k)
            self.kerr.append(my_core.kerr)
            self.cr.append(my_core.cr)
            self.crerr.append(my_core.crerr)

    def fit_f_k(self, x:float, p0:float, p1:float, p2:float) -> float:
        'Fit function for KEFF fit'
        return p2*x**2 + p1*x + p0

    def fit_f_cr(self, x:float, p0:float, p1:float) -> float:
        'Fit function for CR fit'
        return p1*x + p0

    def eval_fit_k(self, x:float):
        'Returns the k-fit value'
        if self.fit_done:
            return self.fit_f_k(x, self.k_fit_p[0], self.k_fit_p[1],self.k_fit_p[2])
        else:
            return False

    def eval_fit_cr(self, x:float):
        'Returns the cr-fit value'
        if self.fit_done:
            return self.fit_f_cr(x, self.cr_fit_p[0], self.cr_fit_p[1])
        else:
            return False

    def fit(self):
        'Fit CR(slit) and k(slit)'
        if len(self.k)<4 :
            print("Error, needs at least 4 blanket configurations for the fit!")
            return False
        # Fit k
        result_k = curve_fit(self.fit_f_k, np.array(self.slits), np.array(self.k), p0=self.k_fit_p, \
            sigma=np.array(self.kerr), absolute_sigma=True, epsfcn=0.0001, full_output=True )
        self.k_fit_p = result_k[0]
        # Get reduced Chi2
        self.k_chi2 = (result_k[2]['fvec']**2).sum() / (len(result_k[2]['fvec'])-len(result_k[0]))
        # Fit CR
        result_cr = curve_fit(self.fit_f_cr, np.array(self.slits), np.array(self.cr), p0=self.cr_fit_p, \
            sigma=np.array(self.crerr), absolute_sigma=True, epsfcn=0.0001, full_output=True )
        self.cr_fit_p = result_cr[0]
        # Get reduced Chi2
        self.cr_chi2  = (result_cr[2]['fvec']**2).sum() / (len(result_cr[2]['fvec'])-len(result_cr[0]))
        self.fit_done = True

    def fit_results(self) -> str:
        'Print the fit data'
        repr = ""
        if self.fit_done :
            repr += "Core K\n"
            repr += "p0 = %r, p1 = %r, p2 = %r\n" % (self.k_fit_p[0], self.k_fit_p[1], self.k_fit_p[2])
            repr += "chi2 = %r\n" % self.k_chi2
            repr += "CR\n"
            repr += "p0 =  %r, p1 = %r\n" % (self.cr_fit_p[0], self.cr_fit_p[1])
            repr += "chi2 = %r\n" % self.cr_chi2
        else:
            repr += "Error, the lattice was not fitted!\n"
        return repr

    def plot(self, plot_file:str='plot.pdf'):
        'Make plot of K and cr'
        if len(self.k)<2 :
            print("Error, needs at least 2 blanket configurations to plot!")
            return False
        # Plot data points
        plt.errorbar(self.slits, self.k,  self.kerr,  color="red",  marker="+", ls="none")
        plt.errorbar(self.slits, self.cr, self.crerr, color="blue", marker="+", ls="none")
        plt.title( r"Core [cm]: r$_{fuel}$=%5.1f"%self.c.rfuel + r", r$_{core}$=%5.1f"%self.c.rcore + r", z$_{core}$=%5.1f"%self.c.zcore \
            + "\n" + "Lattice [cm]: l=%5.2f"%self.c.l + ", r1=%5.3f"%self.c.r1 + ", r2=%5.3f"%self.c.r2 + \
            ", r3=%5.3f"%self.c.r3  )
        plt.xlabel("slit width [cm]")
        plt.ylabel(r"k$_{eff}$ red, CR blue")
        plt.grid(True)
        if self.fit_done :      # Plot fit lines
            x = np.linspace(min(self.slits), max(self.slits), num=200, endpoint=True)
            plt.plot(x, self.eval_fit_k(x), color="#ee3030", ls="-", \
                label="p0 = %6.4f" % self.k_fit_p[0] + ", p1 = %6.4f" % self.k_fit_p[1] + \
                ", p2 = %6.4f" % self.k_fit_p[2])
            plt.plot(x, self.eval_fit_cr(x),color="#3030ee", ls="-", \
                label="p0 = %6.4f" % self.cr_fit_p[0] + ", p1 = %6.4f" % self.cr_fit_p[1] )
        plt.legend(loc="best", fontsize="medium", title="Parameters of the polynomial fits")
        if plot_file == None:
            plt.show()
        else:
            plt.savefig(plot_file, bbox_inches='tight')
        plt.close()

    def get_k_at_cr(self, cr:float=1.0, extrapolate:bool=True) -> float:
        'Find K corresponding to particular CR using fits'
        if self.get_slit_at_cr(cr, extrapolate) != False:
            return self.eval_fit_k(self.get_slit_at_cr(cr, extrapolate))
        else:
            return -1.0

    def get_l4_at_cr(self, cr:float=1.0, extrapolate:bool=True) -> float:
        'Find l4 [cm] corresponding to particular CR using fits'
        if self.get_slit_at_cr(cr, extrapolate) != False:
            return self.c.l - self.get_slit_at_cr(cr, extrapolate)  # Return corresponding l4
        else:
            return -1.0

    def get_slit_at_cr(self, cr:float=1.0, extrapolate:bool=True) -> float:
        'Find l4 [cm] corresponding to particular CR using fits'
        if self.fit_done :          # extrapolate=1|0: do|not use extrapolated slitss
            # Find slit width for the required cr from the fit
            my_slit = (cr - self.cr_fit_p[0]) / self.cr_fit_p[1]
            if  my_slit < min(self.slits) or my_slit > max(self.slits) :  # Extrapolation warning
                print("Wrn, slit width for CR=%6.4f"%cr, "is %7.4f"%my_slit," - out of interp. range!")
                if not extrapolate:
                    return False
            return my_slit          # Return corresponding slit width based on fit functions
        else:
            print("Error, no fit data found!")
            return -1.0

if __name__ == '__main__':
    print("This is a module to find keff/CR relation for a ORNL-4528 core.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    s = CoreSlitSearch()    # Constructor of the slit search class
    s.set_nominal_core()    # Set nominal core using default parameters
    icr = 5  # how far, in lattice locations, are the 6 additional control rods
    s.c.control_rods=[(0,0), (icr,-icr), (icr,0), (0,icr), (-icr,0), (0,-icr), (-icr,icr)]
    s.c.crod_state  =[ 0, 0, 0, 0, 0, 0, 0 ]   # All control rods out
    s.set_path_from_geometry()      # Make sure searches run in separate directories
    s.build_cores()                 # Build geometries for all cores
    s.qsub_cores()                  # Submit all jobs
    s.get_calculated_values()       # Get all calculated data
    s.fit()                         # Make keff(slit) and CR(slit) fits
    print(s.fit_results())          # Print fit results
    s.plot('./plots/core-results.png')   # Plot fit results
    print("k and l4 at CR=1.00",s.get_k_at_cr(1.00), s.get_l4_at_cr(1.00))
    print("k and l4 at CR=1.01",s.get_k_at_cr(1.01), s.get_l4_at_cr(1.01))
    print("k and l4 at CR=1.02",s.get_k_at_cr(1.02), s.get_l4_at_cr(1.02))

