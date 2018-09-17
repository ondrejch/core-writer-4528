#!/usr/bin/python3
#
# Class to find relations of reactor physics parameters to
# slit with, given the rest of lattice geometry.
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

from latgen4528 import LatGen, std_bf, min_blanket_width

my_debug = True

def nominal_l4(l:float, r1:float) -> float:
    'Return nominal l4 using standard blanket fraction'
    return math.sqrt(std_bf * l**2 )

class SlitCritSearch(object):
    'Establish k_inf(slit_width) and CR(slit_width)'
    def __init__(self, my_path:str='./.runs'):
        'Constructor with default slit criticality search values'
        self.relslit_min:float = 0.80   # Minimum relative slit width
        self.relslit_max:float = 1.20   # Maximum relative slit width
        self.relslit_N:int     = 11     # Number of steps in slit width
        self.latlist           = []     # List of lattices

        # Sensible defaults
        self.l:float        = 7.09      # Hex lattice apothem size / half pitch [cm]
        self.r1:float       = 1.94      # Radius of the inner fuel channel log [cm]
        self.r2:float       = 2.86      # Radius of the first graphite ring [cm]
        self.r3:float       = 3.45      # Radius of the outer fuel channel [cm]
        self.l4:float       = 6.83      # Apothem of the graphite hex [cm]

        self.main_path:str  = my_path   # Main directory to run cases
        self.deck_file_name = 'msbr.inp'# Serpent lattice deck name
        self.sleep_sec:int  = 60        # Sleep timer between results read attempts [s]

        self.slits      = array('d')    # Slit sizes
        self.k          = array('d')    # ANALYTICAL_KEFF
        self.kerr       = array('d')    # ANALYTICAL_KEFF Error
        self.cr         = array('d')    # CONVERSION_RATIO
        self.crerr      = array('d')    # CONVERSION_RATIO Error
        self.fit_done   = False         # Fit flag
        self.k_fit_p    = [1.0, -1.0, 0.1]  # k fit parametrs
        self.cr_fit_p   = [1.0,  1.0]       # CR fit parametrs
        self.k_chi2:float  = -0.1       # Reduced \Chi^2 of the k(slit) fit
        self.cr_chi2:float = -0.1       # Reduced \Chi^2 of the CR(slit) fit

    def set_path_from_geometry(self):
        'Sets path to directory to run cases based on slit geometry'
        self.main_path = self.main_path + "/" + "%08.5f"%self.l + "/%08.5f"%self.r1 + "/%08.5f"%self.r2 + \
            "/%08.5f"%self.l4           # r3 omitted since fixed by flow speed up = down

    def nominal_l4(self) ->float:
        'Return nominal l4 using standard blanket fraction'
        return math.sqrt(std_bf * self.l**2 )

    def build_lattices(self):
        'Builds lattices for the slit study case'
        rel_slit_array = np.linspace(self.relslit_min, self.relslit_max, self.relslit_N)
        if my_debug:
            print("[DEBUG] Nominal l4 ", self.nominal_l4())
        nom_slit:float = self.l - self.nominal_l4()  # nominal slit size using standard blanket fraction
        for rslit in rel_slit_array:
            slit:float = rslit * nom_slit
            if slit < min_blanket_width:
                continue            # skip if slit is too small
            my_l4:float = self.l - slit
            if my_debug:
                print("[DEBUG] Building lattice, rslit {rslit}, slit {slit} cm".format(**locals()))
            lat = LatGen(self.l, self.r1, self.r2, self.r3, my_l4)
            lat.deck_path = self.main_path +'/'+ str(lat.l4)
            self.latlist.append( lat )

    def qsub_lattices(self, rerun=False):
        'Run lattice simulations, rerun if rerun=True'
        for lat in self.latlist:
            try:                    # Create directory for each job
                os.makedirs(lat.deck_path)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise Exception("Unable to create directory: ", lat.deck_path, str(e))
                    return False
            lat.deck_file_name = self.deck_file_name    # Set Serpent input file name
            if (rerun == False) and lat.get_calculated_values():
                pass # Avoid recalculating existing data
            else:
                if my_debug:
                    print("[DEBUG] Submitting deck", lat, " in ", lat.deck_path)
                lat.run_deck()

    def get_calculated_values(self, reread=True):
        'Fill K and CR values from Serpent results file'
        lats_done:int = 0
        if reread:
            for lat in self.latlist:        # Flag all lattices for re-read
                lat.k = -0.1
        while lats_done < len(self.latlist):
            lats_done = 0
            for lat in self.latlist:
                if(lat.k < 0.0):
                    if lat.get_calculated_values():
                        lats_done += 1
                        if my_debug:
                            print("[DEBUG] Got results for ", lat.deck_path+'/'+lat.deck_file_name)
                else: # Lattice already read in
                    lats_done += 1
            if lats_done < len(self.latlist):
                if my_debug:
                    print("[DEBUG] ", lats_done, " done, sleeping ...")
                time.sleep(self.sleep_sec)  # Wait a minute for Serpent ...
        for lat in self.latlist:    # Fill arrays for fit
            self.slits.append(lat.slit())
            self.k.append(lat.k)
            self.kerr.append(lat.kerr)
            self.cr.append(lat.cr)
            self.crerr.append(lat.crerr)

    def fit_f_k(self, x:float, p0:float, p1:float, p2:float) -> float:
        'Fit function for Lattice.fit() of KEFF'
        return p2*x**2 + p1*x + p0

    def fit_f_cr(self, x:float, p0:float, p1:float) -> float:
        'Fit function for Lattice.fit() of CR'
        return p1*x + p0

    def eval_fit_k(self, x:float):
        'Returns k-fit value'
        if self.fit_done:
            return self.fit_f_k(x, self.k_fit_p[0], self.k_fit_p[1],self.k_fit_p[2])
        else:
            return False

    def eval_fit_cr(self, x:float):
        'Returns cr-fit value'
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
            repr += "Lattice K\n"
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
        plt.title("MSBR lat [cm]:l=%5.2f"%self.l + ", r1=%5.3f"%self.r1 + ", r2=%5.3f"%self.r2 + \
            ", r3=%5.3f"%self.r3)
        plt.xlabel("slit width [cm]")
        plt.ylabel(r"k$_{inf}$ red, CR blue")
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
            my_file = self.main_path + '/' + plot_file
            if not os.path.exists(os.path.dirname(my_file)):
                os.makedirs(os.path.dirname(my_file))
            plt.savefig(my_file, bbox_inches='tight')
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
            return self.l - self.get_slit_at_cr(cr, extrapolate) # Return corresponding l4
        else:
            return -1.0

    def get_slit_at_cr(self, cr:float=1.0, extrapolate:bool=True) -> float:
        'Find slit [cm] corresponding to particular CR using fits'
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
    print("This is a module to find keff/CR relation for a LFTR lattice.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    s = SlitCritSearch()        # Call constructor
    s.relslit_N = 21            # Change defaults
    s.set_path_from_geometry()  # Make sure lattices run in proper (and different) directories
    s.build_lattices()          # Build lattice objects
    s.qsub_lattices()           # Run them on the cluster
    s.get_calculated_values()   # Harvest results
    s.fit()                     # Make keff(slit) and CR(slit) fits
    print(s.fit_results())      # Print fit results
    s.plot('./plots/results.png')   # Plot fit results
    print("k and l4 at CR=1.00",s.get_k_at_cr(1.00), s.get_l4_at_cr(1.00))
# k and l4 at CR=1.00 1.04125042367 6.84842491184
    print("k and l4 at CR=0.97",s.get_k_at_cr(0.97), s.get_l4_at_cr(0.97))
# k and l4 at CR=0.97 1.05677315491 6.85672775238


