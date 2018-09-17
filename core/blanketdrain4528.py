#!/usr/bin/python3
#
# Script to test for blanket drain criticality safety and shutdown margin
#
# Ondrej Chvala, ochvala@utk.edu

from ornl4528core import CoreGen
import time

# Core parameters
hpitch:float    = 6.045
r1:float        = 2.550
r2:float        = 3.758
r3:float        = 4.534
l4:float        = 5.697
tempC:float     = 700.
rgr_scale:float = 0.90
rfuel:float     = 200.0
rcore:float     = 230.0
zcore:float     = 400.0
refl_ht:float   = 100.0
n_deckname:str  = 'ORNL-4528 deck for the normal core'
s_deckname:str  = 'ORNL-4528 deck for the scrammed core'
d_deckname:str  = 'ORNL-4528 deck for the core with a drained blanket'
t_deckname:str  = 'ORNL-4528 deck for the scrammed/stopped core with no flow'

# Define control rod positions and states
icr:int         = 8 # This position gets the largest reactivity effects
control_rods    = [(0,0), (icr,-icr), (icr,0), (0,icr), (-icr,0), (0,-icr), (-icr,icr)]
n_crod_state    = [ 0, 0, 0, 0, 0, 0, 0 ]   # All rods out
s_crod_state    = [ 1, 1, 1, 1, 1, 1, 1 ]   # All control clusters in
d_crod_state    = [ 2, 2, 2, 2, 2, 2, 2 ]   # All floating rods in
t_crod_state    = [ 3, 3, 3, 3, 3, 3, 3 ]   # All rods in

# Create ORNL-4528 cores
normal_core     = CoreGen(hpitch, r1, r2, r3, l4, \
       tempC, rgr_scale, rfuel, rcore, zcore, refl_ht, n_deckname)
scrammed_core   = CoreGen(hpitch, r1, r2, r3, l4, \
       tempC, rgr_scale, rfuel, rcore, zcore, refl_ht, s_deckname)
drained_core    = CoreGen(hpitch, r1, r2, r3, l4, \
       tempC, rgr_scale, rfuel, rcore, zcore, refl_ht, d_deckname)
stopped_core    = CoreGen(hpitch, r1, r2, r3, l4, \
       tempC, rgr_scale, rfuel, rcore, zcore, refl_ht, t_deckname)

# List of my cores
my_cores  = [ normal_core, scrammed_core, drained_core, stopped_core]

# Assign control rod lattice positions, make sure to make geometry and mesh plots
for c in my_cores:
    c.geomplots     = True
    c.meshplots     = True
    c.control_rods  = control_rods

# Assign control rod states
normal_core.crod_state      = n_crod_state
scrammed_core.crod_state    = s_crod_state
drained_core.crod_state     = d_crod_state
stopped_core.crod_state     = t_crod_state

# Set sub-directories to run each case
normal_core.deck_path       = './00_normal'
scrammed_core.deck_path     = './01_scrammed'
drained_core.deck_path      = './02_drained'
stopped_core.deck_path      = './03_stopped'

# Drain the drained core
drained_core.blanket_drained = True

# Build all cores and save input decks
for c in my_cores:
    c.save_deck()

# Run all cases
for c in my_cores:
    if not c.get_calculated_values():
        c.run_deck()

# Wait for Serpent jobs and get results
cores_done:int = 0
while cores_done < len(my_cores):
    cores_done = 0
    for c in my_cores:
        if(c.k < 0.0):
            if c.get_calculated_values():
                cores_done += 1
                print("[DEBUG] Got results for ", c.deck_path)
        else: # Core already read in
            cores_done += 1
    if cores_done < len(my_cores):
        print("[DEBUG] ", cores_done, " done, sleeping ...")
        time.sleep(20)  # Wait a minute for Serpent ...

# Calculate reactivities [pcm]
rho_n:float  = 1e5*( normal_core.k   - 1.0) / normal_core.k
rho_d:float  = 1e5*( drained_core.k  - 1.0) / drained_core.k
rho_s:float  = 1e5*( scrammed_core.k - 1.0) / scrammed_core.k
rho_t:float  = 1e5*( stopped_core.k  - 1.0) / stopped_core.k

# Print results
print("\nReactivity [pcm] for the studied cases")
print("normal       : %8.1f " % rho_n)
print("drained      : %8.1f " % rho_d)
print("scrammed     : %8.1f " % rho_s)
print("stopped      : %8.1f " % rho_t)

print("\nReactivity difference [pcm] from normal condition")
print("drained      : %8.1f " % (rho_d-rho_n))
print("scrammed     : %8.1f " % (rho_s-rho_n))
print("stopped      : %8.1f " % (rho_t-rho_n))



''' *** Results for different CR absorber choices ***
--- Enriched B4C ---
Reactivity [pcm] for the studied cases
normal       :   2927.7 
drained      :  41704.2 
scrammed     :  -5175.4 
stopped      :  -5573.0 

Reactivity difference [pcm] from normal condition
drained      :  38776.5 
scrammed     :  -8103.1 
stopped      :  -8500.7 

--- Natural B4C ---
Reactivity [pcm] for the studied cases

Reactivity difference [pcm] from normal condition

--- Boron metal ---
Reactivity [pcm] for the studied cases

Reactivity difference [pcm] from normal condition

'''
