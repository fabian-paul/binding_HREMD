description
===========
This plugin implements Hamiltonian replica exchange MD [1] and
boost potentials [2]. Boost potentials are applied to the Lennard-
Jones interaction between two set of atoms (normally two molecules).

This plugin needs ACEMD pro.

plugin parameters
=================
* energyfreq: this is acutally not a plugin paremater but a normal
  ACEMD parameter. It specifies after how many integration steps an
  exchange will be attempted [might change in future versions]
* E: array of energy threshold parameters (one for each replica)
  Can be a whitespace-separated floats or a TCL list of floats.
  Refer to [2]. Energies are measured in kcal/mol.
* alpha: array of alpha parameters (one for each replica)
  Must be the same size as E. Refer to [2] and Fabians notes 
  (and future papers and thesis).
* restart "on" | "off": resume from a restart file
* restartfreq N: write restart files every N integration steps
* fullenergyfreq N: calculate the full matrix of exchange energies
  (every replica under every Hamiltonian) every N integration
  steps. (By default, only the energies that are necessary for
  the trial move are calculated) [Note: this might change, once
  the Chodera swapping scheme is implemented]
  This is also the interval with which the energy log file
  and the permutaion log file are written. Typically this
  is set to the same value as xtcfreq.

files
=====
* boost.c/.h  implementation of boost potentials
* hremd-boost.c/.h  implementation of HREMD that is specific to my
              implementation of boost potential.
* hremd.tcl   example run file for ACEMD
* restart.py  script to work around the bothersome implementation
              of ACEMD restart mechanism. Takes the number of the 
              last integration step N from a HREMD restart file 
              and truncates all xtc files at the frame N/2500.
              Adapt source to match your value of xtcfreq.
              [This might change, when ACEMD gets a sensible restart mechanism.]
* truncate-xtc.c used by restart.py; needs to be compiled using the 
              [XTC libs](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library)
              
files generated by hremd-boost
==============================

* perm.log: contains one (cumulative) permutation per line
  Let R_i be the value located at the i'th column. The the
  map from Hamiltonian index i to replica index is i->R_i.
  Initially the permuatation is the identity (1 2 3 4 5 6 ...)
  This is never written, because perm.log is updated after the
  (first) swap.
  
* N.*.xtc: replica trajectory. Trajectories are continuous in
           replica (i.e. in spatial coordinates).

* energy.log: time series of energy matrices. Successive matrices
              are separated by a line containing only one & character.
              Each matrix element e[i,j] is the energy of replica i
              evaluated under the Hamiltonian j. Units are kcal/mol.

inconveniences to be adressed in future versions
================================================

* In the currect versions the information about which atom
  belongs to which interaction partner is hard-coded into the
  plugin code. Future versions should introduce a plugin
  argument to pass the atom indices or use the beta/occupancy/
  chain fields of the pdb to assign this membership.

* In the current version, the temperature T is hardcoded into
  the program. Future versions should read its value via the
  ACEMD api.

* Currently much of the force calculation is done on the CPU.
  This might change (hopefully) with future version of ACEMD.
  An OpenMP implementation is also on the way.

* Future versions should implement the swapping scheme by Chodera
  (DOI 10.1063/1.3660669) There swapping is only done by one
  agent. The energies are send to a central agent which does all the
  swapping and updates the permutation. The permutation is then
  communicated to all members of the swarm that change their
  Hamiltonian accordingly. How should this communication be organized?

* The plugin executable is only found when it is in the working
  directory of ACEMD (or possibly in the global plugin directory).
  This is a bit inconveniant. Find out why.
  
* Right now the calculation of each trial energy contains a full
  calculation of the physical energy (because V'=f(V)). But when
  trying different energies (for the same confirguation), the 
  physical energy needs to be evaluated only once. Implement this!
  (minor performance improvement)
  
* We could also attempt exchanges every step, because the overhead
  is small. I have to change the calculation of the physical energy
  V such that the integration step can use V from the energy part
  of the last step.
  
* There are some problems with simulation restarts. ACEMD will always
  introduce some small error into the forces at every restart. This
  was explained by Matt and will hopefully be changed in the next versions:

    I can move the point at which plugin_init is called. It's arguably not in
    the right place just now.
    
    Matt
    
    On 15 August 2013 17:53, Fabian Paul <fab@physik.tu-berlin.de> wrote:
    
    >> Hi Fabian,
    
    >>> 2)
    
    >>> When I look into the log file after a restart, there is a line with the energies
    >>> of the zero'th step that comes even before the plugin is initialized. Do these
    >>> energies come from ACEMD's restart file or are they recalculated?
    
    >> ACEMD does one force evaluation which is used to set positions at dt-.5,
    >> for shake. They are evaluated using the coords in the restart file. They
    >> are printed out for diagnostic purposes, and should match the energies of
    >> the last step in the previous sim. plugin_init is called immediately
    >> after, before the first actual iteration.
    
    > Ok, but I think this would introduce some small error. In the last step of
    > the old simulation the Hamiltonians might have been permuted. But
    > immediately after the restart the plugin has no chance to restore the old
    > permutation of Hamiltonians. So this force is evaluated with the wrong
    > combination of coordinates and Hamiltonian.
    > So maybe it's better to exchange coordinates instead of Hamiltonian after
    > all? Or permute the topology files?

* There is some effect after that happens after a restart and which may look
  odd to the user: after a restart, the files *.perm.log and *.energy.log
  are not truncated at the last step but the file pointer is repositioned
  to the last step. This means that immedeately after a restart the number
  of lines in the *.perm.log and *.energy.log files will be larger than
  timestep/fullenergyfreq. As the simulation continues, the old data in
  these files will be overwritten and when the simulation passed the time
  step where the previous simulation was terminated, the files will 
  start to grow again. [I might fix this in the future using unistd.h's
  truncate().]

References
==========
[1] Yuji Sugita and Yuko Okamoto. Replica-exchange molecular dynamics
    method for protein folding. Chemical Physics Letters 314(1-2):141-151, 1999.

[2] Hamelberg et al. Accelerated molecular dynamics: A promising and 
    effecient simulation method for biomolecules. J. Chem. Phys., 120(24):11919, 2004.

author and license
==================
(C) 2014 Fabian Paul <fabian.paul@mpikg.mpg.de>

hremd-boost is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

hremd-boost is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with hremd-boost.  If not, see <http://www.gnu.org/licenses/>.
