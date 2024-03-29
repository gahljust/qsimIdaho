# qsimIdaho

This is the modified qsim-master directory, with changes made to suit the SLAC test beam data

Justin Gahley
gahljust@isu.edu
November 14, 2019

qsim installation and running instructions

Seamus Riordan
sriordan@physics.umass.edu
September 26, 2013

Updated June 26, 2015

---------------------------------------------------
Packages to build this:

>=Geant4.10
cmake > 2.6
root

---------------------------------------------------

Instructions:

To build, create the directory you would like to
build in, say build

mkdir build
cd build
cmake <path to qsim>
make


It just needs to be downloaded and for the example, present in the directory
you are running in.  Otherwise, specify it in the macro as in the
exmaple provided in macros/  An example macro will be found in

macros/runexmaple.mac

To run in batch mode, execute with a macro such as:

./qsim runexample.mac

Ensure that all macros include a /run/initialization command or else they will
not work. 

==== Visualization ===============

Visualization macros are found in vis/

To run, execute

./qsim

which should bring up a graphical command interface

To see the geometry:

/control/execute macros/vis.mac


---------------------------------------------------

==== CLI User Commands ===========

Using the Geant4 CLI it is possible to pass commands to modify behavior
and utilize the vis.mac macro from the command line.
These are all visible from the menu on the left.


---------------------------------------------------

==== Operational Mode Switches ===

There are CLI User Commands that allow the user to change the stand design and
incident particle characteristics to test different configurations and 
experimental expectations. These commands must be passed before a visualization
macro is used or the initialization command have been passed.

*******************
*** Source mode *** 
*******************

Set by: /qsim/fSourceMode <0, 1, 2>

0 = cosmic muons
Generates primary particles following cosmic muon angular distribution and energy spectrum.  
Energy spectrum obtained from fit to PDG data for muons with 0 deg incidence (good to 25% out to 36 GeV).  
Note that this does not automatically change the primary particle type to muons; this must be set by "/gun/particle mu-"

1 = beam
Generates perfectly straight, monoenergetic beam.
Current implementation generates particles at pinpoint, but beam spot size can be changed in qsimPrimaryGeneratorAction.cc 
Energy of beam can be changed in qsimPrimaryGeneratorAction.cc

2 = PREX
Generates 1.063 GeV particles following position and angular distribution observed at VDCs during PREX-I.
The z position of primary vertex can be changed in qsimPrimaryGeneratorAction.cc, to effectively move detector closer/farther from VDC.
The distributions are stored in file primaryDistribution.root (copied to build directory when qsim is made), which has 2e6 events.
 
******************
*** Stand mode *** 
******************

Set by: /qsim/fStandMode <0, 1>

0 = beam/PREX
Detector only.

1 = cosmic
Detector, top/bottom scintillators, and lead.
Scintillator size/separation and lead size can be adjusted in qsimDetectorConstruction.cc


*********************
*** Detector mode *** 
*********************

Set by: /qsim/fDetMode <0, 1>

0 = PREX-I detector
Detector geometry follows PREX-I detector design.
Detector rotated such that primary particles in beam mode are incident on quartz at exactly 45 deg.

1 = PREX-II detector prototype
Detector geometry follows PREX-II prototype.
Detector rotated such that primary particles in beam mode are incident on quartz at exactly 90 deg.
Current implementation follows Stony Brook prototype (small quartz-PMT separation, no light guide).


**************************
*** Configuration mode *** 
**************************

Set by: /qsim/gConfMode <0, 1, 2, 3, 4>

For showermax only

0 = 1 Quartz

1 = 1 Stack

2 = 2 Stack

3 = 3 Stack

4 = 4 Stack

*** A NOTE ON OPTICAL PROPERTIES ***

Index of refraction (quartz): 
Specification sheet for Heraeus Spectrosil 2000 provides >25 data points for n(E).
Fit with polynomial.

Absorption length (quartz):
Specification sheet for Heraeus Spectrosil 2000 provides only 2 data points for L(E).
Current functional form of L(E) in qsim is of unknown origin and is inconsistent with Heraeus data points.

Reflectivity (mirror):
Currently defined as a function of photon energy only (and possibly incorrect).
Needs to be defined as a function of photon energy AND photon angle (supposedly possible in Geant 4.10 but not yet implemented in qsim).
