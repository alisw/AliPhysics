// $Id$

/*! 

\page README_fast MUON Fast simulation


The macro fastMUONGen.C allows to generate events using the 
AliGenMUONCocktailpp generator. In the current implementation, the 
generator is set up in order to force all generated particles to decay
within Pythia, including the muonic decay of pions and kaons:
gener->SetDecayModePythia(kAllMuonic); 
By changing the argument from "kAllMuonic" back to "kAll" one can
leave it up to GEANT to do these decays. 
In addition, we require that the "trigger muons" - whose number 
(single- or dimuon trigger or MB trigger),
specified as the third argument of fastMUONGen.C - must have their origin
before the hadron absorber in order to "simulate" GEANT:
gener->SetMuonOriginCut(-130.);

The output of this generation ("galice.root" and "Kinematic.root")
can then be processed by the macro fastMUONSim.C that produces a
root file, called "fastSim_pp.root". It uses the "fast generator" to
simulate the detector response.  The output file contains
TClonesArrays holding AliMUONTrackLight and AliMUONPairLight objects
that are built on the basis of the surviving muons.

Both macros can be run from fastSim.sh test script.

This chapter is defined in the READMEfast.txt file.

*/


