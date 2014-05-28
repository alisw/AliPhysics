//------------------------------------------------------------------------------//
   EvtGen is a particle decay simulator specifically designed for the needs of 
   B-physics studies. Complete manual of EvtGen is available at:               
   http://robbep.home.cern.ch/robbep/EvtGen/GuideEvtGen.pdf. Further information
   can be found on the EvtGen official web site  http://evtgen.warwick.ac.uk/src/          
   The EvtGen release installed here is taken from 
   http://svnweb.cern.ch/guest/evtgen/tags/R01-03-00
//------------------------------------------------------------------------------//

[23/05/2014] 
 1)This module contains: 

   - EvtGen original code organized in four directories:
     EvtGenBase-EvtGenModels-EvtGen-EvtGenExternal. 
     EvtGenModels directory contains all models avilable in EvtGen to decay 
     particles. EvtGenExternal contains the "interface" class with the external
     packages
 
   - The interface classes AliGenEvtGen and AliDecayerEvtGen for AliRoot users:
     the methods to decay particles are in AliDecayerEvtGen which represents 
     the implementation of AliDecayer using EvtGen package

 2)The external packages, interfaced with EvtGen by EvtGenExternal class,
   are used to generate specific particle decays and  
   are located in the TEvtGen subdirectories (except for Pythia8):

   - HepMC (http://lcgapp.cern.ch/project/simu/HepMC/download/HepMC-2.06.08.tar.gz)
     The HepMC package is an object oriented event record written in C++ for High Energy 
     Physics Monte Carlo Generators.    
 
   - Tauola (http://tauolapp.web.cern.ch/tauolapp/resources/TAUOLA.1.1.4/TAUOLA.1.1.4.tar.gz)
     C++ interface of Tauola code, specifically designed to generate tau decays. 

   - Photos (http://photospp.web.cern.ch/photospp/resources/PHOTOS.3.54/PHOTOS.3.54.tar.gz)
     C++ interface of Photos code, specifically designed to generate Final State Radiation.

   - Pythia8 (available in $ALICE_ROOT/PYTHIA8/) 


 3)The "default" decay table with all particles and all decay channels 
   can be found in TEvtGen/EvtGen/DECAY.DEC. The definition of all 
   particles is done in the table "evt.pdl". Both tables (DECAY.DEC and evt.pdl)
   are loaded during the initialization of EvtGen. Several decay tables are available 
   in the directory TEvtGen/EvtGen/DecayTable and are used to generate "forced" decay modes.
   All parameters, models and BR used in these tables are not modified with respect 
   to the official EvtGen release. 

 4)Polarization is still to be implemented in the interface classes.

 5)Below an example of the usage of EvtGen in AliRoot to decay beauty particles. 
   The following lines should be added in the Config.C:

 //----- load libraries
  gSystem->Load("libHepMC.so");
  gSystem->Load("libTauola.so");
  gSystem->Load("libpythia8.so");
  gSystem->Load("libPhotos.so");
  gSystem->Load("libEvtGen");
  gSystem->Load("libEvtGenExternal");
  gSystem->Load("libTEvtGen");

 //----- declare an AliGenCocktail 
  AliGenCocktail *generCock=new AliGenCocktail();
  generCock->UsePerEventRates();

 //----- declare Pythia configuration: switch-off beauty decays in Pythia
  AliGenPythia *pythia= .....
  .... 
  pythia->SetForceDecay(kNoDecayBeauty);
 
 //----- declare EvtGen configuration and put the two generators in the cocktail
  AliGenEvtGen *gene = new AliGenEvtGen();
  gene->SetForceDecay(kBJpsiDiElectron);
  gene->SetParticleSwitchedOff(AliGenEvtGen::kBeautyPart);
  generCock->AddGenerator(pythia, "Pythia", 1.);
  generCock->AddGenerator(gene,"gene",1.);
