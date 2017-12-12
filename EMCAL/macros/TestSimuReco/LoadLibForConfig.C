///
/// \file LoadLibForConfig.C
/// \ingroup EMCAL_TestSimRec
/// \brief Macro loading needed libraries for simulation.
///
/// Example macro to be executed right before TestEMCALSimulation.C in case of Root6
/// Cppy paste from LoadLibs() method in the Config.C
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS). 
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TSystem.h>

#endif

///
/// Main execution method
///
void LoadLibForConfig()
{
  // Load Pythia related libraries                                                                
  gSystem->Load("liblhapdf");      // Parton density functions                                 
  gSystem->Load("libEGPythia6");   // TGenerator interface                                     
  gSystem->Load("libpythia6");     // Pythia                                                   
  gSystem->Load("libAliPythia6");  // ALICE specific
  
  // Load Geant3 related libraries                                                                
  gSystem->Load("libgeant321");
}
