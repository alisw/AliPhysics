/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/// \ingroup macros
/// \file g4Config.C
/// \brief Configuration macro for MUON spectrometer Geant4 simulation
///
/// \author I. Hrivnacova, IPN Orsay
 	
void Config(const char* directory="", 
            const char* option="param", 
            const char* digitstore="AliMUONDigitStoreV2S",
            bool forEmbedding=kFALSE)
{
  cout << "Running g4Config.C ... " << endl;

  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/MUON/commonConfig.C");
  commonConfig(directory, digitstore, forEmbedding);

  // Load Geant4 + Geant4 VMC libraries
  //
  if (gClassTable->GetID("TGeant4") == -1) {
    // Load Geant4 libraries 
    if (!gInterpreter->IsLoaded("$ALICE/geant4_vmc/examples/macro/g4libs.C")) {
      gROOT->LoadMacro("$ALICE/geant4_vmc/examples/macro/g4libs.C");
      gInterpreter->ProcessLine("g4libs()");
    }
  }    

  // Create Geant4 VMC
  //  
  TGeant4 *geant4 = 0;
  if ( ! gMC ) {
    TG4RunConfiguration* runConfiguration 
      = new TG4RunConfiguration("geomRoot", 
                                "QGSP_BERT_EMV", 
                                "specialCuts+stepLimiter",
                                 true);

    geant4 = new TGeant4("TGeant4", "The Geant4 Monte Carlo", runConfiguration);
    cout << "Geant4 has been created." << endl;
  } 
  else {
    cout << "Monte Carlo has been already created." << endl;
  }  

  // Customization of Geant4 VMC
  //

  geant4->ProcessGeantCommand("/mcPhysics/rangeCuts 0.01 mm"); 
  geant4->ProcessGeantCommand("/mcVerbose/all 1");  
  geant4->ProcessGeantCommand("/mcTracking/loopVerbose 0");     
  geant4->ProcessGeantCommand("/mcTracking/skipNeutrino true");

  // Uncomment these lines when running with G4 native navigation
  // (geomRootToGeant4)
  //geant4->ProcessGeantCommand("/vgm/setNameSeparator /");
  //geant4->ProcessGeantCommand("/mcControl/accountAssemblies true");

  // Uncomment this line to get a detail info from each step 
  // geant4->ProcessGeantCommand("/tracking/verbose 1");  
  
  // More info from the physics list
  // the verbosity level is passed to all contained physics lists and their
  // physics builders
  //geant4->ProcessGeantCommand("/mcVerbose/composedPhysicsList 2");  
  
  // More info from optical processes
  //geant4->ProcessGeantCommand("/mcVerbose/opticalPhysicsList 3");  
  
  // More info from geometry building
  //geant4->ProcessGeantCommand("/mcVerbose/geometryManager 2");  

  // More info from setting geometry properties (in materials and surfaces)
  // for optical physics
  //geant4->ProcessGeantCommand("/mcVerbose/opGeometryManager 1");  
  
  // More info about regions construction 
  // and conversion of VMC cuts in cuts in range per regions 
  // geant4->ProcessGeantCommand("/mcVerbose/regionsManager 2");
  
  // Suppress verbose info from tracks which reached maximum number of steps
  // (default value is 30000)  
  //geant4->ProcessGeantCommand("/mcTracking/loopVerbose 0"); 
    
  // AliRoot event generator
  // (it has to be created after MC, as it may use decayer via VMC)
  //
  gROOT->LoadMacro("$ALICE_ROOT/MUON/genTestConfig.C");
  genConfig(option);

  // From external file
  //
  //gROOT->LoadMacro("$ALICE_ROOT/MUON/genExtFileConfig.C");
  //genConfig();

  cout << "Running g4Config.C finished ... " << endl;
}
