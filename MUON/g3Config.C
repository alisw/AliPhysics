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
/// \file g3Config.C
/// \brief Configuration macro for MUON spectrometer Geant3 simulation

void Config(const char* directory="", 
            const char* option="param", 
           const char* digitstore="AliMUONDigitStoreV2S",
           bool forEmbedding=kFALSE)
{
  cout << "Running g3Config.C ... " << endl;

  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/MUON/commonConfig.C");
  commonConfig(directory, digitstore, forEmbedding);

  // Load Geant3 + Geant3 VMC libraries
  //
#if defined(__CINT__)
    gSystem->Load("libgeant321");
#endif

  // Create TGeant3
  //  
  new  TGeant3TGeo("C++ Interface to Geant3");

  // AliRoot event generator
  // (it has to be created after MC, as it may use decayer via VMC)
  //
  gROOT->LoadMacro("$ALICE_ROOT/MUON/genTestConfig.C");
  genConfig(option);

  // From external file
  //
  //gROOT->LoadMacro("$ALICE_ROOT/MUON/genExtFileConfig.C");
  //genConfig();


  cout << "Running g3Config.C finished ... " << endl;
}  
