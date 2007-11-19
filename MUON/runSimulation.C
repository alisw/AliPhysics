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

// Macro extracted from MUON test script
// By Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliSimulation.h"
#include <TRandom.h>
#endif

void runSimulation(int run, int seed, int nevents, const char* config)
{ 
// Uncoment following lines to run simulation with local residual mis-alignment
// (generated via MUONGenerateGeometryData.C macro)
// AliCDBManager* man = AliCDBManager::Instance();
// man->SetDefaultStorage("local://$ALICE_ROOT");
// man->SetSpecificStorage("MUON/Align/Data","local://$ALICE_ROOT/MUON/ResMisAlignCDB");
  AliSimulation MuonSim(config);
  MuonSim.SetRunNumber(run);
  MuonSim.SetSeed(seed);
  MuonSim.SetMakeTrigger("MUON");
  MuonSim.SetWriteRawData("MUON","raw.root",kTRUE);
  MuonSim.Run(nevents);
  //gObjectTable->Print();
}
