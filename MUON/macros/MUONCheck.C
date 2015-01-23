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
/// \file MUONCheck.C
/// \brief Macro for data quality control
///  
/// \author: Frederic Yermia, INFN Torino

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONCheck.h"

#include "AliCDBManager.h"

#include <TSystem.h>

#endif

void MUONCheck(Int_t firstEvent, Int_t lastEvent,
               TString fileNameSim="$ALICE_ROOT/MUON/test_out/galice_sim.root",
               TString fileName="$ALICE_ROOT/MUON/test_out/galice.root",
               TString esdsFileName="$ALICE_ROOT/MUON/test_out/AliESDs.root",
               TString outDir="$ALICE_ROOT/MUON/test_out/DataQualityControl")
{
   // Set default CDB storage
   AliCDBManager* man = AliCDBManager::Instance();
   man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

   gSystem->Load("libMUONevaluation");

   AliMUONCheck* check
     = new AliMUONCheck(fileName.Data(), fileNameSim.Data(), esdsFileName.Data(),
                        firstEvent, lastEvent, outDir.Data());
      
   check->CheckESD();
   check->CheckKine();
   check->CheckTrackRef();
   check->CheckOccupancy();

   // delete check;  
}

