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
/// \file MakeMUONSingleRecoParam.C
/// \brief Macro to set reconstruction parameters and put them in the OCDB
///
/// \author Philippe Pillot, SUBATECH

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONRecoParam.h"
#include "AliMUONCDB.h"

#include "AliCDBManager.h"
#include "AliRecoParam.h"

#include <Riostream.h>

#endif


//-----------------------------------------------------------------------
void MakeMUONSingleRecoParam(Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity(),
			     AliRecoParam::EventSpecie_t eventType = AliRecoParam::kLowMult)
{
  /// set the reconstruction parameters and store them in the OCDB ($ALICE_ROOT/MUON/Calib/RecoParam/).
  /// - make a CDB entry for the run range [startRun, endRun]
  /// - "eventType" specifies the set of parameters to be stored
  
  // init CDB
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(startRun);
  
  // choose desired set of parameters
  AliMUONRecoParam* param;
  switch (eventType) {
  case AliRecoParam::kLowMult: // set of parameters for p-p runs
    param = AliMUONRecoParam::GetLowFluxParam();
    break;
  case AliRecoParam::kHighMult: // set of parameters for Pb-Pb runs
    param = AliMUONRecoParam::GetHighFluxParam();
    break;
  case AliRecoParam::kCosmic: // set of parameters for cosmic runs
    param = AliMUONRecoParam::GetCosmicParam();
    break;
  default: // unknown species
    cout<<"No set of parameters defined for the desired event type! Exiting..."<<endl;
    return;
    break;
  }
  param->Print("FULL");
  
  // save RecoParam in CDB
  AliMUONCDB cdb;
  cdb.WriteToCDB(param, "MUON/Calib/RecoParam", startRun, endRun, "reconstruction parameters for MUON", "Philippe Pillot");
  
}

