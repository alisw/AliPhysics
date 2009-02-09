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
/// \file MakeMUONRecoParamArray.C
/// \brief Macro to set reconstruction parameters and put them in the OCDB
///
/// \author Philippe Pillot, SUBATECH

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONRecoParam.h"
#include "AliMUONCDB.h"

#include "AliCDBManager.h"
#include "AliRecoParam.h"

#include <TObjArray.h>
#include <TIterator.h>

#include <Riostream.h>

#endif


//-----------------------------------------------------------------------
void MakeMUONRecoParamArray(Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity(),
		       AliRecoParam::EventSpecie_t defaultParam = AliRecoParam::kLowMult)
{
  /// set the reconstruction parameters and store them in the OCDB ($ALICE_ROOT/MUON/Calib/RecoParam/).
  /// - make a CDB entry for the run range [startRun, endRun]
  /// - "defaultParam" specifies the parameters to be used as default
  
  // init CDB
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(startRun);
  
  // set RecoParams
  AliMUONRecoParam *param;
  TObjArray recoParams;
  
  // set parameters for p-p runs
  param = AliMUONRecoParam::GetLowFluxParam();
  recoParams.AddLast(param);
  
  // set parameters for Pb-Pb runs
  param = AliMUONRecoParam::GetHighFluxParam();
  recoParams.AddLast(param);
  
  // set parameters for cosmic runs
  param = AliMUONRecoParam::GetCosmicParam();
  recoParams.AddLast(param);
  
  // identify default parameters (exit if identification failed)
  Bool_t defaultIsSet = kFALSE;
  TIter next(recoParams.MakeIterator());
  while ( (param = static_cast<AliMUONRecoParam*>(next())) ) {
    if (param->GetEventSpecie() == defaultParam) {
      param->SetAsDefault();
      defaultIsSet = kTRUE;
    }
    param->Print("FULL");
  }
  if (!defaultIsSet) {
    cout<<"The default reconstruction parameters are not set! Exiting..."<<endl;
    return;
  }
  
  // save RecoParam in CDB
  AliMUONCDB cdb;
  cdb.WriteToCDB(&recoParams, "MUON/Calib/RecoParam", startRun, endRun, "reconstruction parameters for MUON", "Philippe Pillot");
  
}

