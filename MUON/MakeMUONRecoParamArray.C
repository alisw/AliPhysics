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
#include "AliCDBManager.h"
#include "AliMUONCDB.h"

#include <Riostream.h>

#endif


//-----------------------------------------------------------------------
void MakeMUONRecoParamArray(Int_t startRun = 0, 
                            Int_t endRun = AliCDBRunRange::Infinity(),
                            const char* settings="ppIdeal")
{
  /// set the reconstruction parameters and store them in the OCDB ($ALICE_ROOT/OCDB/MUON/Calib/RecoParam/).
  ///
  /// - make a CDB entry for the run range [startRun, endRun]
  ///
  /// for the possible values of settings, please see AliMUONRecoParam::Create
  
  // init CDB
  AliCDBManager* man = AliCDBManager::Instance();
  
  if (!man->IsDefaultStorageSet()) 
  {
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");    
  }
  
  man->SetRun(startRun);
  
  TObjArray* recoParams = AliMUONRecoParam::Create(settings);
  
  if (recoParams)
  {
    // save RecoParam in CDB
    AliMUONCDB::WriteToCDB(recoParams, "MUON/Calib/RecoParam", startRun, endRun, 
                           "reconstruction parameters for MUON", "MakeMUONRecoParamArray $Id$");
  }
  
  delete recoParams;
}

