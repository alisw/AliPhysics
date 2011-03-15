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
/// \file MUONStatusMap.C
/// \brief Macro to check/test pad status and pad status map makers
///
/// \author Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONLogger.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONRecoParam.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVStore.h"
#include "AliMUONCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include "Riostream.h"
#endif

//AliMUONVStore* MUONStatusMap(const TString& cdbStorage = "alien://folder=/alice/data/2009/OCDB",
AliMUONVStore* MUONStatusMap(const TString& cdbStorage = "alien://folder=/alice/data/2011/OCDB",
                             Int_t runNumber=145000, Bool_t statusOnly=kTRUE)
{  

  AliCDBManager::Instance()->SetDefaultStorage(cdbStorage);
  AliCDBManager::Instance()->SetRun(runNumber);
  
  
  AliMUONCDB::LoadMapping();
  
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  
  AliCDBManager* man = AliCDBManager::Instance();
  
  man->SetDefaultStorage(cdbStorage.Data());
  
//  man->SetSpecificStorage("MUON/Calib/OccupancyMap","local://$ALICE_ROOT/OCDB");
//  man->SetSpecificStorage("MUON/Calib/OccupancyMap","alien://folder=/alice/cern.ch/user/l/laphecet/OCDB");
//  man->SetSpecificStorage("MUON/Calib/RejectList","alien://folder=/alice/cern.ch/user/l/laphecet/OCDB");
//  man->SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/cern.ch/user/l/laphecet/OCDB");

  AliMUONCalibrationData cd(runNumber);
  
  AliMUONPadStatusMaker statusMaker(cd);
  
  statusMaker.SetLimits(*recoParam);
  
  UInt_t mask = recoParam->PadGoodnessMask();

  //  delete recoParam;
  
  statusMaker.Report(mask);
  
  if ( statusOnly ) return statusMaker.StatusStore();
  
  const Bool_t deferredInitialization = kFALSE;
  
  AliMUONPadStatusMapMaker statusMapMaker(cd,mask,deferredInitialization);
    
  return statusMapMaker.StatusMap();
}
