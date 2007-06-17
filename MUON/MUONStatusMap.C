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

/// Macro to check/test pad status and pad status map makers
///
// Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONVStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMpIntPair.h"
#include "Riostream.h"
#endif


void findBad(const AliMUONVStore& status)
{
  TIter next(status.CreateIterator());
  AliMUONVCalibParam* param;
  
  while ( ( param = dynamic_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t detElemId = param->ID0();
    Int_t manuId = param->ID1();
    Bool_t bad(kFALSE);
    for ( Int_t i = 0; i < param->Size(); ++i ) 
    {
      if ( param->ValueAsInt(0) ) bad = kTRUE;
    }
    if (bad)
    {
      cout << Form("DE %4d ManuId %4d",detElemId,manuId) << endl;
    }
  }
}

AliMUONVStore* MUONStatusMap(const TString& cdbStorage = "local://$ALICE_ROOT",
                               Int_t runNumber=0, Bool_t statusOnly=kFALSE, Int_t mask=0)
{  
  AliCDBManager::Instance()->SetDefaultStorage(cdbStorage.Data());

  AliMUONCalibrationData cd(runNumber);
  
  AliMUONPadStatusMaker statusMaker(cd);
  
//  statusMaker.SetPedMeanLimits(50,200);
  statusMaker.SetPedSigmaLimits(0.5,2);
  
  AliMUONVStore* status = statusMaker.MakeStatus();
 
  if ( status )
  {
    findBad(*status);
  }
  else
  {
    cout << "ERROR. Could not get status from CDB" << endl;
    return 0;
  }
  
  if ( statusOnly ) return status;
  
  AliMUONPadStatusMapMaker statusMapMaker(cd);
  
  return statusMapMaker.MakePadStatusMap(*status,mask);
}
