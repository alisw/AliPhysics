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
// Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONObjectPair.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONV2DStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDataIterator.h"
#include "AliMpIntPair.h"
#include "Riostream.h"
#endif

void findBad(const AliMUONV2DStore& status)
{
  AliMUONVDataIterator* it = status.Iterator();
  AliMUONObjectPair* pair;
  
  while ( ( pair = static_cast<AliMUONObjectPair*>(it->Next()) ) )
  {
    AliMpIntPair* p = static_cast<AliMpIntPair*>(pair->First());
    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();
    AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(pair->Second());
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

AliMUONV2DStore* MUONStatusMap(Int_t runNumber=0, Bool_t statusOnly=kFALSE, Int_t mask=0)
{
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB");
  
  AliMUONCalibrationData cd(runNumber);
  
  AliMUONPadStatusMaker statusMaker(cd);
  
//  statusMaker.SetPedMeanLimits(50,200);
  statusMaker.SetPedSigmaLimits(0.5,2);
  
  AliMUONV2DStore* status = statusMaker.MakeStatus();
 
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
  
  AliMUONPadStatusMapMaker statusMapMaker;
  
  return statusMapMaker.MakePadStatusMap(*status,mask);
}
