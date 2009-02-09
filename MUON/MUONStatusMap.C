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
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVStore.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpIntPair.h"
#include "AliMpManuIterator.h"
#include "Riostream.h"
#endif

void FindBad(AliMUONPadStatusMaker& statusMaker, Int_t mask, Int_t& nBadPads, Int_t& nPads)
{
  AliMpManuIterator it;
  
  nBadPads = nPads = 0;
  
  Int_t detElemId;
  Int_t manuId;
  
  while ( it.Next(detElemId,manuId) )
  {
    Bool_t bad(kFALSE);
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    Int_t nb(0);
    Int_t n(0);

    for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i )
    {
      if ( de->IsConnectedChannel(manuId,i) )
      {
        ++n;
        ++nPads;
        Int_t status = statusMaker.PadStatus(detElemId,manuId,i);
        if ( ( status & mask) || (!mask && status) )
        {
          bad = kTRUE;
          ++nBadPads;
          ++nb;
        }
      }
    }
    
    if (bad)
    {
      cout << Form("DE %4d ManuId %4d %2d bad pads over %2d pads",
                   detElemId,manuId,nb,n) << endl;
    }
  }
}

AliMUONVStore* MUONStatusMap(const TString& cdbStorage = "local://$ALICE_ROOT/OCDB",
                             Int_t runNumber=0, Bool_t statusOnly=kFALSE, 
                             Int_t mask=0x8080)
{  
  AliCDBManager::Instance()->SetDefaultStorage(cdbStorage.Data());
  AliCDBManager::Instance()->SetRun(runNumber);

  AliMpCDB::LoadDDLStore();
  
  AliMUONCalibrationData cd(runNumber);
  
  AliMUONPadStatusMaker statusMaker(cd);
  
//  statusMaker.SetPedMeanLimits(50,200);
//  statusMaker.SetPedSigmaLimits(0.5,2);
  
  Int_t nbad;
  Int_t ntotal;
  
  FindBad(statusMaker,mask,nbad,ntotal);

  if (ntotal<=0) 
  {
    cout << "Error : got no pad at all ?!" << endl;
    return 0x0;
  }  
  
  cout << Form("Nbad = %6d over %6d pads (%7.2f %%)",
               nbad,ntotal,100.0*nbad/ntotal) << endl;
  
  if ( statusOnly ) return statusMaker.StatusStore();
  
  const Bool_t deferredInitialization = kFALSE;
  
  AliMUONPadStatusMapMaker statusMapMaker(cd,mask,deferredInitialization);
    
  return statusMapMaker.StatusMap();
}
