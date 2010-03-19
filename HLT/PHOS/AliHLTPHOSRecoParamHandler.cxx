/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSRecoParamHandler.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliPHOSRecoParam.h"
#include "AliPHOSPIDv1.h"
#include "TObjArray.h"

AliHLTPHOSRecoParamHandler::AliHLTPHOSRecoParamHandler() :
AliHLTCaloRecoParamHandler("PHOS")
,fRecoParamPtr(0),
fPHOSPidPtr(0)
{
   // See header file for class documentation
   
   fPHOSPidPtr = new AliPHOSPIDv1();
   
}


AliHLTPHOSRecoParamHandler::~AliHLTPHOSRecoParamHandler()
{
   // See header file for class documentation
   if(fRecoParamPtr) delete fRecoParamPtr; fRecoParamPtr = 0;
   
}

Int_t AliHLTPHOSRecoParamHandler::GetParametersFromCDB()
{
   // See header file for documentation
   AliCDBPath path("PHOS","Calib","RecoParam");
   if(path.GetPath())
    {
//      HLTInfo("configure from entry %s", path.GetPath());
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
      if (pEntry) 
	{
	  if(!fRecoParamPtr) 
	    {
	      delete fRecoParamPtr;
	      fRecoParamPtr = 0;
	    }
	    TObjArray *paramArray = dynamic_cast<TObjArray*>(pEntry->GetObject());
	    fRecoParamPtr = dynamic_cast<AliPHOSRecoParam*>(paramArray->At(0));
	    if(!fRecoParamPtr)
	    {
	       return -1;
	    }
	    fLogWeight = fRecoParamPtr->GetEMCLogWeight();
	    fRecPointMemberThreshold = fRecoParamPtr->GetEMCMinE();
	    fRecPointThreshold = fRecoParamPtr->GetEMCClusteringThreshold();
	}
      else
	{
//	    HLTError("can not fetch object \"%s\" from OCDB", path);
	    return -1;
	}
    }
    return 0;
}

Float_t AliHLTPHOSRecoParamHandler::GetCorrectedEnergy ( Float_t e )
{
   // See header file for class documentation
   return fPHOSPidPtr->GetCalibratedEnergy(e);
}

