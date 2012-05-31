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

#include "AliHLTCaloRecoParamHandler.h"
#include "AliDetectorRecoParam.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TObjArray.h"

ClassImp(AliHLTCaloRecoParamHandler);

AliHLTCaloRecoParamHandler::AliHLTCaloRecoParamHandler ( TString det ) : 
AliHLTCaloConstantsHandler(det)
,AliHLTLogging()
,fLogWeight(4.5)
,fRecPointMemberThreshold(0.01)
,fRecPointThreshold(0.1)
,fRecoParamPtr(0)
,fRecoParamPath(det, "Calib", "RecoParam")
{
   // See header file for class documentation

}


AliHLTCaloRecoParamHandler::~AliHLTCaloRecoParamHandler()
{
   // See header file for class documentation
   
}



Int_t AliHLTCaloRecoParamHandler::GetParametersFromCDB()
{
   // See header file for documentation

   if(fRecoParamPath.GetPath())
    {
//      HLTInfo("configure from entry %s", path.GetPath());
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(fRecoParamPath/*,GetRunNo()*/);
      if (pEntry) 
	{
	    
	    TObjArray *paramArray = dynamic_cast<TObjArray*>(pEntry->GetObject());
	    if(paramArray)
	      {
		fRecoParamPtr = dynamic_cast<AliDetectorRecoParam*>((paramArray)->At(0));
	      }
	    if(!fRecoParamPtr)
	      {
		HLTError("can not fetch object reconstruction parameters from \"%s\"", fRecoParamPath.GetPath().Data());
		return -1;
	      }
	}
      else
	{
	    HLTError("can not fetch object \"%s\" from OCDB", fRecoParamPath.GetPath().Data());
	    return -1;
	}
    }
    FillParameters();
    return 0;
}
