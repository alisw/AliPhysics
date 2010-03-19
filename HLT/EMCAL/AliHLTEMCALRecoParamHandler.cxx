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

#include "AliHLTEMCALRecoParamHandler.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

AliHLTEMCALRecoParamHandler::AliHLTEMCALRecoParamHandler() :
AliHLTCaloRecoParamHandler("EMCAL")
{
   // See header file for class documentation
   
   fEMCALPidPtr = new AliEMCALPIDv1();
   
}


AliHLTEMCALRecoParamHandler::~AliHLTEMCALRecoParamHandler()
{
   // See header file for class documentation
   if(fRecoParamPtr) delete fRecoParamPtr; fRecoParamPtr = 0;
   
}

Int_t AliHLTEMCALRecoParamHandler::GetParametersFromCDB()
{
   // See header file for documentation
   AliCDBPath path("EMCAL","Calib","RecoParam");
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
	    fRecoParamPtr = dynamic_cast<AliEMCALRecoParam*>(paramArray)->At(0);
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

Float_t AliHLTEMCALRecoParamHandler::GetCorrectedEnergy ( Float_t e )
{
   // See header file for class documentation
   return e;
}

