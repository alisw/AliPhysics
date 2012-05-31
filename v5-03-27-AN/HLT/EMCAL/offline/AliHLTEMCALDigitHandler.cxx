/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland <oysteind@ift.uib.no>               *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTEMCALDigitHandler.h"
#include "AliRunLoader.h"
#include "AliEMCALLoader.h"
#include "AliHLTEMCALGeometry.h"
#include "TTree.h"
#include "AliEMCALDigit.h"
#include "AliEMCALCalibData.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"

AliHLTEMCALDigitHandler *AliHLTEMCALDigitHandler::fgkInstance = NULL;

AliHLTEMCALDigitHandler::AliHLTEMCALDigitHandler() : AliHLTCaloDigitHandler("EMCAL")
,fCalibData(0)
{

}

AliHLTEMCALDigitHandler::~AliHLTEMCALDigitHandler()
{

}

AliHLTEMCALDigitHandler* AliHLTEMCALDigitHandler::Instance()
{
    if (!fgkInstance)
    {
        fgkInstance = new AliHLTEMCALDigitHandler;
    }
    return fgkInstance;
}

Int_t AliHLTEMCALDigitHandler::Init(AliRunLoader* runLoader)
{

    fGeometry = new AliHLTEMCALGeometry();
    if (fGeometry) fGeometry->InitialiseGeometry();
    if(GetGainsFromCDB())
    {
      HLTFatal("Could not get gains from CDB");
      return -3;
    }

    Int_t nev = AliHLTCaloDigitHandler::Init(runLoader);
    if (nev > 0)
    {
        if (fRunLoader)
        {
            fDetLoader = dynamic_cast<AliEMCALLoader*>(fRunLoader->GetDetectorLoader("EMCAL"));
            if (!fDetLoader)
            {
                HLTFatal("Could not get EMCAL loader");
                return -1;
            }
        }
        else
        {
            return -2;
        }
    }

    return nev;
}

Int_t AliHLTEMCALDigitHandler::ConvertDigit(AliDigitNew *digit)
{
    AliEMCALDigit *dig = dynamic_cast<AliEMCALDigit*>(digit);

    if(!dig)
    {
      HLTError("Wrong data, cannot create digits");
      return -1;
    }
    
    Int_t module = 0;
    Int_t x = 0;
    Int_t z = 0;

    fGeometry->GetLocalCoordinatesFromAbsId(dig->GetId(), module, x, z);
    AliHLTCaloDigitDataStruct *hDig = &(fDigits[module][fDigitsInModule[module]]);

    hDig->fID = dig->GetId();
    hDig->fX = x;
    hDig->fZ = z;
    hDig->fModule = module;
    hDig->fEnergy = dig->GetAmplitude()*(fCalibData->GetADCchannel(module, z, x));
    hDig->fTime = dig->GetTime();
    hDig->fAmplitude = 0;
    hDig->fOverflow = false;
    hDig->fHgPresent = true;
    hDig->fAssociatedCluster = -1;
    fDigitsInModule[module]++;
    
  return 0;
}

int AliHLTEMCALDigitHandler::GetGainsFromCDB()
{
   // See header file for class documentation

  AliCDBPath path("EMCAL","Calib","Data");
  if(path.GetPath())
    {
      //      HLTInfo("configure from entry %s", path.GetPath());
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
	if (pEntry) 
	{
	    fCalibData = (AliEMCALCalibData*)pEntry->GetObject();
	}
      else
	{
//	    HLTError("can not fetch object \"%s\" from CDB", path);
	    return -1;
	}
    }
   if(!fCalibData) return -1;
   return 0;
}
