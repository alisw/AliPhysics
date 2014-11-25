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

#include "AliHLTPHOSDigitHandler.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliHLTPHOSGeometry.h"
#include "TTree.h"
#include "AliPHOSDigit.h"

AliHLTPHOSDigitHandler *AliHLTPHOSDigitHandler::fgkInstance = NULL;

AliHLTPHOSDigitHandler::AliHLTPHOSDigitHandler() : AliHLTCaloDigitHandler("PHOS")
{

}

AliHLTPHOSDigitHandler::~AliHLTPHOSDigitHandler()
{

}

AliHLTPHOSDigitHandler* AliHLTPHOSDigitHandler::Instance()
{
    if (!fgkInstance)
    {
        fgkInstance = new AliHLTPHOSDigitHandler;
    }
    return fgkInstance;
}

Int_t AliHLTPHOSDigitHandler::Init(AliRunLoader* runLoader)
{
  
    fGeometry = new AliHLTPHOSGeometry();
    if(fGeometry) fGeometry->InitialiseGeometry();
    
    Int_t nev = AliHLTCaloDigitHandler::Init(runLoader);
    if(nev > 0)
    {
      if (fRunLoader)
      {
        fDetLoader = fRunLoader->GetDetectorLoader("PHOS");
        if (!fDetLoader)
        {
            HLTFatal("Could not get PHOS loader");
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

Int_t AliHLTPHOSDigitHandler::ConvertDigit(AliDigitNew* digit)
{
  AliPHOSDigit *dig = dynamic_cast<AliPHOSDigit*>(digit);
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
    hDig->fEnergy = dig->GetEnergy();
    hDig->fTime = dig->GetTime();
    hDig->fAmplitude = 0;
    hDig->fOverflow = false;
    hDig->fHgPresent = true;
    hDig->fAssociatedCluster = -1;
    fDigitsInModule[module]++;
    
    return 0;
}
