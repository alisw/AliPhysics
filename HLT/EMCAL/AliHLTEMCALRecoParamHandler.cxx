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
#include "AliEMCALRecParam.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

AliHLTEMCALRecoParamHandler::AliHLTEMCALRecoParamHandler() :
AliHLTCaloRecoParamHandler("EMCAL")
{
   // See header file for class documentation
   
}

AliHLTEMCALRecoParamHandler::~AliHLTEMCALRecoParamHandler()
{
   // See header file for class documentation
   
}

Float_t AliHLTEMCALRecoParamHandler::GetCorrectedEnergy ( Float_t e )
{
   // See header file for class documentation
   return e;
}


Int_t 
AliHLTEMCALRecoParamHandler::GetParametersFromCDB()
{
  // Avoiding linking error in EMCAL
  return 0;
}



void AliHLTEMCALRecoParamHandler::FillParameters()
{
   //See header file for class documentation
   if(fRecoParamPtr)
   {
      fLogWeight = dynamic_cast<AliEMCALRecParam*>(fRecoParamPtr)->GetW0(); 
      fRecPointMemberThreshold = dynamic_cast<AliEMCALRecParam*>(fRecoParamPtr)->GetMinECut();
      fRecPointThreshold = dynamic_cast<AliEMCALRecParam*>(fRecoParamPtr)->GetClusteringThreshold();
   }
}
