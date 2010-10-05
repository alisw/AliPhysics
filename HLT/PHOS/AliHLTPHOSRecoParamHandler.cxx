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
#include "TMatrixF.h"
#include "TVector3.h"
#include "AliPHOSReconstructor.h"
#include "TObjArray.h"

AliHLTPHOSRecoParamHandler::AliHLTPHOSRecoParamHandler() :
AliHLTCaloRecoParamHandler("PHOS")
{
   // See header file for class documentation


}

AliHLTPHOSRecoParamHandler::~AliHLTPHOSRecoParamHandler()
{
   // See header file for class documentation
}

Float_t AliHLTPHOSRecoParamHandler::GetCorrectedEnergy ( Float_t e )
{
   // See header file for class documentation
   return AliPHOSReconstructor::CorrectNonlinearity(e) ;
}

void AliHLTPHOSRecoParamHandler::FillParameters()
{
   //See header file for class documentation
   if(fRecoParamPtr)
   {
      fLogWeight = dynamic_cast<AliPHOSRecoParam*>(fRecoParamPtr)->GetEMCLogWeight();
      fRecPointMemberThreshold = dynamic_cast<AliPHOSRecoParam*>(fRecoParamPtr)->GetEMCMinE();
      fRecPointThreshold = dynamic_cast<AliPHOSRecoParam*>(fRecoParamPtr)->GetEMCClusteringThreshold();
   }
}

