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

AliHLTCaloRecoParamHandler::AliHLTCaloRecoParamHandler ( TString det ) : 
AliHLTCaloConstantsHandler(det)
,fLogWeight(4.5)
,fRecPointMemberThreshold(0.01)
,fRecPointThreshold(0.1)
{
   // See header file for class documentation
   
}


AliHLTCaloRecoParamHandler::~AliHLTCaloRecoParamHandler()
{
   // See header file for class documentation
   
}

