 
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSBaseline.h"
#include "AliHLTPHOSConstants.h"

ClassImp(AliHLTPHOSBaseline);
         
AliHLTPHOSBaseline::AliHLTPHOSBaseline() :
  TObject(),
  fBaseline(-1),
  fX(-1),
  fZ(-1),
  fGain(-1),
  fEntries(0)
{
}

AliHLTPHOSBaseline::~AliHLTPHOSBaseline()
{
}


