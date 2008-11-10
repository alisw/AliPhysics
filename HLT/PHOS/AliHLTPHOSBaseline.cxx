// $Id$

 
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

/**
 * Baseline class for PHOS HLT
 *
 * @file   AliHLTPHOSBaseline.cxx
 * @author Oystein Djuvsland
 * @date
 * @brief  Baseline class for PHOS HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

ClassImp(AliHLTPHOSBaseline);
         
AliHLTPHOSBaseline::AliHLTPHOSBaseline() :
  TObject(),
  fBaseline(-1),
  fX(-1),
  fZ(-1),
  fGain(-1),
  fEntries(0)
{
  //See header file for documentation
}

AliHLTPHOSBaseline::~AliHLTPHOSBaseline()
{
  //See header file for documentation
}
