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

#ifndef ALIHLTPHOSSANITYINSPECTOR_H
#define ALIHLTPHOSSANITYINSPECTOR_H

#include "AliHLTPHOSBase.h"
#include "AliHLTAltroData.h"
#include "Rtypes.h"

class AliHLTPHOSSanityInspector : public AliHLTPHOSBase
{
  
public:
  AliHLTPHOSSanityInspector();
  
  virtual ~AliHLTPHOSSanityInspector();
   
  Int_t CheckInsanity(UInt_t*, Int_t);
  Int_t CheckInsanity(Int_t*, Int_t);
  
  Int_t CheckAndHealInsanity(UInt_t*, Int_t);  //Not completely reliable
  Int_t CheckAndHealInsanity(Int_t*, Int_t);  //Not completely reliable

  void SetMaxDifference(Int_t maxDiff) { fMaxDifference = maxDiff; }
    
private:
  Int_t fMaxDifference;
 
  ClassDef(AliHLTPHOSSanityInspector, 1);
};

#endif
