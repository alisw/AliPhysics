/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class to set HMPID reconstruction parameters (normal, HTA, UserCut ...    //
//                                                                           //  
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//
//Email: Levente.Molnar@ba.infn.it
//

#include "AliHMPIDRecoParam.h"
#include "AliHMPIDParam.h"

ClassImp(AliHMPIDRecoParam)

//_____________________________________________________________________________
AliHMPIDRecoParam::AliHMPIDRecoParam():TNamed(),
  fRecoMode(kTRUE),
  fUserCutMode(kTRUE)
{
  //
  // ctor
  //
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) fUserCut[iCh]=3;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecoParam::AliHMPIDRecoParam(const AliHMPIDRecoParam &p):TNamed(p),
    fRecoMode(kTRUE),
    fUserCutMode(kTRUE)
{ 
   //copy Ctor

   fRecoMode= p.fRecoMode;
   for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) fUserCut[iCh]=p.fUserCut[iCh];
   fUserCutMode = p.fUserCutMode;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecoParam& AliHMPIDRecoParam::operator=(const AliHMPIDRecoParam &p)
{
//
// assign. operator
//
  this->fRecoMode= p.fRecoMode;
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) this->fUserCut[iCh] = p.fUserCut[iCh];       
  this->fUserCutMode = p.fUserCutMode;       
  return *this;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecoParam *AliHMPIDRecoParam::GetUserModeParam(){
  //
  // Provide access to reconstruction parameters fro  the rec.C
  //
  AliHMPIDRecoParam *hmpidRecoParam = new AliHMPIDRecoParam;
  return hmpidRecoParam;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecoParam::~AliHMPIDRecoParam() 
{
  //
  // dtor
  //  
}
