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
#include "AliT0LookUpValue.h"

ClassImp(AliT0LookUpValue)

  AliT0LookUpValue:: AliT0LookUpValue(): TObject(),
				     fTRM(-1),
				     fTDC(-1),
				     fChain(-1),
				     fChannel(-1)


{

  ///
}

 AliT0LookUpValue:: AliT0LookUpValue(Int_t trm, Int_t tdc, Int_t chain, Int_t channel) : TObject(),
  fTRM(trm),
  fTDC(tdc),
  fChain(chain),
  fChannel(channel)
{
 
}
 Bool_t AliT0LookUpValue:: IsEqual(const TObject* obj) const
{
     AliT0LookUpValue *b = ( AliT0LookUpValue*)obj;
/*
printf("IsEqual %i %i\n %i %i\n %i %i\n %i %i\n \n",b->GetTRM(), fTRM, 
b->GetTDC(), fTDC, b->GetChain(), fChain, b->GetChannel(), fChannel );
*/
     if ( b->GetTRM() == fTRM && 
	  b->GetTDC() == fTDC && 
	  b->GetChain() == fChain &&
	  b->GetChannel() == fChannel  ) return kTRUE;
     else return kFALSE;
}
ULong_t AliT0LookUpValue:: Hash(void) const 
{
  ULong_t hash = 1000*fTRM+100*fTDC+10*fChain+1*fChannel;
  //printf("hash=%i\n",hash);
  return hash; 
}
/*
void AliT0LookUpValue:: Clear()  
{
  fTRM=-1;
  fTDC=-1;
  fChain=-1;
  fChannel=-1;
}
*/

ClassImp(AliT0LookUpKey)

 AliT0LookUpKey::AliT0LookUpKey():TObject(),
 fKey(0)  
 {
 
 }

 AliT0LookUpKey::AliT0LookUpKey(Int_t key):TObject(),
 fKey(key)  
 {
 
 }



  Bool_t AliT0LookUpKey:: IsEqual(const TObject* obj) const
{
  Int_t k;
//  printf("IsEqual\n");
    k=((AliT0LookUpKey*)obj)->GetKey();
    if (k==fKey) return kTRUE; else return kFALSE;
}
