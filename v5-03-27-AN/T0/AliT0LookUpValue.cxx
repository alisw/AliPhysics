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
//#include "iostream.h"

ClassImp(AliT0LookUpValue)

  AliT0LookUpValue:: AliT0LookUpValue(): TObject(),
				     fTRM(-1),
				     fTDC(-1),
				     fChain(-1),
				     fChannel(-1)


{

  ///
}
//--------------------------------------------------------------------------------------
 AliT0LookUpValue:: AliT0LookUpValue(Int_t trm, Int_t tdc, Int_t chain, Int_t channel) : TObject(),
  fTRM(trm),
  fTDC(tdc),
  fChain(chain),
  fChannel(channel)
{

//--------------------------------------------------------------------------------------
 
}
 Bool_t AliT0LookUpValue:: IsEqual(const TObject* obj) const
{
     AliT0LookUpValue *b = ( AliT0LookUpValue*)obj;
     if ( b->GetTRM() == fTRM && 
	  b->GetTDC() == fTDC && 
	  b->GetChain() == fChain &&
	  b->GetChannel() == fChannel  ) return kTRUE;
     else return kFALSE;
}

//------------------------------------------------------------------------------------ 
void AliT0LookUpValue::Print(Option_t *) const
{
  printf("AliT0LookUpValue TRM %i, TDC %i, chain %i channel %i\n",fTRM,fTDC,fChain,fChannel);
}
//--------------------------------------------------------------------------------------
