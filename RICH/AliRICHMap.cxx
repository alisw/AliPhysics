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

#include "AliRICHMap.h"
#include "AliRICH.h"

ClassImp(AliRICHMap)

AliRICHMap::AliRICHMap(TClonesArray *pDig)
{
// main ctor
  fDigits=pDig;
  fNdigits=fDigits->GetEntries();
  fMap=new TMatrix(1,AliRICHParam::NpadsX(),1,AliRICHParam::NpadsY());
  Clear("");
  FillHits();
}
//__________________________________________________________________________________________________	
void  AliRICHMap::FillHits()
{
// Loops over the list of digits filling the "pad fired by digits" structure
  if(!fNdigits) return;    
  for(Int_t iDigN=0;iDigN<fNdigits;iDigN++){
    AliRICHDigit *pDig= (AliRICHDigit*)fDigits->At(iDigN);
    SetHit(pDig->PadX(),pDig->PadY(),iDigN);
  }
}
//__________________________________________________________________________________________________
