//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************
#include "AliRICHParam.h"
#include "AliRICHChamber.h"

ClassImp(AliRICHParam)
Bool_t   AliRICHParam::fgIsWireSag            =kTRUE;
Bool_t   AliRICHParam::fgIsResolveClusters    =kTRUE;
Bool_t   AliRICHParam::fgIsRadioSrc           =kFALSE;
Double_t AliRICHParam::fgAngleRot             =-60;
Int_t    AliRICHParam::fgHV[kNsectors]        ={2050,2050,2050,2050,2050,2050};
Int_t    AliRICHParam::fgNsigmaTh             =4;
Float_t  AliRICHParam::fgSigmaThMean          =1.132; //QDC 
Float_t  AliRICHParam::fgSigmaThSpread        =0.035; //     

//__________________________________________________________________________________________________
void AliRICHParam::Print(Option_t*)
{
  ::Info("","Pads in chamber (%3i,%3i) in sector (%2i,%2i)",NpadsX(),NpadsY(),NpadsXsec(),NpadsYsec());
  fpChambers->Print();
}
//__________________________________________________________________________________________________
void AliRICHParam::CreateChambers()
{
//Create all RICH Chambers on each call. Previous chambers deleted.
  if(fpChambers) delete fpChambers;
  fpChambers=new TObjArray(kNchambers);
  fpChambers->SetOwner();
  for(int iChamberN=0;iChamberN<kNchambers;iChamberN++)  fpChambers->AddAt(new AliRICHChamber(iChamberN+1),iChamberN);  
}//void AliRICH::CreateChambers()
