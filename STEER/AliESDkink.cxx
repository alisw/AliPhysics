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

/* $Id$ */

//-------------------------------------------------------------------------
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include <Riostream.h>
#include <TMath.h>
#include <TPDGCode.h>
#include "AliESDkink.h"


ClassImp(AliESDkink)

//____________________________________________________________________
AliESDkink::AliESDkink() :
  TObject(),
  fID(0),
  fParamDaughter(),
  fParamMother(),
  fDist1(-1),
  fDist2(-1),
  fRr(-1),
  fShapeFactor(0),
  fRow0(-1)
{
  //
  //Dafault constructor
  //
  for (Int_t i=0;i<12;i++) fStatus[i]=0;
  for (Int_t i=0;i<2;i++)
    for (Int_t j=0;j<2;j++){
      fTPCdensity[i][j]=-1;
      fTPCdensity2[i][j]=-1;
    }
  fTPCncls[0]=fTPCncls[1]=0;

  for (Int_t i=0; i<3; i++) {
    fPdr[i] = 0;
    fXr[i] = 0;
    fPm[i] = 0;
    fAngle[i] = 0;
  }
  
  fLab[0]=fLab[1]=-1;
  fIndex[0]=fIndex[1]=-1;
  fMultiple[0]=fMultiple[1]=0;
}

void AliESDkink::SetMother(const AliExternalTrackParam & pmother)  {
  //
  // set mother
  //
  fParamMother   = pmother;
}

void AliESDkink::SetDaughter(const AliExternalTrackParam & pdaughter){
  //
  //set daughter
  //
  fParamDaughter = pdaughter;

}
  
Float_t AliESDkink::GetTPCDensityFactor() const
{
  //
  //
  return fTPCdensity[0][0]+fTPCdensity[1][1]-TMath::Max(fTPCdensity[0][1],Float_t(0.0))-TMath::Max(fTPCdensity[1][0],Float_t(0.0)); 
}

Float_t AliESDkink::GetQt() const
{
  Float_t dmomentum = TMath::Sqrt(fPdr[0]*fPdr[0]+fPdr[1]*fPdr[1]+fPdr[2]*fPdr[2]);
  return TMath::Sin(fAngle[2])*dmomentum;
}
