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

#include "AliESDRun.h"
#include "AliESDVertex.h"

//-------------------------------------------------------------------------
//                     Implementation Class AliESDRun
//   Run by run data
//   for the ESD   
//   Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------

ClassImp(AliESDRun)  
 
//______________________________________________________________________________
AliESDRun::AliESDRun() :
  TObject(),
  fMagneticField(0),
  fPeriodNumber(0),
  fRunNumber(0),
  fRecoVersion(0) 
{
  for (Int_t i=0; i<2; i++) fDiamondXY[i]=0.;
  for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=0.;
}

//______________________________________________________________________________
AliESDRun::AliESDRun(const AliESDRun &esd) :
  TObject(esd),
  fMagneticField(esd.fMagneticField),
  fPeriodNumber(esd.fPeriodNumber),
  fRunNumber(esd.fRunNumber),
  fRecoVersion(esd.fRecoVersion)
{ 
  // Copy constructor
  for (Int_t i=0; i<2; i++) fDiamondXY[i]=esd.fDiamondXY[i];
  for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=esd.fDiamondCovXY[i];
}

//______________________________________________________________________________
AliESDRun& AliESDRun::operator=(const AliESDRun &esd)
{ 
  // assigment operator
  if(this!=&esd) {
    TObject::operator=(esd);
    fRunNumber=esd.fRunNumber;
    fPeriodNumber=esd.fPeriodNumber;
    fRecoVersion=esd.fRecoVersion;
    fMagneticField=esd.fMagneticField;
    for (Int_t i=0; i<2; i++) fDiamondXY[i]=esd.fDiamondXY[i];
    for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=esd.fDiamondCovXY[i];
  } 
  return *this;
}

void AliESDRun::SetDiamond(const AliESDVertex *vertex) {
  // set the interaction diamond
  fDiamondXY[0]=vertex->GetXv();
  fDiamondXY[1]=vertex->GetYv();
  Double32_t cov[6];
  vertex->GetCovMatrix(cov);
  fDiamondCovXY[0]=cov[0];
  fDiamondCovXY[1]=cov[1];
  fDiamondCovXY[2]=cov[2];
}


//______________________________________________________________________________
void AliESDRun::Print(const Option_t *) const
{
  // Print some data members
  printf("Mean vertex in RUN %d: X=%.4f Y=%.4f cm\n",
	 GetRunNumber(),GetDiamondX(),GetDiamondY());
  printf("Magnetic field = %f T\n",
	 GetMagneticField());
  printf("Event from reconstruction version %d \n",fRecoVersion);
}

void AliESDRun::Reset() 
{
  // reset data members
  fRunNumber = 0;
  fPeriodNumber = 0;
  fRecoVersion = 0;
  fMagneticField = 0;
  for (Int_t i=0; i<2; i++) fDiamondXY[i]=0.;
  for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=0.;
}

