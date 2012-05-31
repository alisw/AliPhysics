/***************************************************************************
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
//
// Author : Z. Ahmed
//
#include "TNamed.h"
#include "AliCDBEntry.h"
#include "AliPMDMeanSm.h"


ClassImp(AliPMDMeanSm)

AliPMDMeanSm::AliPMDMeanSm()
{
  // Default constructor
  Reset();
}
// ----------------------------------------------------------------- //
AliPMDMeanSm::AliPMDMeanSm(const char* name)
{
  //constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  
}
// ----------------------------------------------------------------- //
AliPMDMeanSm::AliPMDMeanSm(const AliPMDMeanSm& meanda) :
  TNamed(meanda)
{
  // copy constructor
  SetName(meanda.GetName());
  SetTitle(meanda.GetName());
  Reset();
  for(Int_t det = 0; det < kDet; det++)
    {
      for(Int_t smn = 0; smn < kModule; smn++)
	{
	  fMeanSm[det][smn] = meanda.GetMeanSm(det,smn);
	}
    }
}
// ----------------------------------------------------------------- //
AliPMDMeanSm &AliPMDMeanSm::operator =(const AliPMDMeanSm& meanda)
{
  //asignment operator
  SetName(meanda.GetName());
  SetTitle(meanda.GetName());
  Reset();
  for(Int_t det = 0; det < kDet; det++)
    {
      for(Int_t smn = 0; smn < kModule; smn++)
	{
	  fMeanSm[det][smn] = meanda.GetMeanSm(det,smn);
	}
    }
  return *this;
}
// ----------------------------------------------------------------- //
AliPMDMeanSm::~AliPMDMeanSm()
{
  //destructor
}
// ----------------------------------------------------------------- //
void AliPMDMeanSm::Reset()
{
  //memset(fgainfact ,1,2*24*48*96*sizeof(Float_t));
  
  for(Int_t det = 0; det < kDet; det++)
    {
      for(Int_t smn = 0; smn < kModule; smn++)
	{
	  fMeanSm[det][smn] = 1.;
	}
    }
}
// ----------------------------------------------------------------- //
Float_t AliPMDMeanSm:: GetMeanSm(Int_t det, Int_t smn) const
{
  return fMeanSm[det][smn];
}
// ----------------------------------------------------------------- //
void AliPMDMeanSm::SetMeanSm(Int_t det, Int_t smn, Float_t smmean)
{
  fMeanSm[det][smn]= smmean;
}

// ----------------------------------------------------------------- //
void AliPMDMeanSm::Print(Option_t *) const
{
  printf("\n ######Mean Values of Super Modules are ####\n");
  for(Int_t det = 0; det < kDet; det++)
    {
      for(Int_t smn = 0; smn < kModule; smn++)
	{
	  {
	    printf("Gain[%d,%d]= %2.1f \n",det,smn,fMeanSm[det][smn]);
	  }
	  printf("\n");
	}
    }
}
