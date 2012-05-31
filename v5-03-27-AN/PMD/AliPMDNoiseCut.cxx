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
// Author : 
//
#include "TNamed.h"
#include "AliCDBEntry.h"
#include "AliPMDNoiseCut.h"


ClassImp(AliPMDNoiseCut)

AliPMDNoiseCut::AliPMDNoiseCut()
{
  // Default constructor
  for(Int_t imod = 0; imod < 48; imod++)
    {
      fNoiseCut[imod] = 0.;
    }
}
// ----------------------------------------------------------------- //
AliPMDNoiseCut::AliPMDNoiseCut(const char* name)
{
  //constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  for(Int_t imod = 0; imod < 48; imod++)
    {
      fNoiseCut[imod] = 0.;
    }
}
// ----------------------------------------------------------------- //
AliPMDNoiseCut::AliPMDNoiseCut(const AliPMDNoiseCut& noisecut) :
  TNamed(noisecut)
{
  // copy constructor
  SetName(noisecut.GetName());
  SetTitle(noisecut.GetName());
  for(Int_t imod = 0; imod < 48; imod++)
    {
      fNoiseCut[imod] = noisecut.GetNoiseCut(imod);
    }

}
// ----------------------------------------------------------------- //
AliPMDNoiseCut &AliPMDNoiseCut::operator =(const AliPMDNoiseCut& noisecut)
{
  //asignment operator
  SetName(noisecut.GetName());
  SetTitle(noisecut.GetName());

  for(Int_t imod = 0; imod < 48; imod++)
    {
      fNoiseCut[imod] = noisecut.GetNoiseCut(imod);
    }

  return *this;
}
// ----------------------------------------------------------------- //
AliPMDNoiseCut::~AliPMDNoiseCut()
{
  //destructor
}
// ----------------------------------------------------------------- //

void AliPMDNoiseCut::Print(Option_t *) const
{
  printf("\n ######Noise Cut for each Module ####\n");


  for(Int_t imod = 0; imod < 48; imod++)
    {
      printf("%d %f\n",imod,fNoiseCut[imod]);
      printf("\n");
    }

}
