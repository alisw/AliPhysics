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
//
#include "TNamed.h"
#include "AliCDBEntry.h"
#include "AliPMD.h"
#include "AliPMDPedestal.h"


ClassImp(AliPMDPedestal)

AliPMDPedestal::AliPMDPedestal()
{
  // Default constructor
  Reset();
}
// ----------------------------------------------------------------- //
AliPMDPedestal::AliPMDPedestal(const char* name)
{
  //constructor
  TString namst = "Pedestal_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  
}
// ----------------------------------------------------------------- //
AliPMDPedestal::AliPMDPedestal(const AliPMDPedestal &pedestal) :
  TNamed(pedestal)
{
  // copy constructor
  SetName(pedestal.GetName());
  SetTitle(pedestal.GetName());
  Reset();
  for(Int_t det = 0; det < kDet; det++)
  {
      for(Int_t smn = 0; smn < kModule; smn++)
      {
	  for(Int_t row = 0; row < kRow; row++)
	  {
	      for(Int_t col = 0; col < kCol; col++)
	      {
		  fPedMeanRms[det][smn][row][col] =
		      pedestal.GetPedMeanRms(det,smn,row,col);
	      }
	  }
      }
  }
}
// ----------------------------------------------------------------- //
AliPMDPedestal &AliPMDPedestal::operator =(const AliPMDPedestal &pedestal)
{
  //asignment operator
  SetName(pedestal.GetName());
  SetTitle(pedestal.GetName());
  Reset();
  for(Int_t det = 0; det < kDet; det++)
  {
      for(Int_t smn = 0; smn < kModule; smn++)
      {
	  for(Int_t row = 0; row < kRow; row++)
	  {
	      for(Int_t col = 0; col < kCol; col++)
	      {
		  fPedMeanRms[det][smn][row][col] = 
		      pedestal.GetPedMeanRms(det,smn,row,col);
	      }
	  }
      }
  }
  return *this;
}
// ----------------------------------------------------------------- //
AliPMDPedestal::~AliPMDPedestal()
{
  //destructor
}
// ----------------------------------------------------------------- //
void AliPMDPedestal::Reset()
{
  //memset(fgainfact ,1,2*24*96*96*sizeof(Float_t));
    Float_t mean = 0.0;
    Float_t rms  = 0.0;
    for(Int_t det = 0; det < kDet; det++)
    {
	for(Int_t smn = 0; smn < kModule; smn++)
	{
	    for(Int_t row = 0; row < kRow; row++)
	    {
		for(Int_t col = 0; col < kCol; col++)
		{
		    Int_t mean1 = (Int_t) (mean*10.);
		    Int_t rms1  = (Int_t) (rms*10.);
		    fPedMeanRms[det][smn][row][col] = mean1*100 + rms1;
		}
	    }
	}
    }
}
// ----------------------------------------------------------------- //
Int_t AliPMDPedestal:: GetPedMeanRms(Int_t det, Int_t smn, Int_t row, Int_t col) const
{
  return fPedMeanRms[det][smn][row][col];
}
// ----------------------------------------------------------------- //

void AliPMDPedestal::SetPedMeanRms(Int_t det, Int_t smn, Int_t row,
				   Int_t col,Float_t pedmean, Float_t pedrms)
{
    Int_t mean1 = (Int_t) (pedmean*10.);
    Int_t rms1  = (Int_t) (pedrms*10.); 
    fPedMeanRms[det][smn][row][col] = mean1*100 + rms1;
}
// ----------------------------------------------------------------- //
void AliPMDPedestal::Print(Option_t *) const
{
  printf("\n ###### Pedestal values for each cells are ####\n");
  for(Int_t det = 0; det < kDet; det++)
    {
      for(Int_t smn = 0; smn < kModule; smn++)
	{
	  for(Int_t row = 0; row < kRow; row++)
	    {
	      for(Int_t col = 0; col < kCol; col++)
		{
		    printf("Pedestal[%d,%d,%d,%d]= %d \n",det,smn,row,col,
			   fPedMeanRms[det][smn][row][col]);

		}
	      printf("\n");
	    }
	}
    }
}
// ----------------------------------------------------------------- //
