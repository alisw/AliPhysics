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
#include "AliPMD.h"
#include "AliPMDCalibData.h"


ClassImp(AliPMDCalibData)

AliPMDCalibData::AliPMDCalibData()
{
  // Default constructor
  Reset();
}
// ----------------------------------------------------------------- //
AliPMDCalibData::AliPMDCalibData(const char* name)
{
  //constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  
}
// ----------------------------------------------------------------- //
AliPMDCalibData::AliPMDCalibData(const AliPMDCalibData& calibda) :
  TNamed(calibda)
{
  // copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(Int_t det=0;det<2;det++){
    for(Int_t smn=0;smn<24;smn++) {
      for(Int_t row=0;row<96;row++) {
	for(Int_t col=0;col<96;col++) {
	  fGainFact[det][smn][row][col]=calibda.GetGainFact(det,smn,row,col);
	}
      }
    }
  }
}
// ----------------------------------------------------------------- //
AliPMDCalibData &AliPMDCalibData::operator =(const AliPMDCalibData& calibda)
{
  //asignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(Int_t det=0;det<2;det++){
    for(Int_t smn=0;smn<24;smn++) {
      for(Int_t row=0;row<96;row++) {
	for(Int_t col=0;col<96;col++) {
	  fGainFact[det][smn][row][col]=calibda.GetGainFact(det,smn,row,col);
	}
      }
    }
  }
  return *this;
}
// ----------------------------------------------------------------- //
AliPMDCalibData::~AliPMDCalibData()
{
  //destructor
}
// ----------------------------------------------------------------- //
void AliPMDCalibData::Reset()
{
  //memset(fgainfact ,1,2*24*96*96*sizeof(Float_t));
  for(Int_t det=0;det<2;det++){
    for(Int_t smn=0;smn<24;smn++) {
      for(Int_t row=0;row<96;row++) {
	for(Int_t col=0;col<96;col++) {
	  fGainFact[det][smn][row][col]=1.0;
	}
      }
    }
  }
}
// ----------------------------------------------------------------- //
Float_t AliPMDCalibData:: GetGainFact(Int_t det, Int_t smn, Int_t row, Int_t col) const
{
  return fGainFact[det][smn][row][col];
}
// ----------------------------------------------------------------- //
void AliPMDCalibData::SetGainFact(Int_t det, Int_t smn, Int_t row, Int_t col, Float_t gain)
{
  fGainFact[det][smn][row][col]= gain;
}

// ----------------------------------------------------------------- //
void AliPMDCalibData::Print(Option_t *) const
{
  printf("\n ######gain factors for each cells are ####\n");
  for(Int_t det=0;det<2;det++)
    {
      for(Int_t smn=0;smn<24;smn++)
	{
	  for(Int_t row=0;row<96;row++)
	    {
	      for(Int_t col=0;col<96; col++)
		{
		  printf("%4.1f",fGainFact[det][smn][row][col]);
		}
	      printf("\n");
	    }
	}
    }
}
