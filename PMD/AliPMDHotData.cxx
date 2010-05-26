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
#include "AliPMDHotData.h"


ClassImp(AliPMDHotData)

AliPMDHotData::AliPMDHotData()
{
  // Default constructor
  Reset();
}
// ----------------------------------------------------------------- //
AliPMDHotData::AliPMDHotData(const char* name)
{
  //constructor
  TString namst = "hot_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  
}
// ----------------------------------------------------------------- //
AliPMDHotData::AliPMDHotData(const AliPMDHotData& hotda) :
  TNamed(hotda)
{
  // copy constructor
  SetName(hotda.GetName());
  SetTitle(hotda.GetName());
  Reset();
  for(Int_t det = 0; det < kDet; det++)
    {
      for(Int_t smn = 0; smn < kModule; smn++)
	{
	  for(Int_t row = 0; row < kRow; row++)
	    {
	      for(Int_t col = 0; col < kCol; col++)
		{
		  fHotChannel[det][smn][row][col] = hotda.GetHotChannel(det,smn,row,col);
		}
	    }
	}
    }
}
// ----------------------------------------------------------------- //
AliPMDHotData &AliPMDHotData::operator =(const AliPMDHotData& hotda)
{
  //asignment operator
  SetName(hotda.GetName());
  SetTitle(hotda.GetName());
  Reset();
  for(Int_t det = 0; det < kDet; det++)
    {
      for(Int_t smn = 0; smn < kModule; smn++)
	{
	  for(Int_t row = 0; row < kRow; row++)
	    {
	      for(Int_t col = 0; col < kCol; col++)
		{
		  fHotChannel[det][smn][row][col] = hotda.GetHotChannel(det,smn,row,col);
		}
	    }
	}
    }
  return *this;
}
// ----------------------------------------------------------------- //
AliPMDHotData::~AliPMDHotData()
{
  //destructor
}
// ----------------------------------------------------------------- //
void AliPMDHotData::Reset()
{

  for(Int_t det = 0; det < kDet; det++)
    {
      for(Int_t smn = 0; smn < kModule; smn++)
	{
	  for(Int_t row = 0; row < kRow; row++)
	    {
	    for(Int_t col = 0; col < kCol; col++)
	      {
		fHotChannel[det][smn][row][col] = 0.;
	      }
	    }
	}
    }
}
// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //
Float_t AliPMDHotData:: GetHotChannel(Int_t det, Int_t smn, Int_t row, Int_t col) const
{
  return fHotChannel[det][smn][row][col];
}
void AliPMDHotData::SetHotChannel(Int_t det, Int_t smn, Int_t row, Int_t col, Float_t flag)
{
  fHotChannel[det][smn][row][col] = flag;
}
//------------------------------------------------------------------------------ //
void AliPMDHotData::Print(Option_t *) const
{
  printf("\n ######Flag for each cells ####\n");
  for(Int_t det = 0; det < kDet; det++)
    {
      for(Int_t smn = 0; smn < kModule; smn++)
	{
	  for(Int_t row = 0; row < kRow; row++)
	    {
	      for(Int_t col = 0; col < kCol; col++)
		{
		  printf("Flag[%d,%d,%d,%d]= %4.1f \n",det,smn,row,col,
			 fHotChannel[det][smn][row][col]);
		}
	      printf("\n");
	    }
	}
    }
}
