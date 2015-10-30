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
#include "TArrayI.h"
#include "AliITSMFTCalibrationPix.h"
#include "AliLog.h"

///////////////////////////////////////////////////////////////////////////
//
//  Calibration class for set:ITS                   
//  Specific subdetector implementation for         
//  Silicon pixels                                  
//
//  Modified by D. Elia, G.E. Bruno, H. Tydesjo
//  Adapted for upgrade ruben.shahoyan@cern.ch
//
///////////////////////////////////////////////////////////////////////////

ClassImp(AliITSMFTCalibrationPix)

//______________________________________________________________________
AliITSMFTCalibrationPix::AliITSMFTCalibrationPix() 
: fDataType(), fNChips(0)
  ,fNColPerChip(0)
  ,fNCol(0)
  ,fNRow(0)
  ,fNrBadSingle(0)
  ,fBadChips(0)
  ,fBadChannels(0)
{
  // constructor
   SetDataType("simulated");
   ClearBad();
}

//______________________________________________________________________
AliITSMFTCalibrationPix::AliITSMFTCalibrationPix(Short_t nChips,Short_t nColPerChip,Short_t nRow) 
: fDataType(), fNChips(0)
  ,fNColPerChip(0)
  ,fNCol(0)
  ,fNRow(0)
  ,fNrBadSingle(0)
  ,fBadChips(0)
  ,fBadChannels(0)
{
  // constructor
   SetDataType("simulated");
   SetColRowData(nChips,nColPerChip,nRow);
   ClearBad();
}

//______________________________________________________________________
AliITSMFTCalibrationPix::AliITSMFTCalibrationPix(const AliITSMFTCalibrationPix &src) : 
  TObject(src),fDataType()
  ,fNChips(src.fNChips)
  ,fNColPerChip(src.fNColPerChip)
  ,fNCol(src.fNCol)
  ,fNRow(src.fNRow)
  ,fNrBadSingle(src.fNrBadSingle)
  ,fBadChips(src.fBadChips)
  ,fBadChannels(src.fBadChannels)
{
}

//____________________________________________________________________________
void AliITSMFTCalibrationPix::ClearBad() 
{
  // clear all bad pixels (single+chips)
  fBadChannels.Reset();
  fNrBadSingle=0;
  fBadChips = 0;
  //
}

//____________________________________________________________________________
void AliITSMFTCalibrationPix::AddBad(Int_t col, Int_t row) 
{
  // add single bad pixel 
  fBadChannels.Set(fNrBadSingle*2+2);
  fBadChannels.AddAt(col,fNrBadSingle*2);
  fBadChannels.AddAt(row,fNrBadSingle*2+1);
  fNrBadSingle++;
}

//____________________________________________________________________________
void AliITSMFTCalibrationPix::SetChipBad(Int_t chip) 
{
  // set full chip bad
  if ((int)chip>=fNChips) {AliError(Form("chip number %d exceeds allowed limit %d",chip,fNChips)); return;}
  fBadChips |= 0x1<<chip;
  //
}

//____________________________________________________________________________
void AliITSMFTCalibrationPix::UnSetChipBad(Int_t chip) 
{
  // unset full chip bad
  if (chip>=fNChips) {AliError(Form("chip number %d exceeds allowed limit %d",chip,fNChips)); return;}
  fBadChips &= ~(0x1<<chip);
  //
}

//____________________________________________________________________________
Int_t AliITSMFTCalibrationPix::GetBadColAt(Int_t index) const 
{
  // Get column of index-th bad pixel
  int nrc = fNColPerChip*fNRow;
  if (nrc<1) AliFatal("Number of colums and rows is not set");
  //
  if ((Int_t)index<GetNrBadSingle()) return fBadChannels.At(index*2);
  else {
    Int_t badChipIndex=(index-GetNrBadSingle())/nrc;
    Int_t badChipsFound =0;
    for (int chip=fNChips; chip--;) {
      if (IsChipMarkedBad(chip)) badChipsFound++;
      if (badChipIndex==badChipsFound-1) {
	Int_t badPixelIndex=(index-GetNrBadSingle())%(nrc);
	return chip*fNColPerChip + badPixelIndex/fNRow;
      }
    }
  }
  AliError(Form("Index %d is out of bounds - returning -1",index));
  return -1;
}

//____________________________________________________________________________
Int_t AliITSMFTCalibrationPix::GetBadRowAt(Int_t index) const 
{
  // Get row of index-th bad pixel
  int nrc = fNColPerChip*fNRow;
  if (nrc<1) AliFatal("Number of colums and rows is not set");
  //
  if ((Int_t)index<GetNrBadSingle()) return fBadChannels.At(index*2+1);
  else {
    Int_t badChipIndex=(index-GetNrBadSingle())/nrc;
    Int_t badChipsFound =0;
    for (int chip=fNChips; chip--;) {
      if (IsChipMarkedBad(chip)) badChipsFound++;
      if (badChipIndex==badChipsFound-1) {
	Int_t badPixelIndex=(index-GetNrBadSingle())%nrc;
	return badPixelIndex%fNRow;
      }
    }
  }
  AliError(Form("Index %d is out of bounds - returning -1",index));
  return -1;
}

//____________________________________________________________________________
void AliITSMFTCalibrationPix::GetBadPixel(Int_t index, Int_t &row, Int_t &col) const 
{
  // i: is the i-th bad pixel in single bad pixel list
  // row: is the corresponding row (-1 if i is out of range)
  // col: is the corresponding column (-1 if i is out of range)
  row = -1;
  col = -1;
  if(index>=0 && index<GetNrBadSingle()){
    col = GetBadColAt(index);
    row = GetBadRowAt(index);
    return;
  }
  else {
    if (index>=0) {
      int nrc = fNColPerChip*fNRow;
      if (nrc<1) AliFatal("Number of colums and rows is not set");
      Int_t badChipIndex=(index-GetNrBadSingle())/nrc;
      Int_t badChipsFound =0;
      for (int chip=fNChips; chip--;) {
	if (IsChipMarkedBad(chip)) badChipsFound++;
	if (badChipIndex==badChipsFound-1) {
	  Int_t badPixelIndex=(index-GetNrBadSingle())%nrc;
	  col = chip*fNColPerChip + badPixelIndex/fNRow;
	  row = badPixelIndex%fNRow;
	  return;
	}
      }
    }
  }
  AliError(Form("Index %d is out of bounds - nothing done",index));
}

//____________________________________________________________________________
void AliITSMFTCalibrationPix::GetBadPixelSingle(Int_t index, UInt_t &row, UInt_t &col) const 
{
  // i: is the i-th bad pixel in single bad pixel list
  // row: is the corresponding row (-1 if i is out of range)
  // col: is the corresponding column (-1 if i is out of range)
  if(index<0 && index>=GetNrBadSingle())  AliFatal(Form("Index %d >= NrBadSingle=%d",index,fNrBadSingle));
  col = fBadChannels.At(index*2);
  row = fBadChannels.At(index*2+1);
}

//___________________________________________________________________________
Int_t  AliITSMFTCalibrationPix::GetNrBad() const 
{
  // Total number of bad pixels (including bad chips) in a given chip
  Int_t bad=0;
  // single pixels:
  bad += fNrBadSingle;
  // whole chips:
  for (int chip=fNChips; chip--;) if (IsChipMarkedBad(chip)) bad += fNColPerChip*fNRow;
  return bad;
}

//___________________________________________________________________________
Int_t  AliITSMFTCalibrationPix::GetNrBadInChip(Int_t chip) const 
{
  // Total number of bad pixels (including bad chips) in a given chip
  if(chip<0 || chip>=fNChips) {AliError("Wrong chip number"); return -1;}
  if (IsChipMarkedBad(chip)) return fNColPerChip*fNRow;
  else {
    Int_t bad=0;
    for (int i=fNrBadSingle; i--;) {
      Int_t col = GetBadColAt(i);
      if (col!=-1) if (GetChipIndexFromCol(col)==chip) bad++;
    }
    return bad;
  }
}

//___________________________________________________________________________
Int_t  AliITSMFTCalibrationPix::GetNrBadInColumn(Int_t col) const 
{
  // Total number of bad pixels (including bad chips) in a given column: col. range
  if(col<0 || col>=fNCol) {AliError("Wrong column number"); return -1;}
  if (IsChipMarkedBad(GetChipIndexFromCol(col))) return fNRow;
  else {
    Int_t bad=0;
    for (int i=fNrBadSingle; i--;) if (GetBadColAt(i)==col) bad++;
    return bad;
  }
}

//______________________________________________________________________
Bool_t AliITSMFTCalibrationPix::IsBad() const 
{
  // Are all chips of this chip bad?
  for (Int_t chip=fNChips; chip--;) if (!IsChipMarkedBad(chip)) return kFALSE;
  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSMFTCalibrationPix::IsChipBad(Int_t chip) const 
{
  // Is the full chip bad?
  return (GetNrBadInChip(chip)==fNColPerChip*fNRow);
}

//______________________________________________________________________
Bool_t AliITSMFTCalibrationPix::IsColumnBad(Int_t col) const 
{
  // Is the full column bad?
  return (GetNrBadInColumn(col)==fNRow);
}

//____________________________________________________________________________
Bool_t AliITSMFTCalibrationPix::IsPixelBad(Int_t col, Int_t row) const 
{
  // Is this pixel bad?
  if(col<0 || col>=fNCol) {AliError("Wrong column number"); return kFALSE;}
  Int_t chip = GetChipIndexFromCol(col);
  if (IsChipMarkedBad(chip)) return kTRUE;
  for (Int_t i=fNrBadSingle; i--;) if (GetBadColAt(i)==col && GetBadRowAt(i)==row) return kTRUE;
  return kFALSE;
}

//______________________________________________________________________
Int_t AliITSMFTCalibrationPix::GetChipIndexFromCol(Int_t col) const 
{
  // returns chip index for specific column
  if(col>=fNCol) {AliWarning("Wrong column number"); return -1;}
  return col/fNColPerChip;
}

//______________________________________________________________________
void AliITSMFTCalibrationPix::SetNrBad(Int_t /*nr*/) 
{
  // should not be used anymore !!!
  AliError("This method should not be used anymore. Use SetNrBadSingle instead!!!");
}

//____________________________________________________________________________
void  AliITSMFTCalibrationPix::SetColRowData(Short_t nchip, Short_t ncolperchip, Short_t nrow) 
{
  // set segmentation data
  fNChips = nchip; 
  fNCol   = nchip*ncolperchip; 
  fNRow   = nrow; 
  fNColPerChip = ncolperchip;
}
