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
#include "TBuffer.h"
#include "AliITSCalibrationSPD.h"

///////////////////////////////////////////////////////////////////////////
//  Calibration class for set:ITS                   
//  Specific subdetector implementation for         
//  Silicon pixels                                  
//
//  Modified by D. Elia, G.E. Bruno, H. Tydesjo
///////////////////////////////////////////////////////////////////////////

ClassImp(AliITSCalibrationSPD)

//______________________________________________________________________
AliITSCalibrationSPD::AliITSCalibrationSPD():
AliITSCalibration(),
fNrBad(0),
fBadChannels(0){
  // constructor

   SetDataType("simulated");
   ClearBad();
}
//____________________________________________________________________________
void AliITSCalibrationSPD::ClearBad() {
  // clear all bad pixels (single+chips)
  fBadChannels.Reset();
  fNrBad=0;
  for (UInt_t chip=0; chip<5; chip++) {
    fBadChip[chip]=kFALSE;
  }
}
//____________________________________________________________________________
void AliITSCalibrationSPD::AddBad(UInt_t col, UInt_t row) {
  // add single bad pixel 
  fBadChannels.Set(fNrBad*2+2);
  fBadChannels.AddAt(col,fNrBad*2);
  fBadChannels.AddAt(row,fNrBad*2+1);
  fNrBad++;
}
//____________________________________________________________________________
void AliITSCalibrationSPD::SetChipBad(UInt_t chip) {
  // set full chip bad
  if (chip>=5) {AliError("Wrong chip number");
  }
  else {
    fBadChip[chip]=kTRUE;
  }
}
//____________________________________________________________________________
void AliITSCalibrationSPD::UnSetChipBad(UInt_t chip) {
  // unset full chip bad
  if (chip>=5 ) {AliError("Wrong chip number");
  }
  else {
    fBadChip[chip]=kFALSE;
  }
}
//____________________________________________________________________________
Int_t AliITSCalibrationSPD::GetBadColAt(UInt_t index) const {
  // Get column of index-th bad pixel
  if ((Int_t)index<GetNrBadSingle()) {
    return fBadChannels.At(index*2);
  }
  else {
    Int_t badChipIndex=(index-GetNrBadSingle())/(32*256);
    Int_t badChipsFound =0;
    for (UInt_t chip=0; chip<5; chip++) {
      if (fBadChip[chip]) badChipsFound++;
      if (badChipIndex==badChipsFound-1) {
	Int_t badPixelIndex=(index-GetNrBadSingle())%(32*256);
	return chip*32 + badPixelIndex/256;
      }
    }
  }
  AliError(Form("Index %d is out of bounds - returning -1",index));
  return -1;
}
//____________________________________________________________________________
Int_t AliITSCalibrationSPD::GetBadRowAt(UInt_t index) const {
  // Get row of index-th bad pixel
  if ((Int_t)index<GetNrBadSingle()) {
    return fBadChannels.At(index*2+1);
  }
  else {
    Int_t badChipIndex=(index-GetNrBadSingle())/(32*256);
    Int_t badChipsFound =0;
    for (UInt_t chip=0; chip<5; chip++) {
      if (fBadChip[chip]) badChipsFound++;
      if (badChipIndex==badChipsFound-1) {
	Int_t badPixelIndex=(index-GetNrBadSingle())%(32*256);
	return badPixelIndex%256;
      }
    }
  }
  AliError(Form("Index %d is out of bounds - returning -1",index));
  return -1;
}
//____________________________________________________________________________
void AliITSCalibrationSPD::GetBadPixel(Int_t index, Int_t &row, Int_t &col) const {
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
      Int_t badChipIndex=(index-GetNrBadSingle())/(32*256);
      Int_t badChipsFound =0;
      for (UInt_t chip=0; chip<5; chip++) {
	if (fBadChip[chip]) badChipsFound++;
	if (badChipIndex==badChipsFound-1) {
	  Int_t badPixelIndex=(index-GetNrBadSingle())%(32*256);
	  col = chip*32 + badPixelIndex/256;
	  row = badPixelIndex%256;
	  return;
	}
      }
    }
  }
  AliError(Form("Index %d is out of bounds - nothing done",index));
}
//___________________________________________________________________________
Int_t  AliITSCalibrationSPD::GetNrBad() const {
  // Total number of bad pixels (including bad chips) in a given module
  Int_t bad=0;
  // single pixels:
  bad+=fNrBad;
  // whole chips:
  for (UInt_t chip=0; chip<5; chip++) {
    bad+=fBadChip[chip]*32*256;
  }
  return bad;
}
//___________________________________________________________________________
Int_t  AliITSCalibrationSPD::GetNrBadInChip(Int_t chip) const {
  // Total number of bad pixels (including bad chips) in a given chip: chip range [0,4]
  if(chip<0 || chip>4) {AliError("Wrong chip number"); return -1;}
  if (fBadChip[chip]) return 32*256;
  else {
    Int_t bad=0;
    for (UInt_t i=0; i<fNrBad; i++) {
      Int_t col = GetBadColAt(i);
      if (col!=-1) {
	if (GetChipIndexFromCol(col)==chip) bad++;
      }
    }
    return bad;
  }
}
//___________________________________________________________________________
Int_t  AliITSCalibrationSPD::GetNrBadInColumn(Int_t col) const {
  // Total number of bad pixels (including bad chips) in a given column: col. range [0,159]
  if(col<0 || col>159) {AliError("Wrong column number"); return -1;}
  if (fBadChip[GetChipIndexFromCol(col)]) return 256;
  else {
    Int_t bad=0;
    for (UInt_t i=0; i<fNrBad; i++) {
      if (GetBadColAt(i)==col) bad++;
    }
    return bad;
  }
}
//______________________________________________________________________
Bool_t AliITSCalibrationSPD::IsBad() const {
  // Are all chips of this module bad?
  for (UInt_t chip=0; chip<5; chip++) {
    if (!fBadChip[chip]) return kFALSE;
  }
  return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSCalibrationSPD::IsChipBad(Int_t chip) const {
  // Is the full chip bad?
  return (GetNrBadInChip(chip)==32*256);
}
//______________________________________________________________________
Bool_t AliITSCalibrationSPD::IsColumnBad(Int_t col) const {
  // Is the full column bad?
  return (GetNrBadInColumn(col)==256);
}
//____________________________________________________________________________
Bool_t AliITSCalibrationSPD::IsPixelBad(Int_t col, Int_t row) const {
  // Is this pixel bad?
  if(col<0 || col>159) {AliError("Wrong column number"); return kFALSE;}
  Int_t chip = GetChipIndexFromCol(col);
  if (fBadChip[chip]) return kTRUE;
  for (UInt_t i=0; i<fNrBad; i++) { 
    if (GetBadColAt(i)==col && GetBadRowAt(i)==row) {
      return kTRUE;
    }
  }
  return kFALSE;
}
//______________________________________________________________________
Int_t AliITSCalibrationSPD::GetChipIndexFromCol(UInt_t col) const {
  // returns chip index for specific column
  if(col>=160) {AliWarning("Wrong column number"); return -1;}
  return col/32;
}
//______________________________________________________________________
void AliITSCalibrationSPD::SetNrBad(UInt_t /*nr*/) {
  // should not be used anymore !!!
  AliError("This method should not be used anymore. Use SetNrBadSingle instead!!!");
}
//______________________________________________________________________
void AliITSCalibrationSPD::Streamer(TBuffer &ruub) {
  // Stream an object of class AliITSCalibrationSPD.
  UInt_t ruus, ruuc;
  if (ruub.IsReading()) {
    Version_t ruuv = ruub.ReadVersion(&ruus, &ruuc); if (ruuv) { }
    AliITSCalibration::Streamer(ruub);
    if (ruuv >= 8) {
      ruub >> fNrBad;
      fBadChannels.Streamer(ruub);
      ruub.ReadStaticArray((bool*)fBadChip);
    }
    else {
      Double_t dummy;
      ruub >> dummy;
      ruub >> dummy;
      ruub >> dummy;
      ruub >> dummy;
      ruub >> dummy;
      ruub >> dummy;
      ruub >> dummy;
      ruub >> fNrBad;
      if (ruuv == 7) {
	fBadChannels.Streamer(ruub);
	ruub.ReadStaticArray((bool*)fBadChip);
      }
      else {
	if (ruuv == 6) {
	  fBadChannels.Streamer(ruub);
	}
	else {
	  TArrayI fBadChannelsV1;
	  fBadChannelsV1.Streamer(ruub);
	  fBadChannels.Set(fNrBad*2);
	  for (UInt_t i=0; i<fNrBad*2; i++) {
	    fBadChannels[i] = fBadChannelsV1[i];
	  }
	}
	for (UInt_t i=0; i<5; i++) {
	  fBadChip[i]=kFALSE;
	}
      }
    }
    ruub.CheckByteCount(ruus, ruuc, AliITSCalibrationSPD::IsA());
  }
  else {
    ruuc = ruub.WriteVersion(AliITSCalibrationSPD::IsA(), kTRUE);
    AliITSCalibration::Streamer(ruub);
    ruub << fNrBad;
    fBadChannels.Streamer(ruub);
    ruub.WriteArray(fBadChip, 5);
    ruub.SetByteCount(ruuc, kTRUE);
  }
}
