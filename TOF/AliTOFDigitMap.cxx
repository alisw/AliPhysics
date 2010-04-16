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

/*
$Log$
Revision 1.12  2007/02/20 15:57:00  decaro
Raw data update: to read the TOF raw data defined in UNPACKED mode

*/

////////////////////////////////////////////////////////////////////////
//
// AliTOFDigitMap class
//
// digitmap enables fast check if the pad was already digit.

// The index of a AliTOFdigit is saved in the each digitmap "cell"
// (there is an offset +1, because the index can be zero and zero
// means empty cell).
// In TOF, number of strips varies according plate type, the highest
// number is in plate C. For all plates is used this number, so the
// size of the digitmap is a little bit greater than necessary, but it
// simplifies the access algorithm.
// 
//
// Author: F. Pierella based on AliTOFHitMap
//
// Modified by A. De Caro
//
///////////////////////////////////////////////////////////////////////

#include "AliLog.h"

#include "AliTOFDigitMap.h"
#include "AliTOFGeometry.h"

ClassImp(AliTOFDigitMap)

AliTOFDigitMap::AliTOFDigitMap():
  fNSector(AliTOFGeometry::NSectors()),
  fNplate(AliTOFGeometry::NPlates()),
  fNstrip(AliTOFGeometry::NStripC()),
  fNpx(AliTOFGeometry::NpadX()),
  fNpz(AliTOFGeometry::NpadZ()),
  fMaxIndex(-1),
  fDigitMap(0x0)
{
//
// Default ctor
//

  fMaxIndex=fNSector*fNplate*fNstrip*fNpx*fNpz;
  fDigitMap = new Int_t*[fMaxIndex];

  for (Int_t i=0; i<fMaxIndex; i++) fDigitMap[i] = new Int_t[kMaxDigitsPerPad];
  Clear();
}

////////////////////////////////////////////////////////////////////////
AliTOFDigitMap::AliTOFDigitMap(const AliTOFDigitMap & digitMap):
  TObject(digitMap),
  fNSector(digitMap.fNSector),
  fNplate(digitMap.fNplate),
  fNstrip(digitMap.fNstrip),
  fNpx(digitMap.fNpx),
  fNpz(digitMap.fNpz),
  fMaxIndex(digitMap.fMaxIndex),
  fDigitMap(0x0)
{
//
// dummy copy constructor
//
  
  fMaxIndex=fNSector*fNplate*fNstrip*fNpx*fNpz;
  fDigitMap = new Int_t*[fMaxIndex];

  for (Int_t i=0; i<fMaxIndex; i++) fDigitMap[i] = new Int_t[kMaxDigitsPerPad];
}

////////////////////////////////////////////////////////////////////////
AliTOFDigitMap &
AliTOFDigitMap::operator=(const AliTOFDigitMap & /*digitMap*/)
{
//
// dummy copy const
//
    return *this;
}

 
////////////////////////////////////////////////////////////////////////
AliTOFDigitMap::~AliTOFDigitMap()
{
//
// Destructor
//
  if (fDigitMap) {
    for (Int_t i=0; i<fMaxIndex; i++)  delete[] fDigitMap[i];
    delete [] fDigitMap;
  }


}

////////////////////////////////////////////////////////////////////////
void AliTOFDigitMap::Clear(const Option_t*)
{
  //
  // Clear digitmap
  //

  for(Int_t ii=0; ii<fMaxIndex; ii++) {
    for (Int_t jj=0; jj<kMaxDigitsPerPad; jj++) {
      fDigitMap[ii][jj] = 0;
    }
  }
 
}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFDigitMap::CheckedIndex(Int_t * const vol) const
{
  //
  // Return checked index for vol
  //

  Int_t index =
    vol[0]*fNplate*fNstrip*fNpx*fNpz+             // sector
    vol[1]*fNstrip*fNpx*fNpz+                     // plate
    vol[2]*fNpx*fNpz+                             // strip
    vol[3]*fNpz+                                  // padx
    vol[4];                                       // padz

    if (index >= fMaxIndex || index < 0) {
      AliError("CheckedIndex - input outside bounds");
      return -1;
    } else {
      return index;
    }
}

////////////////////////////////////////////////////////////////////////
void AliTOFDigitMap::AddDigit(Int_t *vol, Int_t idigit)
{
  //
  // Assign digit to pad vol
  //
  // 0 means empty pad, we need to shift indeces by 1

  if (fDigitMap[CheckedIndex(vol)][kMaxDigitsPerPad-1]!=0) {
    AliDebug(1,Form("In the volume (Se%i, Pl%i, St%i, PadR%i, Pad%i) there is not more possibility to add other digits.", vol[0], vol[1], vol[2], vol[4], vol[3]));
    AliDebug(1,Form("Then, the digit number %i will be not inserted in the digit map, i.e. it will be lost.", idigit));
    AliDebug(1,Form("Please, check the possibility to increase the digit map size (succently set to %i)", kMaxDigitsPerPad));
    return;
  }

  for (Int_t slot=0; slot<kMaxDigitsPerPad; slot++) {

    if (fDigitMap[CheckedIndex(vol)][slot]==0) {
      fDigitMap[CheckedIndex(vol)][slot]=idigit+1;
      break;
    }
    //else continue;

  }

}

////////////////////////////////////////////////////////////////////////
void AliTOFDigitMap::GetDigitIndex(Int_t *vol, Int_t *digitLabels) const
{
  //
  // Get all contents (digitLabels) of pad volume (vol)
  //

  // 0 means empty pad, we need to shift indeces by 1

  Int_t dummy;
    for (Int_t j=0; j<kMaxDigitsPerPad; j++) {
      dummy = GetDigitIndex(vol,j);
      if (dummy>=0) digitLabels[j] = dummy;
      else break;
    }
}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFDigitMap::GetDigitIndex(Int_t *vol, Int_t label) const
{
  //
  // Get one of the contents (label) of pad volume (vol)
  //

  // 0 means empty pad, we need to shift indeces by 1

  if (!(label<kMaxDigitsPerPad)) {
    AliWarning(Form("label (=%i) >= kMaxDigitsPerPad (=%i)", label, kMaxDigitsPerPad));
    return -1;
  }

  Int_t ci = CheckedIndex(vol);
  if (ci==-1) return -1;
  
  Int_t dummy = fDigitMap[ci][label];
  
  if (dummy>0) return dummy-1;
  else return -1;

}

////////////////////////////////////////////////////////////////////////
FlagType AliTOFDigitMap::TestDigit(Int_t *vol) const
{
//
// Check if hit cell is empty, used or unused
//
  Int_t inf=fDigitMap[CheckedIndex(vol)][0]; // to be modified
    if (inf > 0) {
	return kUsed;
    } else if (inf == 0) {
	return kEmpty;
    } else {
	return kUnused;
    }
}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFDigitMap::GetFilledCellNumber() const
{
  //
  // Returns the number of filled cells of the TOF digit map
  //

  Int_t counter = 0;

  for (Int_t index = 0; index < fMaxIndex; ++index)
  {
    for (Int_t label = 0; label < kMaxDigitsPerPad; ++label)
    {
      if (fDigitMap[index][label] > 0)
      {
	++counter;
	break;
      }
    }
  }

  return counter;
}

////////////////////////////////////////////////////////////////////////
Bool_t AliTOFDigitMap::StripDigitCheck(Int_t iSector, Int_t iPlate, Int_t iStrip) const
{
  //
  // Returns:
  //           kFALSE if the strip doesn't contain digits
  //           kTRUE  if the strip contains at least one digit
  //

  Int_t volume[5] = {iSector, iPlate, iStrip, -1, -1};
  Bool_t counter = kFALSE;

  for (Int_t iPadX=0; iPadX<fNpx; iPadX++)
    for (Int_t iPadZ=0; iPadZ<fNpz; iPadZ++)
      {
	volume[3] = iPadX;
	volume[4] = iPadZ;
	for (Int_t label=0; label<kMaxDigitsPerPad; label++) {
	  if (GetDigitIndex(volume, label)>=0) {
	    counter = kTRUE;
	    break;
	  }
	}
      }

  return counter;

}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFDigitMap::DigitInStrip(Int_t iSector, Int_t iPlate, Int_t iStrip) const
{
  //
  // Returns number of digits in the strip iStrip,
  //         in the plate iPlate of the sector iSector
  //

  Int_t volume[5] = {iSector, iPlate, iStrip, -1, -1};
  Int_t counter = 0;

  for (Int_t iPadX=0; iPadX<fNpx; iPadX++)
    for (Int_t iPadZ=0; iPadZ<fNpz; iPadZ++)
      for (Int_t label=0; label<kMaxDigitsPerPad; label++) {
	volume[3] = iPadX;
	volume[4] = iPadZ;
	if (GetDigitIndex(volume, label)>=0)
	  counter++;
      }

  return counter;

}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFDigitMap::FilledCellsInStrip(Int_t iSector, Int_t iPlate, Int_t iStrip) const
{
  //
  // Returns number of digits in the strip iStrip,
  //         in the plate iPlate of the sector iSector
  //

  Int_t volume[5] = {iSector, iPlate, iStrip, -1, -1};
  Int_t counter = 0;

  for (Int_t iPadX=0; iPadX<fNpx; iPadX++)
    for (Int_t iPadZ=0; iPadZ<fNpz; iPadZ++) {
      volume[3] = iPadX;
      volume[4] = iPadZ;
      if (GetDigitIndex(volume, 0)>=0)
	counter++;
    }

  return counter;

}

////////////////////////////////////////////////////////////////////////
void AliTOFDigitMap::ResetDigitNumber(Int_t *vol, Int_t dig)
{
  //
  // Reset digit into pad vol
  //

  for (Int_t slot=0; slot<kMaxDigitsPerPad; slot++) {
    if (fDigitMap[CheckedIndex(vol)][slot]-1==dig) {
      fDigitMap[CheckedIndex(vol)][slot] = 0;
    }
  }

}

////////////////////////////////////////////////////////////////////////
void AliTOFDigitMap::ResetDigit(Int_t *vol, Int_t dig)
{
  //
  // Reset digit into pad vol
  //
  // 0 means empty pad, we need to shift indeces by 1

  fDigitMap[CheckedIndex(vol)][dig] = 0;

}

////////////////////////////////////////////////////////////////////////
void AliTOFDigitMap::ResetDigit(Int_t *vol)
{
  //
  // Reset digit into pad vol
  //
  // 0 means empty pad, we need to shift indices by 1

  for (Int_t slot=0; slot<kMaxDigitsPerPad; slot++)
    fDigitMap[CheckedIndex(vol)][slot] = 0;

}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFDigitMap::GetNumberOfDigits(Int_t *vol)
{
  //
  // Returns the number of digit
  //   into pad volume vol
  //
  // 0 means empty pad
  //

  Int_t counter = 0;

  for (Int_t slot=0; slot<kMaxDigitsPerPad; slot++)
    if (GetDigitIndex(vol, slot)>=0) counter++;

  return counter;

}
