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

AliTOFDigitMap::AliTOFDigitMap()
{
//
// Default ctor
//
  fDigitMap = 0;

  fTOFGeometry = new AliTOFGeometry();

  fNSector = AliTOFGeometry::NSectors();
  fNplate = AliTOFGeometry::NPlates();
  fNstrip = fTOFGeometry->NStripC();//fTOFGeometry->NMaxNstrip();
  fNpx  = AliTOFGeometry::NpadX();
  fNpz  = AliTOFGeometry::NpadZ();
  fMaxIndex=fNSector*fNplate*fNstrip*fNpx*fNpz;

  fDigitMap = new Int_t*[fMaxIndex];
  for (Int_t i=0; i<fMaxIndex; i++) fDigitMap[i] = new Int_t[kMaxDigitsPerPad];
  Clear();
}

////////////////////////////////////////////////////////////////////////
AliTOFDigitMap::AliTOFDigitMap(const AliTOFDigitMap & /*digitMap*/)
:TObject()
{
//
// Dummy copy constructor
//
  ;

}

 
////////////////////////////////////////////////////////////////////////
AliTOFDigitMap::~AliTOFDigitMap()
{
//
// Destructor
//
  if (fDigitMap) {
    for (Int_t i=0; i<fMaxIndex; i++)  delete[] fDigitMap[i];
  }

  fTOFGeometry = 0;

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
Int_t AliTOFDigitMap::CheckedIndex(Int_t *vol) const
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

  if (CheckedIndex(vol)==-1) return -1;
  
  Int_t dummy = fDigitMap[CheckedIndex(vol)][label];
  
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
AliTOFDigitMap & AliTOFDigitMap::operator = (const AliTOFDigitMap & /*rhs*/) 
{
// Dummy assignment operator
    return *this;
}
