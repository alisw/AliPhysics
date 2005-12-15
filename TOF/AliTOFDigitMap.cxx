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
// digitmap enables fast check if the pad was already hit
// The index of a AliTOFdigit is saved in the each hitmap "cell"
// (there is an offset +1, because the index can be zero and 
// zero means empty cell. 
// In TOF, number of strips varies according plate type, the highest
// number is in plate C. For all plates is used this number, so
// the size of the digitmap is a little bit greater than necessary, but
// it simplifies the access algorithm. 
// 
//
// Author: F. Pierella based on AliTOFHitMap
//
////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>

#include "AliTOFDigitMap.h"
#include "AliTOFdigit.h"
#include "AliTOFGeometry.h"


#include <TClonesArray.h>

ClassImp(AliTOFDigitMap)

AliTOFDigitMap::AliTOFDigitMap()
{
//
// Default ctor
//
  fDigitMap = 0;
  fDigits = 0;

  fTOFGeometry = new AliTOFGeometry();

}

////////////////////////////////////////////////////////////////////////
AliTOFDigitMap::AliTOFDigitMap(TClonesArray *dig, AliTOFGeometry *tofGeom)
{
  //
  // ctor
  //  
  // of course, these constants must not be hardwired
  // change later
  
  fTOFGeometry = tofGeom;

  fNSector = AliTOFGeometry::NSectors();
  fNplate = AliTOFGeometry::NPlates();
  fNstrip = fTOFGeometry->NMaxNstrip();
  fNpx  = AliTOFGeometry::NpadX();
  fNpz  = AliTOFGeometry::NpadZ();
  fMaxIndex=fNSector*fNplate*fNstrip*fNpx*fNpz;
  fDigitMap = new Int_t[fMaxIndex];
  fDigits =  dig;
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
  if (fDigitMap) delete[] fDigitMap;

  delete fTOFGeometry;

}

////////////////////////////////////////////////////////////////////////
void AliTOFDigitMap::Clear(const char *)
{
//
// Clear hitmap
//
    memset(fDigitMap,0,sizeof(int)*fMaxIndex);
}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFDigitMap::CheckedIndex(Int_t *vol) const
{
//
// Return checked indices for vol
//
  Int_t index=
    vol[0]*fNplate*fNstrip*fNpx*fNpz+             // sector
    vol[1]*fNstrip*fNpx*fNpz+                     // plate
    vol[2]*fNpx*fNpz+                             // strip
    vol[3]*fNpz+                                  // padx
    vol[4];                                       // padz

    if (index >= fMaxIndex) {
      Error("AliTOFDigitMap","CheckedIndex - input outside bounds");
	return -1;
    } else {
	return index;
    }
}

////////////////////////////////////////////////////////////////////////
void  AliTOFDigitMap::SetHit(Int_t *vol, Int_t idigit)
{
//
// Assign digit to pad vol
//

// 0 means empty pad, we need to shift indeces by 1
    fDigitMap[CheckedIndex(vol)]=idigit+1;
}

////////////////////////////////////////////////////////////////////////
void  AliTOFDigitMap::SetHit(Int_t *vol)
{
//
// Assign last digit to pad vol 
//

// 0 means empty pad, we need to shift indeces by 1
    fDigitMap[CheckedIndex(vol)]=fDigits->GetLast()+1;
}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFDigitMap::GetHitIndex(Int_t *vol) const
{
//
// Get contents of pad vol
//

// 0 means empty pad, we need to shift indeces by 1
    return fDigitMap[CheckedIndex(vol)]-1;
}

////////////////////////////////////////////////////////////////////////
TObject* AliTOFDigitMap::GetHit(Int_t *vol) const
{
//
// Get pointer to object at vol
// return 0 if vol out of bounds
    Int_t index=GetHitIndex(vol);
    return (index <0) ? 0 : fDigits->UncheckedAt(index);
}

////////////////////////////////////////////////////////////////////////
FlagType AliTOFDigitMap::TestHit(Int_t *vol) const
{
//
// Check if hit cell is empty, used or unused
//
    Int_t inf=fDigitMap[CheckedIndex(vol)];
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
