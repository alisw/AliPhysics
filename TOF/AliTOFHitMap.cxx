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

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// AliTOFHitMap class
//
// hitmap enables fast check if the pad was already hit
// The index of a AliTOFSDigit is saved in the each hitmap "cell"
// (there is an offset +1, because the index can be zero and 
// zero means empty cell. 
// In TOF, number of strips varies according plate type, the highest
// number is in plate C. For all plates is used this number, so
// the size of the hitmap is a little bit greater than necessary, but
// it simplifies the access algorithm. 
// 
//
// Author: Jiri Chudoba (CERN), based on AliMUONHitMap
//
////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliTOFHitMap.h"
#include "AliTOFGeometry.h"


#include <TClonesArray.h>

ClassImp(AliTOFHitMap)

AliTOFHitMap::AliTOFHitMap():
  fNSector(-1),
  fNplate(-1),
  fNstrip(-1),
  fNpx(-1),
  fNpz(-1),
  fSDigits(0x0),
  fMaxIndex(-1),
  fHitMap(0x0)
{
//
// Default ctor
//
}

////////////////////////////////////////////////////////////////////////
AliTOFHitMap::AliTOFHitMap(TClonesArray *dig):
  fNSector(-1),
  fNplate(-1),
  fNstrip(-1),
  fNpx(-1),
  fNpz(-1),
  fSDigits(dig),
  fMaxIndex(-1),
  fHitMap(0x0)
{
//
// ctor
//

// of course, these constants must not be hardwired
// change later

  fNSector = AliTOFGeometry::NSectors();
  fNplate = AliTOFGeometry::NPlates();
  fNstrip = AliTOFGeometry::NStripC();//fTOFGeometry->NMaxNstrip();
  fNpx  = AliTOFGeometry::NpadX();
  fNpz  = AliTOFGeometry::NpadZ();
  fMaxIndex=fNSector*fNplate*fNstrip*fNpx*fNpz;
  fHitMap = new Int_t[fMaxIndex];
  Clear();
}

////////////////////////////////////////////////////////////////////////
AliTOFHitMap::AliTOFHitMap(const AliTOFHitMap & /*hitMap*/)
  :TObject(),
  fNSector(-1),
  fNplate(-1),
  fNstrip(-1),
  fNpx(-1),
  fNpz(-1),
  fSDigits(0x0),
  fMaxIndex(-1),
  fHitMap(0x0)
{
//
// Dummy copy constructor
//
  ;
}

 
////////////////////////////////////////////////////////////////////////
AliTOFHitMap::~AliTOFHitMap()
{
//
// Destructor
//
  delete[] fHitMap;

}

////////////////////////////////////////////////////////////////////////
void AliTOFHitMap::Clear(const char *)
{
//
// Clear hitmap
//
    memset(fHitMap,0,sizeof(int)*fMaxIndex);
}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFHitMap::CheckedIndex(Int_t *vol) const
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
      AliError("CheckedIndex - input outside bounds");
	return -1;
    } else {
	return index;
    }
}

////////////////////////////////////////////////////////////////////////
void  AliTOFHitMap::SetHit(Int_t *vol, Int_t idigit)
{
//
// Assign digit to pad vol
//

// 0 means empty pad, we need to shift indeces by 1
    fHitMap[CheckedIndex(vol)]=idigit+1;
}

////////////////////////////////////////////////////////////////////////
void  AliTOFHitMap::SetHit(Int_t *vol)
{
//
// Assign last digit to pad vol 
//

// 0 means empty pad, we need to shift indeces by 1
    fHitMap[CheckedIndex(vol)]=fSDigits->GetLast()+1;
}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFHitMap::GetHitIndex(Int_t *vol) const
{
//
// Get contents of pad vol
//

// 0 means empty pad, we need to shift indeces by 1
    return fHitMap[CheckedIndex(vol)]-1;
}

////////////////////////////////////////////////////////////////////////
TObject* AliTOFHitMap::GetHit(Int_t *vol) const
{
//
// Get pointer to object at vol
// return 0 if vol out of bounds
    Int_t index=GetHitIndex(vol);
    return (index <0) ? 0 : fSDigits->UncheckedAt(index);
}

////////////////////////////////////////////////////////////////////////
FlagType AliTOFHitMap::TestHit(Int_t *vol) const
{
//
// Check if hit cell is empty, used or unused
//
    Int_t inf=fHitMap[CheckedIndex(vol)];
    if (inf > 0) {
	return kUsed;
    } else if (inf == 0) {
	return kEmpty;
    } else {
	return kUnused;
    }
}

////////////////////////////////////////////////////////////////////////
AliTOFHitMap & AliTOFHitMap::operator = (const AliTOFHitMap & /*rhs*/) 
{
// Dummy assignment operator
    return *this;
}
