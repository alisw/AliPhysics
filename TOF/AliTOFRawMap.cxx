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

/////////////////////////////////////////////////////////////////////////
//                                                                     //
// AliTOFRawMap class                                                  //
//                                                                     //
// It enables fast check if the TDC channel was already engaged        //
// for a measurement.                                                  //
// The index of a AliTOFrawData is saved in the each rawdatamap "cell" //
// (there is an offset +1, because the index can be zero and           //
// zero means empty cell.                                              //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include "TClonesArray.h"

#include "AliLog.h"

#include "AliTOFGeometry.h"
#include "AliTOFRawMap.h"

ClassImp(AliTOFRawMap)

AliTOFRawMap::AliTOFRawMap():
  TObject(),
  fNtrm(-1),
  fNtrmChain(-1),
  fNtdc(-1),
  fNtdcChannel(-1),
  fRawData(0x0),
  fMaxIndex(-1),
  fRawMap(0x0)
{
//
// Default ctor
//
}

////////////////////////////////////////////////////////////////////////
AliTOFRawMap::AliTOFRawMap(TClonesArray *dig)://, AliTOFGeometry *tofGeom:
  TObject(),
  fNtrm(AliTOFGeometry::NTRM()+2),
  fNtrmChain(AliTOFGeometry::NChain()),
  fNtdc(AliTOFGeometry::NTdc()),
  fNtdcChannel(AliTOFGeometry::NCh()),
  fRawData(dig),
  fMaxIndex(-1),
  fRawMap(0x0)
{
//
// ctor
//

// of course, these constants must not be hardwired
// change later

  fMaxIndex = fNtrm*fNtrmChain*fNtdc*fNtdcChannel;
  fRawMap = new Int_t[fMaxIndex];
  Clear();
}

////////////////////////////////////////////////////////////////////////
AliTOFRawMap::AliTOFRawMap(const AliTOFRawMap & rawMap)
  :TObject(rawMap),
  fNtrm(rawMap.fNtrm),
  fNtrmChain(rawMap.fNtrmChain),
  fNtdc(rawMap.fNtdc),
  fNtdcChannel(rawMap.fNtdcChannel),
  fRawData(rawMap.fRawData),
  fMaxIndex(-1),
  fRawMap(0x0)
{
//
// Dummy copy constructor
//

  fMaxIndex = fNtrm*fNtrmChain*fNtdc*fNtdcChannel;
  fRawMap = new Int_t[fMaxIndex];
}

////////////////////////////////////////////////////////////////////////
AliTOFRawMap &
AliTOFRawMap::operator=(const AliTOFRawMap & /*rawMap*/)
{
//
// Dummy copy constructor
//
  return *this;
}

 
////////////////////////////////////////////////////////////////////////
AliTOFRawMap::~AliTOFRawMap()
{
//
// Destructor
//
  if (fRawMap)
    delete[] fRawMap;

}

////////////////////////////////////////////////////////////////////////
void AliTOFRawMap::Clear(const char *)
{
//
// Clear hitmap
//
    memset(fRawMap,0,sizeof(int)*fMaxIndex);
}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFRawMap::CheckedIndex(Int_t *slot) const
{
//
// Return checked indices for vol
//
  Int_t index =
    slot[0]*fNtrmChain*fNtdc*fNtdcChannel + // TRM
    slot[1]*fNtdc*fNtdcChannel +            // TRM chain
    slot[2]*fNtdcChannel +                  // TDC
    slot[3];                                // TDC channel

    if (index >= fMaxIndex) {
      AliError("CheckedIndex - input outside bounds");
	return -1;
    } else {
	return index;
    }
}

////////////////////////////////////////////////////////////////////////
void  AliTOFRawMap::SetHit(Int_t *slot, Int_t idigit)
{
//
// Assign digit to pad vol
//

// 0 means empty pad, we need to shift indeces by 1
    fRawMap[CheckedIndex(slot)]=idigit+1;
}

////////////////////////////////////////////////////////////////////////
void  AliTOFRawMap::SetHit(Int_t *slot)
{
//
// Assign last digit to channel slot
//

// 0 means empty pad, we need to shift indeces by 1
    fRawMap[CheckedIndex(slot)]=fRawData->GetLast()+1;
}

////////////////////////////////////////////////////////////////////////
Int_t AliTOFRawMap::GetHitIndex(Int_t *slot) const
{
//
// Get contents of channel slot
//

// 0 means empty pad, we need to shift indeces by 1
    return fRawMap[CheckedIndex(slot)]-1;
}

////////////////////////////////////////////////////////////////////////
TObject* AliTOFRawMap::GetHit(Int_t *slot) const
{
//
// Get pointer to object at alot
// return 0 if vol out of bounds
    Int_t index = GetHitIndex(slot);
    return (index <0) ? 0 : fRawData->UncheckedAt(index);
}

////////////////////////////////////////////////////////////////////////
FlagType AliTOFRawMap::TestHit(Int_t *slot) const
{
//
// Check if hit cell is empty, used or unused
//
    Int_t inf = fRawMap[CheckedIndex(slot)];
    if (inf > 0) {
	return kUsed;
    } else if (inf == 0) {
	return kEmpty;
    } else {
	return kUnused;
    }
}

