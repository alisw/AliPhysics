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

///////////////////////////////////////////////////////////////////////////////
//
// This is a class for reading ITS SPD raw data files and providing
// information about digits
//
///////////////////////////////////////////////////////////////////////////////

#include "AliITSRawStreamSPD.h"
#include "AliRawReader.h"

ClassImp(AliITSRawStreamSPD)


AliITSRawStreamSPD::AliITSRawStreamSPD(AliRawReader* rawReader) :
  AliITSRawStream(rawReader)
{
// create an object to read ITS SPD raw digits

  fRawReader->Select(1);
}


Bool_t AliITSRawStreamSPD::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  fPrevModuleID = fModuleID;
  while (fRawReader->ReadNextShort(fData)) {
  
    if ((fData & 0xE000) == 0x6000) {           // header
      fHitCount = 0;
      UShort_t halfStave = (fData >> 4) & 0x007F;
      UShort_t chipAddr = fData & 0x000F;
      fModuleID = 2 * halfStave;
      if (chipAddr >= 5) fModuleID++;
      fOffset = 32 * (chipAddr % 5);
    } else if ((fData & 0xE000) == 0x0000) {    // trailer
      UShort_t hitCount = fData & 0x1FFF;
      if (hitCount != fHitCount) Error("Next", "wrong number of hits!");
    } else if ((fData & 0xC000) == 0x8000) {    // pixel hit
      fHitCount++;
      fCoord1 = (fData & 0x001F) + fOffset;
      fCoord2 = (fData >> 5) & 0x00FF;
      return kTRUE;
    } else {                                    // fill word
      if (fData != 0xFEDC) Error("Next", "wrong fill word!");
    }

  }

  return kFALSE;
}
