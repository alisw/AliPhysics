/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTLogging.h"
#include "AliHLTPHOSDigitReader.h"
#include "AliHLTPHOSDigitDataStruct.h"


AliHLTPHOSDigitReader::AliHLTPHOSDigitReader() :
  fDigitHeader(0),
  fCurrentDigit(0),
  fNextDigit(0),
  fPrevDigit(0),
  fFirstDigit(0)
{
  // See header file for documentation
}

AliHLTPHOSDigitReader::~AliHLTPHOSDigitReader()
{
  // See header file for documentation
}


AliHLTPHOSDigitDataStruct* AliHLTPHOSDigitReader::NextDigit()
{

  fPrevDigit = fCurrentDigit;
  fCurrentDigit = fNextDigit;

  if(fCurrentDigit == 0) return 0;

  if(fCurrentDigit->fMemOffsetNext != 0)
    {
      fNextDigit = reinterpret_cast<AliHLTPHOSDigitDataStruct*>(reinterpret_cast<UChar_t*>(fCurrentDigit) + fCurrentDigit->fMemOffsetNext);
    }
  else
    {
      fNextDigit = 0;
    }
  //  cout << "Digit is (x, z): " << fCurrentDigit->fX << ", " << fCurrentDigit->fZ << endl;
  return fCurrentDigit;
}

void AliHLTPHOSDigitReader::DropDigit()
{
  if(fCurrentDigit == fFirstDigit)
    {
      fFirstDigit = reinterpret_cast<AliHLTPHOSDigitDataStruct*>(reinterpret_cast<UChar_t*>(fFirstDigit) + fFirstDigit->fMemOffsetNext);
      fDigitHeader->fFirstDigitOffset += fCurrentDigit->fMemOffsetNext;
      //      HLTError("Dropping digit (x,z): %d, %d was first in list", fCurrentDigit->fX, fCurrentDigit->fZ);
    }
  else if(fCurrentDigit != 0)
    {
      fPrevDigit->fMemOffsetNext = fPrevDigit->fMemOffsetNext + fCurrentDigit->fMemOffsetNext;
      //      HLTError("Dropping digit (x,z): %d, %d, first digit is (x,z): %d, %d", fCurrentDigit->fX, fCurrentDigit->fZ, fFirstDigit->fX, fFirstDigit->fZ);
    }
  fCurrentDigit = fPrevDigit;
}
