/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: $*/

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SDD digits in beam test raw data
///
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliITSRawStreamSDDBeamTest.h"
#include "AliRawReader.h"


ClassImp(AliITSRawStreamSDDBeamTest)
  
const UInt_t AliITSRawStreamSDDBeamTest::fgkCodeLength[8] =  {8, 18, 2, 3, 4, 5, 6, 7};

//______________________________________________________________________
AliITSRawStreamSDDBeamTest::AliITSRawStreamSDDBeamTest(AliRawReader* rawReader) :
  AliITSRawStream(rawReader),
fData(0),
fSkip(0),
fEventId(0),
fCarlosId(-1),
fChannel(0),
fJitter(0)
{
// create an object to read ITS SDD raw digits
  fRawReader->Reset();
  fRawReader->SelectEquipment(17, 204, 204);


}

//______________________________________________________________________
UInt_t AliITSRawStreamSDDBeamTest::ReadBits()
{
// read bits from the given channel
  UInt_t result = (fChannelData[fCarlosId][fChannel] & ((1<<fReadBits[fCarlosId][fChannel]) - 1));
  fChannelData[fCarlosId][fChannel] >>= fReadBits[fCarlosId][fChannel]; 
  fLastBit[fCarlosId][fChannel] -= fReadBits[fCarlosId][fChannel];
  return result;
}

//______________________________________________________________________
Int_t AliITSRawStreamSDDBeamTest::DecompAmbra(Int_t value) const
{
  // AMBRA decompression (from 8 to 10 bit)
  
  if ((value & 0x80) == 0) {
    return value & 0x7f;
  } else if ((value & 0x40) == 0) {
    return 0x081 + ((value & 0x3f) << 1);
  } else if ((value & 0x20) == 0) {
    return 0x104 + ((value & 0x1f) << 3);
  } else {
    return 0x208 + ((value & 0x1f) << 4);
  }
  
}

//______________________________________________________________________
Bool_t AliITSRawStreamSDDBeamTest::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  // skip the first 8 words
  while (fSkip < 8) {
    if (!fRawReader->ReadNextInt(fData)) return kFALSE;
    if ((fData >> 30) == 0x01) continue;  // JTAG word
    if (fSkip == 4) {
      if (fData != 0) {
	Error("Next", "data not valid: %8.8d", fData);
	return kFALSE;
      }
    }
    fSkip++;
  }

  while (kTRUE) {
    if ((fChannel < 0) || (fLastBit[0][fChannel] < fReadBits[0][fChannel])) {
      if (!fRawReader->ReadNextInt(fData)) return kFALSE;  // read next word

      fChannel = -1;
      if ((fData >> 28) == 0x02) {           // header
	fEventId = (fData >> 3) & 0x07FF;
	fCarlosId = (fData >> 1) & 0x03;
      } else if ((fData >> 28) == 0x03) {    // footer
	// ignored
      } else if ((fData >> 29) == 0x00) {    // error
	if ((fData & 0x1FFFFFFF) != 0) {
	  Error("Next", "error codes = %x, %x\n", 
		(fData >> 0) & 0x3FFF, (fData >> 14) & 0x3FFF);
	  return kFALSE;
	}
      } else if ((fData >> 30) == 0x01) {    // JTAG word
	// ignored
      } else if ((fData >> 30) == 0x02) {    // channel 0 data
	fChannel = 0;
      } else if ((fData >> 30) == 0x03) {    // channel 1 data
	fChannel = 1;
      } else {                               // unknown data format
	Error("Next", "invalid data: %8.8x\n", fData);
	return kFALSE;
      }

      if (fChannel >= 0) {          // add read word to the data
	fChannelData[0][fChannel] += 
	  (ULong64_t(fData & 0x3FFFFFFF) << fLastBit[0][fChannel]);
	fLastBit[0][fChannel] += 30;
      }

    } else {  // decode data
      if (fReadCode[0][fChannel]) {    // read the next code word
	fChannelCode[0][fChannel] = ReadBits();
	fReadCode[0][fChannel] = kFALSE;
	fReadBits[0][fChannel] = fgkCodeLength[fChannelCode[0][fChannel]];

      } else {                      // read the next data word
	UInt_t data = ReadBits();
	fReadCode[0][fChannel] = kTRUE;
	fReadBits[0][fChannel] = 3;
	if (fChannelCode[0][fChannel] == 0) {         // set the time bin
	  fTimeBin[0][fChannel] = data;
	} else if (fChannelCode[0][fChannel] == 1) {  // next anode
	  fTimeBin[0][fChannel] = 0;
	  fAnode[0][fChannel]++;
	} else {                                   // ADC signal data
	  fSignal = DecompAmbra(data + (1 << fChannelCode[0][fChannel]) + 
	    fLowThreshold[fChannel]);
	  fCoord1 = fAnode[0][fChannel];
	  fCoord2 = fTimeBin[0][fChannel];
	  fTimeBin[0][fChannel]++;
	  return kTRUE;
	}
      }
    }
  }

  return kFALSE;
}




