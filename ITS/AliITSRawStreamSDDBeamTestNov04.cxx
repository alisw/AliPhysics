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
///
/// This class provides access to ITS SDD digits in test beam raw data.
//  for beam test of November 2004
///
///////////////////////////////////////////////////////////////////////////////

/*
	Error Flag words: (http://www.bo.infn.it/~falchier/alice.html)
        with multi-event buffer
	   
bits 31-14: all zeros
bit  13   : L0 ack
bit  12   : L1 reject ack
bit  11   : L2 reject ack
bit  10   : prepulse ack
bit   9   : testpulse ack
bit   8   : flush
bit   7   : busy
bit   6   : flag error ch 1
bit   5   : flag error ch 0
bit   4   : disable trigger mismatch ack
bit   3   : parity error right hybrid
bit   2   : parity error left hybrid
bit   1   : parity error CARLOS ch 1
bit   0   : parity error CARLOS ch 2
*/
#include "AliITSRawStreamSDDBeamTestNov04.h"
#include "AliRawReader.h"

ClassImp(AliITSRawStreamSDDBeamTestNov04)





AliITSRawStreamSDDBeamTestNov04::AliITSRawStreamSDDBeamTestNov04(AliRawReader* rawReader) :
  AliITSRawStreamSDDBeamTest(rawReader)
{
// create an object to read ITS SDD raw digits


  fRawReader->Reset();
  fRawReader->SelectEquipment(17, 204, 204);
}


Bool_t AliITSRawStreamSDDBeamTestNov04::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left
  // skip the first 8 words
  while (fSkip < 9) {
    if (!fRawReader->ReadNextInt(fData)) return kFALSE;
    if ((fData >> 30) == 0x01) continue;  // JTAG word
    fSkip++;
  }

  Int_t countFoot=0;	
  while (kTRUE) {
    if ((fChannel < 0) || (fLastBit[0][fChannel] < fReadBits[0][fChannel])) {
      if (!fRawReader->ReadNextInt(fData)) return kFALSE;  // read next word
      fChannel = -1;
      if ((fData >> 28) == 0x02) {           // header
	fEventId = (fData >> 3) & 0x07FF;
	fCarlosId = (fData >> 1) & 0x03;
      } else if ((fData >> 28) == 0x03) {    // footer
        countFoot++; // stop before the last word (last word=jitter)
        if(countFoot==3) return kFALSE;	 
      } else if ((fData >> 29) == 0x00) {    // error

	if ((fData & 0x00000163) != 0) {
	  Error("Next", "error codes = %8.8x",fData);
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


Int_t AliITSRawStreamSDDBeamTestNov04::ReadJitter() {

  // Reads the value of the jitter between L0 and pascal stop
  // written in the last word of the buffer

     if (!fRawReader->ReadNextInt(fData)){
       Error("ReadJitter","Jitter word not found!!");
       return -1;  // read last word
     }
     if ( (fData&0xff000000) != 0xff000000) {
       Error("ReadJitter","wrong mask on Jitter word (0xffxxxxxx): %8.8x",fData);
       return -1;  // read last word
     }
     fJitter = fData&0x000000ff;
     if (fJitter<0x7 || fJitter>0xe) {
       Warning("ReadJitter","Unexpected jitter value %2.2x (%8.8x)",fJitter,fData);
       return fJitter;  // read last word
     }

     if (fRawReader->ReadNextInt(fData)){
       Error("ReadJitter","The equipment payload contains something after jitter");
       return -1;  // read last word
     }
     return fJitter;
}


