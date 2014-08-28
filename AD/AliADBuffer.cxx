/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

// Storing digits in a binary file
// according to the DDL mapping
// Author: B. Cheynis

#include <Riostream.h>
#include <TObjArray.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliRawDataHeaderSim.h"
#include "AliADBuffer.h"
#include "AliADdigit.h"
#include "AliADConst.h"

ClassImp(AliADBuffer)

//_____________________________________________________________________________
AliADBuffer::AliADBuffer():TObject(),
    fRemainingWord(0),
    f()
{
  //
  // default constructor
  //
}
//_____________________________________________________________________________
AliADBuffer::AliADBuffer(const char* fileName):TObject(),
    fRemainingWord(0),
    f()
{
  // Constructor
  f = new AliFstream(fileName);
  AliRawDataHeaderSim header;
  f->WriteBuffer((char*)(&header), sizeof(header));

}

//_____________________________________________________________________________
AliADBuffer::~AliADBuffer(){
  // Destructor, it closes the IO stream
  AliRawDataHeaderSim header;
  header.fSize = f->Tellp();
  header.SetAttribute(0);  // valid data
  f->Seekp(0);
  f->WriteBuffer((char*)(&header), sizeof(header));
  delete f;
}

//_____________________________________________________________________________
void AliADBuffer::WriteTriggerInfo(UInt_t trigger) {
  // The method writes AD trigger information
  // This info is contained in the first two
  // raw-data words following the raw-data header (CDH).

  f->WriteBuffer((char*)(&trigger),sizeof(trigger));

  // By default all the inputs are unmasked... Hopefully
  UInt_t triggerMask = 0xffff;
  f->WriteBuffer((char*)(&triggerMask),sizeof(triggerMask));
}

//_____________________________________________________________________________
void AliADBuffer::WriteTriggerScalers() {
  // The method writes the AD trigger scalers
  // For the moment there is no way to simulate
  // this, so we fill the necessary 16 words with 0

  // First the general trigger scalers (16 of them)
  for(Int_t i = 0; i < 16; i++) {
      UInt_t data = 0;
      f->WriteBuffer((char*)&data,sizeof(data));
  }
}

//_____________________________________________________________________________
void AliADBuffer::WriteBunchNumbers() {
  // The method writes the Bunch Numbers corresponding 
  // to the 10 Minimum Bias events
  // For the moment there is no way to simulate
  // this, so we fill the necessary 10 words with 0

  // First the bunch crossing numbers
  // for these 10 events
  
  for(Int_t i = 0; i < 10; i++) {
      UInt_t data = 0;
      f->WriteBuffer((char*)&data,sizeof(data));
  }

}

//_____________________________________________________________________________
void AliADBuffer::WriteChannel(Int_t channel, Short_t *adc, Bool_t integrator){
  // It writes AD charge information into a raw data file. 
  // Being called by Digits2Raw
  
  UShort_t data = 0;
  for(Int_t i = 0; i < kNClocks; ++i) {
    if (adc[i] > 1023) {
      AliWarning(Form("ADC (channel=%d) saturated: %d. Truncating to 1023",channel,adc[i]));
      adc[i] = 1023;
    }
  }
   
    for(Int_t i = 0; i < (kNClocks); ++i) {
      data =   (adc[i] & 0x3ff);
      data |= ((integrator & 0x1) << 10);      
      f->WriteBuffer((char*)&data,sizeof(data));
    }
    
}

//_____________________________________________________________________________
void AliADBuffer::WriteBeamFlags() {

//void AliADBuffer::WriteBeamFlags(Bool_t *bbFlag, Bool_t *bgFlag) {
  // The method writes information about
  // the Beam-Beam and Beam-Gas flags i.e. 
  // 10  shorts for the 4 channels 
  // of half a CIU card

  // Beam-beam and beam-gas flags are available
  // only for the triggered event-of-interest (sample index = 10)
  // As soon as trigger simulation would become more complex
  // and would allow to simulate neighbouring samples, this code
  // should be extended in order to fill all (or fraction) of the
  // flags
  /*/
  UShort_t data = 0;
  
  for(Int_t iEvOfInt = 0; iEvOfInt < kNClocks; iEvOfInt=iEvOfInt+2) {
        for(Int_t iChannel = 0; iChannel < 4; iChannel++) {
	  data = 0;
	  if (bbFlag[iChannel]) data |= ((bbFlag[iChannel] & 0x1) << 2*iChannel);
	  if (bgFlag[iChannel]) data |= ((bgFlag[iChannel] & 0x1) << 2*iChannel+1);
	  
	  if(iEvOfInt < (kNClocks - 1)) {      
	     if (bbFlag[iChannel]) data |= ((bbFlag[iChannel] & 0x1) << 8 + 2*iChannel);
	     if (bgFlag[iChannel]) data |= ((bgFlag[iChannel] & 0x1) << 8 + 2*iChannel+1);
	  }
	  f->WriteBuffer((char*)&data,sizeof(data));
        }
      }
   data = 0;
   f->WriteBuffer((char*)&data,sizeof(data));//Empty short in the end
   /*/
  for(Int_t i = 0; i < 12; i++) {
  	UShort_t data = 0;
  	f->WriteBuffer((char*)&data,sizeof(data));
	}
}

//_____________________________________________________________________________
void AliADBuffer::WriteMBInfo() {
  // The method writes information about
  // the 10 previous minimum-bias events
  // i.e. channels charge for each of these
  // 10 events (4*10 shorts for the 4 channels 
  // of half a CIU card)
    
  for(Int_t i = 0; i < 40; i++) {
    UShort_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }
}


//_____________________________________________________________________________
void AliADBuffer::WriteMBFlags() {
  // The method writes information about
  // the Minimum Bias flags
  // 5 16-bits words for the 4 channels 
  // of half a CIU card + one empty 16-bit


  for(Int_t i = 0; i < 6; i++) {
    UShort_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }
}

//_____________________________________________________________________________
void AliADBuffer::WriteBeamScalers() {
  // The method writes the AD beam scalers
  // For the moment there is no way to simulate
  // this, so we fill the necessary words with 0

  // Beam-beam and beam-gas scalers for
  // 4 individual channel (4x4 words)
  // (64-bit + 64-bit)*4 = 32bit * 16
  
  for(Int_t i = 0; i < 16; i++) {
    UInt_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }
}



//_____________________________________________________________________________
void AliADBuffer::WriteTiming(Float_t time, Float_t width) {
  // It writes the timing information into a raw data file. 
  // Being called by Digits2Raw

  // Writes the timing information
  UInt_t data = TMath::Nint(time) & 0xfff;
  data |= (TMath::Nint(width) & 0x7f) << 12;
  f->WriteBuffer((char*)&data,sizeof(data));
}
