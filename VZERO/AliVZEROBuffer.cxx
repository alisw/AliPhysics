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
// To be used in Alice Data Challenges
// This class is used by AliVZERODDL.C macro
// Author: B. Cheynis

#include <Riostream.h>
#include <TObjArray.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliRawDataHeaderSim.h"
#include "AliVZEROBuffer.h"
#include "AliVZEROdigit.h"

ClassImp(AliVZEROBuffer)

//_____________________________________________________________________________
AliVZEROBuffer::AliVZEROBuffer():TObject(),
    fRemainingWord(0),
    f()
{
  //
  // default constructor
  //
}
//_____________________________________________________________________________
AliVZEROBuffer::AliVZEROBuffer(const char* fileName):TObject(),
    fRemainingWord(0),
    f()
{
  // Constructor
  f = new AliFstream(fileName);
  // fout=new TFile(fileName,"recreate");
  // tree=new TTree("tree","Values");
  AliRawDataHeaderSim header;
  f->WriteBuffer((char*)(&header), sizeof(header));

}

//_____________________________________________________________________________
AliVZEROBuffer::~AliVZEROBuffer(){
  // Destructor, it closes the IO stream
  AliRawDataHeaderSim header;
  header.fSize = f->Tellp();
  header.SetAttribute(0);  // valid data
  f->Seekp(0);
  f->WriteBuffer((char*)(&header), sizeof(header));
  delete f;
  //delete tree;
  //delete fout;
}

//_____________________________________________________________________________
void AliVZEROBuffer::WriteTriggerInfo(UInt_t trigger) {
  // The method writes VZERO trigger information
  // This info is contained in the first two
  // raw-data words following the raw-data header (CDH).

  f->WriteBuffer((char*)(&trigger),sizeof(trigger));

  // By default all the inputs are unmasked... Hopefully
  UInt_t triggerMask = 0xffff;
  f->WriteBuffer((char*)(&triggerMask),sizeof(triggerMask));
}

//_____________________________________________________________________________
void AliVZEROBuffer::WriteTriggerScalers() {
  // The method writes the VZERO trigger scalers
  // For the moment there is no way to simulate
  // this, so we fill the necessary 16 words with 0

  // First the general trigger scalers (16 of them)
  for(Int_t i = 0; i < 16; i++) {
      UInt_t data = 0;
      f->WriteBuffer((char*)&data,sizeof(data));
  }
}

//_____________________________________________________________________________
void AliVZEROBuffer::WriteBunchNumbers() {
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
void AliVZEROBuffer::WriteChannel(Int_t channel, Short_t *adc, Bool_t integrator){
  // It writes VZERO charge information into a raw data file. 
  // Being called by Digits2Raw
  
  UInt_t data = 0;
  for(Int_t i = 0; i < AliVZEROdigit::kNClocks; ++i) {
    if (adc[i] > 1023) {
      AliWarning(Form("ADC (channel=%d) saturated: %d. Truncating to 1023",channel,adc[i]));
      adc[i] = 1023;
    }
  }
  
  if(channel%2 == 0) {
    for(Int_t i = 0; i < (AliVZEROdigit::kNClocks/2); ++i) {
      data =   (adc[2*i] & 0x3ff);
      data |= ((integrator & 0x1) << 10);

      data |= ((adc[2*i+1] & 0x3ff) << 16);
      data |= ((!integrator & 0x1) << 26);

      f->WriteBuffer((char*)&data,sizeof(data));
    }
    fRemainingWord = (adc[AliVZEROdigit::kNClocks-1] & 0x3ff);
    fRemainingWord |= ((integrator & 0x1) << 10);
  }
  else {
    data = fRemainingWord;
    data |= ((adc[0] & 0x3ff) << 16);
    data |= ((integrator & 0x1) << 26);
    f->WriteBuffer((char*)&data,sizeof(data));

    for(Int_t i = 1; i <= (AliVZEROdigit::kNClocks/2); ++i) {
      data =   (adc[2*i-1] & 0x3ff);
      data |= ((!integrator & 0x1) << 10);

      data |= ((adc[2*i] & 0x3ff) << 16);
      data |= ((integrator & 0x1) << 26);

      f->WriteBuffer((char*)&data,sizeof(data));
    }
  }
    
}

//_____________________________________________________________________________
void AliVZEROBuffer::WriteBeamFlags(Bool_t *bbFlag, Bool_t *bgFlag) {
  // The method writes information about
  // the Beam-Beam and Beam-Gas flags i.e. 
  // 6  words for the 4 channels 
  // of half a CIU card

  // Beam-beam and beam-gas flags are available
  // only for the triggered event-of-interest (sample index = 10)
  // As soon as trigger simulation would become more complex
  // and would allow to simulate neighbouring samples, this code
  // should be extended in order to fill all (or fraction) of the
  // flags
  for(Int_t i = 0; i < 2; i++) {
    UInt_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }
  {
    UInt_t data = 0;
    for(Int_t iChannel = 0; iChannel < 4; ++iChannel) {
      if (bbFlag[iChannel]) data |= (1 << (2*iChannel + 16));
      if (bgFlag[iChannel]) data |= (1 << (2*iChannel + 17));
    }
    f->WriteBuffer((char*)&data,sizeof(data));
  }  
  for(Int_t i = 0; i < 3; i++) {
    UInt_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }

}


//_____________________________________________________________________________
void AliVZEROBuffer::WriteMBInfo() {
  // The method writes information about
  // the 10 previous minimum-bias events
  // i.e. channels charge for each of these
  // 10 events (20 words for the 4 channels 
  // of half a CIU card)
    
  for(Int_t i = 0; i < 20; i++) {
    UInt_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }
}


//_____________________________________________________________________________
void AliVZEROBuffer::WriteMBFlags() {
  // The method writes information about
  // the Minimum Bias flags
  // 3 32-bits words for the 4 channels 
  // of half a CIU card


  for(Int_t i = 0; i < 3; i++) {
    UInt_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }
}

//_____________________________________________________________________________
void AliVZEROBuffer::WriteBeamScalers() {
  // The method writes the VZERO beam scalers
  // For the moment there is no way to simulate
  // this, so we fill the necessary words with 0

  // Beam-beam and beam-gas scalers for
  // 4 individual channel (4x4 words)
  
  for(Int_t i = 0; i < 16; i++) {
    UInt_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }
}

//_____________________________________________________________________________
void AliVZEROBuffer::WriteTiming(Float_t time, Float_t width) {
  // It writes the timing information into a raw data file. 
  // Being called by Digits2Raw

  // Writes the timing information
  UInt_t data = TMath::Nint(time/(25.0/256.0)) & 0xfff;
  data |= (TMath::Nint(width/(25./64.)) & 0x7f) << 12;
  f->WriteBuffer((char*)&data,sizeof(data));
}
