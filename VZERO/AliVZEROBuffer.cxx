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
#include "AliLog.h"
#include "AliRawDataHeaderSim.h"
#include "AliVZEROBuffer.h"

//#include "TFile.h"
//#include "TTree.h"

ClassImp(AliVZEROBuffer)

//_____________________________________________________________________________
AliVZEROBuffer::AliVZEROBuffer():TObject(),
    fVerbose(0),
    f()
{
  //
  // default constructor
  //
}
//_____________________________________________________________________________
AliVZEROBuffer::AliVZEROBuffer(const char* fileName):TObject(),
    fVerbose(0),
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
AliVZEROBuffer::AliVZEROBuffer(const AliVZEROBuffer &source):TObject(source),
   fVerbose(0),
   f()

{
  // Copy Constructor
  this->fVerbose=source.fVerbose;
  return;
}

//_____________________________________________________________________________
AliVZEROBuffer& AliVZEROBuffer::operator=(const AliVZEROBuffer &source)

{
  //Assigment operator
  this->fVerbose=source.fVerbose;
  return *this;
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
void AliVZEROBuffer::WriteChannel(Int_t cell,Int_t ADC, Int_t Time){
  // It writes VZERO digits as a raw data file. 
  // Being called by AliVZERODDL.C

  UInt_t data = 0;
  // Information about previous 10 interaction
  // Not available in the simulation...
  for(Int_t i = 0; i < 5; i++)
    f->WriteBuffer((char*)&data,sizeof(data));

  // Now write the ADC charge for this channel
  if (ADC < 0 || ADC > 1023) {
    AliInfo(Form("ADC saturated: %d. Truncating to 1023",ADC));
    ADC = 1023;
  }
  data = ADC | 0x400; // 'Interaction' flag
  f->WriteBuffer((char*)&data,sizeof(data));

  data = 0;
  // Information about following 10 interaction
  // Not available in the simulation...
  for(Int_t i = 0; i < 5; i++)
    f->WriteBuffer((char*)&data,sizeof(data));

  // Now write the timing information
  data = Time & 0xfff;
  // The signal width is not available the digits!
  // To be added soon
  // data |= (width & 0x7f) << 12;
  f->WriteBuffer((char*)&data,sizeof(data));
}

//_____________________________________________________________________________
void AliVZEROBuffer::WriteScalers() {
  // The method writes the VZERO trigger scalers
  // For the moment there is no way to simulate
  // this, so we fill the necessary words with 0

  // First the general trigger scalers (16 of them)
  for(Int_t i = 0; i < 16; i++) {
    UInt_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }

  // Then beam-beam and beam-gas scalers for
  // each individual channel (4x64 words)
  for(Int_t i = 0; i < 256; i++) {
    UInt_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }
}

//_____________________________________________________________________________
void AliVZEROBuffer::WriteMBInfo() {
  // The method writes information about
  // the 10 previous minimum-bias events

  // First the bunch crossing numbers
  // for these 10 events
  for(Int_t i = 0; i < 10; i++) {
    UInt_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }

  // Then channels charge for each of these
  // 10 events (5 words/channel)
  for(Int_t i = 0; i < 320; i++) {
    UInt_t data = 0;
    f->WriteBuffer((char*)&data,sizeof(data));
  }
}
