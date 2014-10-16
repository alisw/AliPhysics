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
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  From ACORDE digits to Raw data
//
// there are 4 words of 32 bits corresponding to word 9 to 12
// (words up to 8 correspond to the header)
// Word 9: bits 1 to 30 --> Modules 1 to 30
//         bits 31-32 = '00'
// Word 10: bits 1 to 30 --> Modules 31 to 60
//          bits 31-32 = '01'
// Word 11: bits 1 to 30 --> Modules 1 to 30
//          bits 31-32 = '10'
// Word 12: bits 1 to 30 --> Modules 1 to 30
//          bits 31-32 = '11'
// Words 9 and 10 are the single muon trigger
// Words 11 and 12 are the multi muon trigger
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliACORDERawData.h"
#include "AliDAQ.h"
#include "AliFstream.h"
#include "AliRawDataHeaderSim.h"


ClassImp(AliACORDERawData)


AliACORDERawData::AliACORDERawData()
  :TObject(),
   fWord9(0),
   fWord10(0),
   fWord11(0),
   fWord12(0)
{
}

AliACORDERawData::AliACORDERawData(const AliACORDERawData &r)
  :TObject(),
   fWord9(0),
   fWord10(0),
   fWord11(0),
   fWord12(0)
{
  ((AliACORDERawData &) r).Copy(*this);
}

AliACORDERawData::~AliACORDERawData()

{

}

AliACORDERawData &AliACORDERawData::operator=(const AliACORDERawData &r)

{
  if (this != &r)  ((AliACORDERawData &) r).Copy(*this);
  return *this;
}

void AliACORDERawData::WriteACORDERawData(Bool_t *b,Bool_t multi)

{
  // set words
  SetACORDERawWords(b,multi);

  // open output file
  const char *fileName = AliDAQ::DdlFileName("ACORDE",0);
  AliFstream* fFile = new AliFstream(fileName);

  // write header
  AliRawDataHeaderSim header;
  UInt_t header_position = fFile->Tellp();
  fFile->WriteBuffer((char*)(&header), sizeof(header));

  // write digits
  fFile->WriteBuffer((char*)(&fWord9), sizeof(fWord9));
  fFile->WriteBuffer((char*)(&fWord10), sizeof(fWord10));
  fFile->WriteBuffer((char*)(&fWord11), sizeof(fWord11));
  fFile->WriteBuffer((char*)(&fWord12), sizeof(fWord12));
  
  // write header again
  UInt_t current_position = fFile->Tellp();
  fFile->Seekp(header_position);
  header.fSize = current_position-header_position;
  header.SetAttribute(0);  // valid data
  fFile->WriteBuffer((char*)(&header), sizeof(header));
  fFile->Seekp(current_position);
}

void AliACORDERawData::SetACORDERawWords(Bool_t *b,Bool_t multi)

{
  // set modules
  for (Int_t i=0;i<30;i++) {
    if (b[i]) {
      fWord9|=(1<<i);
      if (multi) fWord11|=(1<<i);
    }
    if (b[i+30]) {
      fWord10|=(1<<i);
      if (multi) fWord12|=(1<<i);
    }
  } // end for
  // set labels
  fWord10|=(unsigned int)(1<<30);
  fWord12|=(unsigned int)(1<<30);
  fWord11|=(unsigned int)(1<<31); 
  fWord12|=(unsigned int)(1<<31);
}
