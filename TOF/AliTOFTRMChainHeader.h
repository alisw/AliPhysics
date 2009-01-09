#ifndef ALITOFTRMCHAINHEADER_H
#define ALITOFTRMCHAINHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFTRMChainHeader
{
 public:
  UInt_t GetSlotID() {return fSlotID;};
  UInt_t GetBunchID() {return fBunchID;};
  UInt_t GetPB24Temp() {return fPB24Temp;};
  UInt_t GetPB24ID() {return fPB24ID;};
  UInt_t GetTSBit() {return fTSBit;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fSlotID:    4; // TRM number
  UInt_t fBunchID:  12; // bunch ID
  UInt_t fPB24Temp:  8; // PB24 temperature
  UInt_t fPB24ID:    3; // PB24 ID
  UInt_t fTSBit:     1; // TS bit
  UInt_t fWordType:  4; // word type
};

#endif
