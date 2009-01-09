#ifndef ALITOFTRMGLOBALHEADER_H
#define ALITOFTRMGLOBALHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFTRMGlobalHeader
{
 public:
  UInt_t GetSlotID() {return fSlotID;};
  UInt_t GetEventWords() {return fEventWords;};
  UInt_t GetACQBits() {return fACQBits;};
  UInt_t GetLBit() {return fLBit;};
  UInt_t GetMBZ() {return fMBZ;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fSlotID:      4; // TRM number
  UInt_t fEventWords: 13; // event word
  UInt_t fACQBits:     2; // ACQ bits
  UInt_t fLBit:        1; // L bit
  UInt_t fMBZ:         8; // must-be-zero bits
  UInt_t fWordType:    4; // word type
};

#endif /* ALITOFTRMGLOBALHEADER_H */
