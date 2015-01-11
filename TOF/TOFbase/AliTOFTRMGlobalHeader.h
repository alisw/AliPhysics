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
  UInt_t GetSlotID() const {return fSlotID;};
  UInt_t GetEventWords() const {return fEventWords;};
  UInt_t GetACQBits() const {return fACQBits;};
  UInt_t GetLBit() const {return fLBit;};
  UInt_t GetEBit() const {return fEBit;};
  UInt_t GetMBZ() const {return fMBZ;};
  UInt_t GetWordType() const {return fWordType;};
 private:
  UInt_t fSlotID:     4;
  UInt_t fEventWords: 13;
  UInt_t fACQBits:    2;
  UInt_t fLBit:       1;
  UInt_t fEBit:       1;
  UInt_t fMBZ:        7;
  UInt_t fWordType:   4;
};

#endif /* ALITOFTRMGLOBALHEADER_H */
