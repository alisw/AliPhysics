#ifndef ALITOFTRMCHAINTRAILER_H
#define ALITOFTRMCHAINTRAILER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFTRMChainTrailer
{
 public:
  UInt_t GetStatus() {return fStatus;};
  UInt_t GetMBZ() {return fMBZ;};
  UInt_t GetEventCounter() {return fEventCounter;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fStatus:        4; // status
  UInt_t fMBZ:          12; // must-be-zero bits
  UInt_t fEventCounter: 12; // event counter
  UInt_t fWordType:      4; // word type
};

#endif
