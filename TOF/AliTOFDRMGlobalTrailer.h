#ifndef ALITOFDRMGLOBALTRAILER_H
#define ALITOFDRMGLOBALTRAILER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFDRMGlobalTrailer
{
 public:
  UInt_t GetSlotID() {return fSlotID;};
  UInt_t GetLocalEventCounter() {return fLocalEventCounter;};
  UInt_t GetUNDEFINED() {return fUNDEFINED;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fSlotID:            4;
  UInt_t fLocalEventCounter: 12;
  UInt_t fUNDEFINED:         12;
  UInt_t fWordType:          4;
};

#endif
