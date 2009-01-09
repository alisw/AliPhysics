#ifndef ALITOFTRMTDCERROR_H
#define ALITOFTRMTDCERROR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFTRMTDCError
{
 public:
  UInt_t GetErrorFlags() {return fErrorFlags;};
  UInt_t GetMBZ() {return fMBZ;};
  UInt_t GetTDCID () {return fTDCID;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fErrorFlags: 15; // error flags
  UInt_t fMBZ:         9; // must-be-zero bits
  UInt_t fTDCID:       4; // TDC number
  UInt_t fWordType:    4; // word type
};

#endif
