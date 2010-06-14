#ifndef ALITOFDRMSTATUSHEADER2_H
#define ALITOFDRMSTATUSHEADER2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFDRMStatusHeader2
{
 public:
  UInt_t GetSlotID() {return fSlotID;};
  UInt_t GetSlotEnableMask() {return fSlotEnableMask;};
  UInt_t GetMBZ() {return fMBZ;};
  UInt_t GetFaultID() {return fFaultID;};
  UInt_t GetRTOBit() {return fRTOBit;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fSlotID:         4;
  UInt_t fSlotEnableMask: 11;
  UInt_t fMBZ:            1;
  UInt_t fFaultID:        11;
  UInt_t fRTOBit:        1;
  UInt_t fWordType:       4;
};

#endif
