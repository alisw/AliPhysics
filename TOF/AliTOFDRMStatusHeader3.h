#ifndef ALITOFDRMSTATUSHEADER3_H
#define ALITOFDRMSTATUSHEADER3_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFDRMStatusHeader3
{
 public:
  UInt_t GetSlot() {return fSlotID;};
  UInt_t GetL0BCID() {return fL0BCID;};
  UInt_t GetRunTimeInfo() {return fRunTimeInfo;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fSlotID: 4;
  UInt_t fL0BCID: 12;
  UInt_t fRunTimeInfo: 12;
  UInt_t fWordType: 4;
};

#endif
