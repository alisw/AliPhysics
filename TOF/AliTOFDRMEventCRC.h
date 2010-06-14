#ifndef ALITOFDRMEVENTCRC_H
#define ALITOFDRMEVENTCRC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFDRMEventCRC
{
 public:
  UInt_t GetSlotID() {return fSlotID;};
  UInt_t GetEventCRC() {return fEventCRC;};
  UInt_t GetUNDEFINED() {return fUNDEFINED;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fSlotID:    4;
  UInt_t fEventCRC:  16;
  UInt_t fUNDEFINED: 12;
  UInt_t fWordType:  4;
};

#endif
