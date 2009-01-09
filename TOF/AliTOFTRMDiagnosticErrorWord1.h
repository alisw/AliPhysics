#ifndef ALITOFTRMDIAGNOSTICERRORWORD1_H
#define ALITOFTRMDIAGNOSTICERRORWORD1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFTRMDiagnosticErrorWord1
{
 public:
  UInt_t GetFaultChipFlagID() {return fFaultChipFlagID;};
  UInt_t GetCBit() {return fCBit;};
  UInt_t GetMBZ() {return fMBZ;};
  UInt_t GetMBO() {return fMBO;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fFaultChipFlagID: 15; // fault chip flag ID
  UInt_t fCBit:             1; // C bit
  UInt_t fMBZ:              8; // must-be-zero bits
  UInt_t fMBO:              4; // must-be-zero bits
  UInt_t fWordType:         4; // word type
};

#endif
