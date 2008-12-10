#ifndef ALITOFTRMDIAGNOSTICERRORWORD2_H
#define ALITOFTRMDIAGNOSTICERRORWORD2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFTRMDiagnosticErrorWord2
{
 public:
  UInt_t GetJtagErrorCode() {return fJtagErrorCode;};
  UInt_t GetTDCID() {return fTDCID;};
  UInt_t GetCBit() {return fCBit;};
  UInt_t GetMBZ() {return fMBZ;};
  UInt_t GetMBO() {return fMBO;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fJtagErrorCode: 11;
  UInt_t fTDCID:         4;
  UInt_t fCBit:          1;
  UInt_t fMBZ:           8;
  UInt_t fMBO:           4;
  UInt_t fWordType:      4;
};

#endif
