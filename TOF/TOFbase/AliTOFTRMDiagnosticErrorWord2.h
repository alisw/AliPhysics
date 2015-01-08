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
  UInt_t GetJtagErrorCode() const {return fJtagErrorCode;};
  UInt_t GetTDCID() const {return fTDCID;};
  UInt_t GetCBit() const {return fCBit;};
  UInt_t GetMBZ() const {return fMBZ;};
  UInt_t GetMBO() const {return fMBO;};
  UInt_t GetWordType() const {return fWordType;};
 private:
  UInt_t fJtagErrorCode: 11; // Jtag error code
  UInt_t fTDCID:          4; // TDC ID
  UInt_t fCBit:           1; // C bit
  UInt_t fMBZ:            8; // must-be-zero bits
  UInt_t fMBO:            4; // must-be-zero bits
  UInt_t fWordType:       4; // word type
};

#endif
