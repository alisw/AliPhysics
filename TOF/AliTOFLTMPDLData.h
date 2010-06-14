#ifndef ALITOFLTMPDLDATA_H
#define ALITOFLTMPDLDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFLTMPDLData
{
 public:
  UInt_t GetPDLValue1() {return fPDLValue1;};
  UInt_t GetPDLValue2() {return fPDLValue2;};
  UInt_t GetPDLValue3() {return fPDLValue3;};
  UInt_t GetPDLValue4() {return fPDLValue4;};
 private:
  UInt_t fPDLValue1: 8;
  UInt_t fPDLValue2: 8;
  UInt_t fPDLValue3: 8;
  UInt_t fPDLValue4: 8;
};

#endif
