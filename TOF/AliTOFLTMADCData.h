#ifndef ALITOFLTMADCDATA_H
#define ALITOFLTMADCDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFLTMADCData
{
 public:
  UInt_t GetADCValue1() {return fADCValue1;};
  UInt_t GetADCValue2() {return fADCValue2;};
  UInt_t GetADCValue3() {return fADCValue3;};
  UInt_t GetMBZ() {return fMBZ;};
 private:
  UInt_t fADCValue1: 10;
  UInt_t fADCValue2: 10;
  UInt_t fADCValue3: 10;
  UInt_t fMBZ:       2;
};

#endif
