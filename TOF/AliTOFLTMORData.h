#ifndef ALITOFLTMORDATA_H
#define ALITOFLTMORDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFLTMORData
{
 public:
  UInt_t GetORValue1() {return fORValue1;};
  UInt_t GetORValue2() {return fORValue2;};
  UInt_t GetORValue3() {return fORValue3;};
  UInt_t GetMBZ() {return fMBZ;};
 private:
  UInt_t fORValue1: 10;
  UInt_t fORValue2: 10;
  UInt_t fORValue3: 10;
  UInt_t fMBZ:       2;
};

#endif
