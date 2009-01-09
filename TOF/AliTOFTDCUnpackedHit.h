#ifndef ALITOFTDCUNPACHEDHIT_H
#define ALITOFTDCUNPACKEDHIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFRawDataFormat.h 23881 2008-02-12 16:46:22Z decaro $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFTDCUnpackedHit
{
 public:
  UInt_t GetHitTime() {return fHitTime;};
  UInt_t GetChan() {return fChan;};
  UInt_t GetTDCID() {return fTDCID;};
  UInt_t GetEBit() {return fEBit;};
  UInt_t GetPSBits() {return fPSBits;};
  UInt_t GetMBO() {return fMBO;};
 private:
  UInt_t fHitTime:  21; // leading or trailing edge measurement
  UInt_t fChan:      3; // TDC channel number
  UInt_t fTDCID:     4; // TDC ID
  UInt_t fEBit:      1; // E bit
  UInt_t fPSBits:    2; // PS bits
  UInt_t fMBO:       1; // must-be-zero bits
};

#endif
