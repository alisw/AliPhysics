#ifndef ALITOFRAWDATAFORMAT_H
#define ALITOFRAWDATAFORMAT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOF raw data bit fields.       //
//                                                           //
///////////////////////////////////////////////////////////////

//#include "TROOT.h"
#include "TObject.h"

//TRM global header
class AliTOFTRMGlobalHeader : public TObject
{
 public:
  UInt_t GetSlotID() {return fSlotID;};
  UInt_t GetEventWords() {return fEventWords;};
  UInt_t GetACQBits() {return fACQBits;};
  UInt_t GetLBit() {return fLBit;};
  UInt_t GetMBZ() {return fMBZ;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fSlotID:     4;
  UInt_t fEventWords: 13;
  UInt_t fACQBits:    2;
  UInt_t fLBit:       1;
  UInt_t fMBZ:        8;
  UInt_t fWordType:   4;
};

//TRM global trailer
class AliTOFTRMGlobalTrailer : public TObject
{
 public:
  UInt_t GetSlotID() {return fSlotID;};
  UInt_t GetEventCRC() {return fEventCRC;};
  UInt_t GetEventCounter() {return fEventCounter;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fSlotID:       4;
  UInt_t fEventCRC:     12;
  UInt_t fEventCounter: 12;
  UInt_t fWordType:     4;
};

//TRM chain header
class AliTOFTRMChainHeader : public TObject
{
 public:
  UInt_t GetSlotID() {return fSlotID;};
  UInt_t GetBunchID() {return fBunchID;};
  UInt_t GetPB24Temp() {return fPB24Temp;};
  UInt_t GetPB24ID() {return fPB24ID;};
  UInt_t GetTSBit() {return fTSBit;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fSlotID:   4;
  UInt_t fBunchID:  12;
  UInt_t fPB24Temp: 8;
  UInt_t fPB24ID:   3;
  UInt_t fTSBit:    1;
  UInt_t fWordType: 4;
};

//TRM chain trailer
class AliTOFTRMChainTrailer : public TObject
{
 public:
  UInt_t GetStatus() {return fStatus;};
  UInt_t GetMBZ() {return fMBZ;};
  UInt_t GetEventCounter() {return fEventCounter;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fStatus:       4;
  UInt_t fMBZ:          12;
  UInt_t fEventCounter: 12;
  UInt_t fWordType:     4;
};

//TDC packed hit
class AliTOFTDCPackedHit : public TObject
{
 public:
  UInt_t GetHitTime() {return fHitTime;};
  UInt_t GetTOTWidth() {return fTOTWidth;};
  UInt_t GetChan() {return fChan;};
  UInt_t GetTDCID() {return fTDCID;};
  UInt_t GetEBit() {return fEBit;};
  UInt_t GetPSBits() {return fPSBits;};
  UInt_t GetMBO() {return fMBO;};
 private:
  UInt_t fHitTime:  13;
  UInt_t fTOTWidth: 8;
  UInt_t fChan:     3;
  UInt_t fTDCID:    4;
  UInt_t fEBit:     1;
  UInt_t fPSBits:   2;
  UInt_t fMBO:      1;
};

//TDC unpacked hit
class AliTOFTDCUnpackedHit : public TObject
{
 public:
  UInt_t GetHitTime() {return fHitTime;};
  UInt_t GetChan() {return fChan;};
  UInt_t GetTDCID() {return fTDCID;};
  UInt_t GetEBit() {return fEBit;};
  UInt_t GetPSBits() {return fPSBits;};
  UInt_t GetMBO() {return fMBO;};
 private:
  UInt_t fHitTime:  21;
  UInt_t fChan:     3;
  UInt_t fTDCID:    4;
  UInt_t fEBit:     1;
  UInt_t fPSBits:   2;
  UInt_t fMBO:      1;
};

//TRM TDC error
class AliTOFTRMTDCError : public TObject
{
 public:
  UInt_t GetErrorFlags() {return fErrorFlags;};
  UInt_t GetMBZ() {return fMBZ;};
  UInt_t GetTDCID () {return fTDCID;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fErrorFlags: 15;
  UInt_t fMBZ:        9;
  UInt_t fTDCID:      4;
  UInt_t fWordType:   4;
};

//TRM diagnostic error word 1
class AliTOFTRMDiagnosticErrorWord1 : public TObject
{
 public:
  UInt_t GetFaultChipFlagID() {return fFaultChipFlagID;};
  UInt_t GetCBit() {return fCBit;};
  UInt_t GetMBZ() {return fMBZ;};
  UInt_t GetMBO() {return fMBO;};
  UInt_t GetWordType() {return fWordType;};
 private:
  UInt_t fFaultChipFlagID: 15;
  UInt_t fCBit:            1;
  UInt_t fMBZ:             8;
  UInt_t fMBO:             4;
  UInt_t fWordType:        4;
};

//TRM diagnostic error word 2
class AliTOFTRMDiagnosticErrorWord2 : public TObject
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

#endif /* ALITOFRAWDATAFORMAT_H */
