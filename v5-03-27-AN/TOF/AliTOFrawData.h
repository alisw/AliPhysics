#ifndef ALITOFRAWDATA_H
#define ALITOFRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////
//                                                //
//   This class provides the TOF raw data object  //
//                                                //
////////////////////////////////////////////////////

#include "TObject.h"

class AliTOFrawData : public TObject {
  // TOF rawData class
 public:
  AliTOFrawData(); // default ctr
  AliTOFrawData(Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g, Int_t h, Int_t l); // ctr
  AliTOFrawData(Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t ee, Int_t ff, Int_t g, Int_t h, Int_t l, Int_t deltaBC = 0, Int_t l0l1 = 0); // ctr
  ~AliTOFrawData() {}; // default dtr
  AliTOFrawData(const AliTOFrawData& r);     // dummy copy constructor
  AliTOFrawData& operator=(const AliTOFrawData& r); // dummy assignment operator
  void Update(Int_t tof, Int_t tot, Int_t leading, Int_t trailing, Int_t psBit, Int_t acq, Int_t errorFlag);

  Int_t GetTRM()        const {return fTRM;};
  Int_t GetTRMchain()   const {return fTRMchain;};
  Int_t GetTDC()        const {return fTDC;};
  Int_t GetTDCchannel() const {return fTDCchannel;};
  
  Int_t GetTOF() const {return fTime;};
  Int_t GetTOT() const;
  Int_t GetLeading() const {return fLeading;};
  Int_t GetTrailing() const {return fTrailing;};

  Int_t GetDeltaBC() const {return fDeltaBC;};
  Int_t GetL0L1Latency() const {return fL0L1Latency;};

  void SetDeltaBC(Int_t value) {fDeltaBC = value;};
  void SetL0L1Latency(Int_t value) {fL0L1Latency = value;};
  
 private:
  Int_t fACQflag;    // ACQ flag
  Int_t fPSbit;      // Packing bit
  
  Int_t fTRM;        // TRM ID
  Int_t fTRMchain;   // TRM Chain ID
  Int_t fTDC;        // TDC number 
  Int_t fTDCchannel; // TDC channel number
  
  Int_t fLeading;  // Leading Edge
  Int_t fTrailing; // Trailing Edge
  Int_t fToT;      // Time-Over-Threashould
  Int_t fTime;     // Time

  Int_t fError;      // Error flag
  
  Int_t fDeltaBC; // delta BC
  Int_t fL0L1Latency; // L0-L1 latency
  
  ClassDef(AliTOFrawData, 2)  // class for TOF raw data
};

#endif
