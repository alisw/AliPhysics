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
  AliTOFrawData(Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t ee, Int_t ff, Int_t g, Int_t h, Int_t l); // ctr
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
  
  ClassDef(AliTOFrawData, 1)  // class for TOF raw data
};

#endif
