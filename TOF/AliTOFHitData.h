#ifndef ALITOFHITDATA_H
#define ALITOFHITDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides the key-reading for TOF raw data.   //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"

class AliTOFHitData : public TObject{
 public:
  AliTOFHitData();
  ~AliTOFHitData()
    {};
  
  
  AliTOFHitData(const AliTOFHitData &source);
  
  AliTOFHitData& operator=(const AliTOFHitData & source); //ass. operator
  
  Int_t  *GetVolume()  {return fVolume;};
  Short_t GetDDLID() const {return fDDLID;};
  Short_t GetSlotID() const {return fSlotID;};
  Short_t GetACQ() const {return fACQ;};
  Short_t GetChain() const {return fChain;};
  Short_t GetPS() const {return fPS;};
  Short_t GetTDC() const {return fTDC;};
  Short_t GetChan() const {return fChan;};
  Float_t GetTime() const {return fTime;};
  Int_t GetTimeBin() const {return fTimeBin;};
  Float_t GetTOT() const {return fTOT;};
  Int_t GetTOTBin() const {return fTOTBin;};
  Int_t GetDeltaBunchID() const {return fDeltaBunchID;};
  Int_t GetDeltaEventCounter() const {return fDeltaEventCounter;};

  void SetVolume(Int_t *Volume);

  void SetDDLID(Short_t DDLID)    { fDDLID=DDLID;};
  void SetSlotID(Short_t slotID)  { fSlotID=slotID;};
  void SetACQ(Short_t ACQ)        { fACQ=ACQ;};
  void SetChain(Short_t chain)    { fChain=chain;};
  void SetPS(Short_t PS)          { fPS=PS;};
  void SetTDC(Short_t TDC)        { fTDC=TDC;};
  void SetChan(Short_t chan)      { fChan=chan;};
  void SetTime(Float_t time)      { fTime=time;};
  void SetTimeBin(Int_t timeBin) {fTimeBin=timeBin;};
  void SetTOT(Float_t TOT)        { fTOT=TOT;};
  void SetTOTBin(Int_t TOTBin) {fTOTBin=TOTBin;};
  void SetDeltaBunchID(Int_t Value) {fDeltaBunchID=Value;};
  void SetDeltaEventCounter(Int_t Value) {fDeltaEventCounter=Value;};
  
 private:
  Int_t   fVolume[5];  // TOF volume index
  Short_t fDDLID;      // DDL index
  Short_t fSlotID;     // slot index
  Short_t fACQ;        // ACQ flag
  Short_t fChain;      // chain index
  Short_t fPS;         // PS bit
  Short_t fTDC;        // TDC index
  Short_t fChan;       // channel index
  Float_t fTime;      // time [ns]
  Int_t fTimeBin;      // time [TDC bin = 24.4ps]
  Float_t fTOT;       // tot [ns]
  Int_t fTOTBin;       // TOT [TOT bin = 48.4ps]
  Int_t fDeltaBunchID; // TRM bunchID - miniEventID
  Int_t fDeltaEventCounter; // TRM event counter - DRM local event counter

  ClassDef(AliTOFHitData, 1);
};

#endif
