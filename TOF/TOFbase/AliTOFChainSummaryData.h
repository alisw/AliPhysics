#ifndef ALITOFCHAINSUMMARYDATA_H
#define ALITOFCHAINSUMMARYDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides a summary for TRM chain data.       //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliTOFTDCHitBuffer.h"
#include "AliTOFTDCErrorBuffer.h"

class AliTOFChainSummaryData : public TObject
{
 public:
  AliTOFChainSummaryData(); //default constructor
  AliTOFChainSummaryData(const AliTOFChainSummaryData &source); //copy constructor
  AliTOFChainSummaryData &operator = (const AliTOFChainSummaryData &source); //operator =
  virtual ~AliTOFChainSummaryData(); //destructor
  /* getters */
  Bool_t    GetHeader() const {return fHeader;}; //get header
  Bool_t    GetTrailer() const {return fTrailer;}; //get trailer
  UShort_t  GetChain() const {return fChain;}; //get chain
  UShort_t  GetBunchID() const {return fBunchID;}; //get bunch ID
  UShort_t  GetPB24Temp() const {return fPB24Temp;}; //get PB24 temp
  UShort_t  GetPB24ID() const {return fPB24ID;}; //get PB24 ID
  UShort_t  GetTSBit() const {return fTSBit;}; //get TS bit
  UShort_t  GetStatus() const {return fStatus;}; //get status
  UShort_t  GetEventCounter() const {return fEventCounter;}; //get event counter
  AliTOFTDCHitBuffer *GetTDCHitBuffer() const {return fTDCHitBuffer;}; //get TDC hit buffer
  AliTOFTDCHitBuffer *GetTDCPackedHitBuffer() const {return fTDCPackedHitBuffer;}; //get TDC packed hit buffer
  AliTOFTDCErrorBuffer *GetTDCErrorBuffer() const {return fTDCErrorBuffer;}; //get TDC error buffer
  /* setters */
  void SetHeader(Bool_t Header) {fHeader = Header;}; //set header
  void SetTrailer(Bool_t Trailer) {fTrailer = Trailer;}; //set trailer
  void SetChain(UShort_t Chain) {fChain = Chain;}; //set chain
  void SetBunchID(UShort_t BunchID) {fBunchID = BunchID;}; //set bunch ID
  void SetPB24Temp(UShort_t PB24Temp) {fPB24Temp = PB24Temp;}; //set PB24 temp
  void SetPB24ID(UShort_t PB24ID) {fPB24ID = PB24ID;}; //set PB24 ID
  void SetTSBit(UShort_t TSBit) {fTSBit = TSBit;}; //set TS bit
  void SetStatus(UShort_t Status) {fStatus = Status;}; //set status
  void SetEventCounter(UShort_t EventCounter) {fEventCounter = EventCounter;}; //set event counter
  /* methods */
  void Reset(); //reset 
 private:
  Bool_t    fHeader; //header detected
  Bool_t    fTrailer; //trailer detected
  UShort_t  fChain; //chain number
  UShort_t  fBunchID; //bunch ID
  UShort_t  fPB24Temp; //piggy-back temperature
  UShort_t  fPB24ID; //piggy-back ID
  UShort_t  fTSBit; //I2C reading of temperature sensor success
  UShort_t  fStatus; //status
  UShort_t  fEventCounter; //event counter
  AliTOFTDCHitBuffer   *fTDCHitBuffer; //TDC hit buffer
  AliTOFTDCHitBuffer   *fTDCPackedHitBuffer; //TDC packed hit buffer
  AliTOFTDCErrorBuffer *fTDCErrorBuffer; //TDC error buffer

  ClassDef(AliTOFChainSummaryData, 1);
};

#endif /* ALITOFCHAINSUMMARYDATA_H */
