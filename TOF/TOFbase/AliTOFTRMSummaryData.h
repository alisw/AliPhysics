#ifndef ALITOFTRMSUMMARYDATA_H
#define ALITOFTRMSUMMARYDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides a summary for TRM data.             //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliTOFChainSummaryData.h"

#define N_CHAIN 2
class AliTOFTRMSummaryData : public TObject
{
 public:
  AliTOFTRMSummaryData(); //default contructor
  AliTOFTRMSummaryData(const AliTOFTRMSummaryData &source); //copy contructor
  AliTOFTRMSummaryData &operator = (const AliTOFTRMSummaryData &source); //operator =
  virtual ~AliTOFTRMSummaryData(); //destructor
  /* getters */
  Bool_t   GetHeader() const {return fHeader;}; //get header
  Bool_t   GetTrailer() const {return fTrailer;}; //get trailer
  UShort_t GetSlotID() const {return fSlotID;}; //get slot ID
  UShort_t GetEventWords() const {return fEventWords;}; //get event words
  UShort_t GetACQBits() const {return fACQBits;}; //get ACQ bits
  UShort_t GetLBit() const {return fLBit;}; //get L bit
  UShort_t GetEBit() const {return fEBit;}; //get E bit
  UShort_t GetEventCRC() const {return fEventCRC;}; //get event CRC
  UShort_t GetEventCounter() const {return fEventCounter;}; //get event counter 
  UShort_t GetDecoderCRC() const {return fDecoderCRC;}; //get decoder CRC
  AliTOFChainSummaryData *GetChainSummaryData(Int_t Chain) const {return Chain < N_CHAIN ? fChainSummaryData[Chain] : 0x0;}; //get chain summary data
  /* setters */
  void SetHeader(Bool_t Header) {fHeader = Header;}; //set header
  void SetTrailer(Bool_t Trailer) {fTrailer = Trailer;}; //set trailer
  void SetSlotID(UShort_t SlotID) {fSlotID = SlotID;}; //set slot ID
  void SetEventWords(UShort_t EventWords) {fEventWords = EventWords;}; //set event words
  void SetACQBits(UShort_t ACQBits) {fACQBits = ACQBits;}; //set ACQ bits
  void SetLBit(UShort_t LBit) {fLBit = LBit;}; //set L bit
  void SetEBit(UShort_t EBit) {fEBit = EBit;}; //set E bit
  void SetEventCRC(UShort_t EventCRC) {fEventCRC = EventCRC;}; //set event CRC
  void SetEventCounter(UShort_t EventCounter) {fEventCounter = EventCounter;}; //set event counter
  void SetDecoderCRC(UShort_t DecoderCRC) {fDecoderCRC = DecoderCRC;}; //set decoder CRC
  /* methods */
  void Reset(); //reset
 private:
  Bool_t   fHeader; //header detected
  Bool_t   fTrailer; //trailer detected
  UShort_t fSlotID; //slot ID [3-12]
  UShort_t fEventWords; //number of TRM words
  UShort_t fACQBits; //HPTDC aquisition mode [0-3]
  UShort_t fLBit; //SEU detected inside LUT tables
  UShort_t fEBit; //empty event inserted while fixing SEU 
  UShort_t fEventCRC; //TRM computed CRC
  UShort_t fEventCounter; //event counter
  UShort_t fDecoderCRC; //decoder computed CRC
  AliTOFChainSummaryData *fChainSummaryData[N_CHAIN]; //chain summary data

  ClassDef(AliTOFTRMSummaryData, 1);
};

#endif /* ALITOFTRMSUMMARYDATA_H */
