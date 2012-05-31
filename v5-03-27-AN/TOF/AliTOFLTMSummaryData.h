#ifndef ALITOFLTMSUMMARYDATA_H
#define ALITOFLTMSUMMARYDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides a summary for LTM data.             //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"

#define LTM_N_PDL 48
#define LTM_N_ADC 60
#define LTM_N_OR  48

class AliTOFLTMSummaryData : public TObject
{
 public:
  AliTOFLTMSummaryData(); //default constructor
  AliTOFLTMSummaryData(const AliTOFLTMSummaryData &source); //copy contructor
  AliTOFLTMSummaryData &operator = (const AliTOFLTMSummaryData &source); //operator =
  virtual ~AliTOFLTMSummaryData(); //destructor
  /* getters */
  Bool_t  GetHeader() const {return fHeader;}; //get header
  Bool_t  GetTrailer() const {return fTrailer;}; //get trailer
  Short_t GetSlotID() {return fSlotID;}; //get slot ID
  Short_t GetEventWords() {return fEventWords;}; //get event words
  Short_t GetCBit() {return fCBit;}; //get C bit
  Short_t GetFault() {return fFault;}; //get fault
  Short_t GetPDL(Int_t i) {return fPDL[i];}; //get PDL
  Short_t GetADC(Int_t i) {return fADC[i];}; //get ADC
  Short_t GetOR(Int_t i) {return fOR[i];}; //get OR
  Short_t GetEventCRC() {return fEventCRC;}; //get event CRC
  Short_t GetEventNumber() {return fEventNumber;}; //set event number
  Short_t GetDecoderCRC() {return fDecoderCRC;}; //decoder CRC
  /* setters */
  void SetHeader(Bool_t Header) {fHeader = Header;}; //set header
  void SetTrailer(Bool_t Trailer) {fTrailer = Trailer;}; //set trailer
  void SetSlotID(Short_t SlotID) {fSlotID = SlotID;}; //set slot ID
  void SetEventWords(Short_t EventWords) {fEventWords = EventWords;}; //set event words
  void SetCBit(Short_t CBit) {fCBit = CBit;}; //set C bit
  void SetFault(Short_t Fault) {fFault = Fault;}; //set fault
  void SetPDL(Int_t i, Short_t PDL) {fPDL[i] = PDL;}; //set PDL
  void SetADC(Int_t i, Short_t ADC) {fADC[i] = ADC;}; //set ADC
  void SetOR(Int_t i, Short_t OR) {fOR[i] = OR;}; //set OR
  void SetEventCRC(Short_t EventCRC) {fEventCRC = EventCRC;}; //set event CRC
  void SetEventNumber(Short_t EventNumber) {fEventNumber = EventNumber;}; //set event number
  void SetDecoderCRC(Short_t DecoderCRC) {fDecoderCRC = DecoderCRC;}; //decoder CRC
  /* methods */
  void Reset(); //reset
 private:
  Bool_t  fHeader; //header
  Bool_t  fTrailer; //trailer
  Short_t fSlotID; //slot ID
  Short_t fEventWords; //event words
  Short_t fCBit; // C bit
  Short_t fFault; //fault
  Short_t fPDL[LTM_N_PDL]; //PDL setting
  Short_t fADC[LTM_N_ADC]; //ADC measurement
  Short_t fOR[LTM_N_OR]; //OR measurement
  Short_t fEventCRC; //event CRC
  Short_t fEventNumber; //event number
  Short_t fDecoderCRC; //decoder CRC

  ClassDef(AliTOFLTMSummaryData, 1);
};

#endif /* ALITOFLTMSUMMARYDATA_H */
