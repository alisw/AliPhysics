#ifndef ALITOFDRMSUMMARYDATA_H
#define ALITOFDRMSUMMARYDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides a summary for DRM data.             //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliTOFLTMSummaryData.h"
#include "AliTOFTRMSummaryData.h"

#define N_TRM 10
class AliTOFDRMSummaryData : public TObject
{
 public:
  AliTOFDRMSummaryData(); //default contructor
  AliTOFDRMSummaryData(const AliTOFDRMSummaryData &source); //copy contructor
  AliTOFDRMSummaryData &operator = (const AliTOFDRMSummaryData &source); //operator =
  virtual ~AliTOFDRMSummaryData(); //destructor
  /* getters */
  Bool_t   GetHeader() const {return fHeader;}; //get header
  Bool_t   GetTrailer() const {return fTrailer;}; //get trailer
  UShort_t GetSlotID() const {return fSlotID;}; //get slot ID
  UInt_t   GetEventWords() const {return fEventWords;}; //get event words
  UShort_t GetDRMID() const {return fDRMID;}; //get DRM ID
  UShort_t GetLocalEventCounter() const {return fLocalEventCounter;}; //get local event counter
  UShort_t GetPartecipatingSlotID() const {return fPartecipatingSlotID;}; //get partecipating slot ID
  UShort_t GetCBit() const {return fCBit;}; //get C bit
  UShort_t GetVersID() const {return fVersID;}; //get vers ID
  UShort_t GetDRMhSize() const {return fDRMhSize;}; //get DRM header size
  UShort_t GetSlotEnableMask() const {return fSlotEnableMask;}; //get slot enable mask
  UShort_t GetFaultID() const {return fFaultID;}; //get fault ID
  UShort_t GetRTOBit() const {return fRTOBit;}; //get RTO bit
  UShort_t GetL0BCID() const {return fL0BCID;}; //get L0 bunch ID
  UShort_t GetRunTimeInfo() const {return fRunTimeInfo;}; //get run time info
  UShort_t GetTemperature() const {return fTemperature;}; //get temperature
  UShort_t GetACKBit() const {return fACKBit;}; //get ACK bit
  UShort_t GetSensAD() const {return fSensAD;}; //get sens AD
  UInt_t   GetEventCRC() const {return fEventCRC;}; //get event CRC
  UInt_t   GetDecoderCRC() const {return fDecoderCRC;}; //get decoder CRC
  UShort_t GetDecoderSlotEnableMask() const {return fDecoderSlotEnableMask;}; //get decoder slot enable mask
  UShort_t GetDecoderSlotEnableMaskBit(UInt_t iBit) const {return fDecoderSlotEnableMask & (1 << iBit);}; //get decoder slot enable mask bit
  AliTOFLTMSummaryData *GetLTMSummaryData() const {return fLTMSummaryData;}; //get LTM summary data
  AliTOFTRMSummaryData *GetTRMSummaryData(Int_t TRM) const {return TRM < N_TRM ? fTRMSummaryData[TRM] : 0x0;}; //get TRM summary data
  /* setters */
  void SetHeader(Bool_t Header) {fHeader = Header;}; //set header
  void SetTrailer(Bool_t Trailer) {fTrailer = Trailer;}; //set trailer
  void SetSlotID(UShort_t SlotID) {fSlotID = SlotID;}; //set slot ID
  void SetEventWords(UInt_t EventWords) {fEventWords = EventWords;}; //set event words
  void SetDRMID(UShort_t DRMID) {fDRMID = DRMID;}; //set DRM ID
  void SetLocalEventCounter(UShort_t LocalEventCounter) {fLocalEventCounter = LocalEventCounter;}; //set local event counter
  void SetPartecipatingSlotID(UShort_t PartecipatingSlotID) {fPartecipatingSlotID = PartecipatingSlotID;}; //set partecipating slot ID
  void SetCBit(UShort_t CBit) {fCBit = CBit;}; //set C bit
  void SetVersID(UShort_t VersID) {fVersID = VersID;}; //set vers ID
  void SetDRMhSize(UShort_t DRMhSize) {fDRMhSize = DRMhSize;}; //set DRM header size
  void SetSlotEnableMask(UShort_t SlotEnableMask) {fSlotEnableMask = SlotEnableMask;}; //set slot enable mask
  void SetFaultID(UShort_t FaultID) {fFaultID = FaultID;}; //set fault ID
  void SetRTOBit(UShort_t RTOBit) {fRTOBit = RTOBit;}; //set RTO bit
  void SetL0BCID(UShort_t L0BCID) {fL0BCID = L0BCID;}; //set L0 bunch ID
  void SetRunTimeInfo(UShort_t RunTimeInfo) {fRunTimeInfo = RunTimeInfo;}; //set run time info
  void SetTemperature(UShort_t Temperature) {fTemperature = Temperature;}; //set temperature
  void SetACKBit(UShort_t ACKBit) {fACKBit = ACKBit;}; //set ACK bit
  void SetSensAD(UShort_t SensAD) {fSensAD = SensAD;}; //set sens ID
  void SetEventCRC(UInt_t EventCRC) {fEventCRC = EventCRC;}; //set event CRC
  void SetDecoderCRC(UInt_t DecoderCRC) {fDecoderCRC = DecoderCRC;}; //set decoder CRC
  void SetDecoderSlotEnableMask(UShort_t DecoderSlotEnableMask) {fDecoderSlotEnableMask = DecoderSlotEnableMask;}; //set decoder slot enable mask
  void SetDecoderSlotEnableMaskBit(UInt_t iBit) {fDecoderSlotEnableMask |= (1 << iBit);}; //get decoder slot enable mask
  /* methods */
  void Reset(); //reset
 private:
  Bool_t   fHeader; //header
  Bool_t   fTrailer; //trailer
  UShort_t fSlotID; //slot ID
  UInt_t   fEventWords; //event words
  UShort_t fDRMID; //DRM ID
  UShort_t fLocalEventCounter; //local event counter
  UShort_t fPartecipatingSlotID; //partecipating slot ID
  UShort_t fCBit; //C bit
  UShort_t fVersID; //vers ID
  UShort_t fDRMhSize; //DRM header size
  UShort_t fSlotEnableMask; //slot enable mask
  UShort_t fFaultID; //fault ID
  UShort_t fRTOBit; //RTO bit
  UShort_t fL0BCID;  //L0 bunch ID
  UShort_t fRunTimeInfo; //run time info
  UShort_t fTemperature; //temperature
  UShort_t fACKBit; //ACK bit
  UShort_t fSensAD; //sens ID
  UInt_t   fEventCRC; //event CRC
  UInt_t   fDecoderCRC; //decoder CRC
  UShort_t fDecoderSlotEnableMask; //decoder slot enable mask
  AliTOFLTMSummaryData *fLTMSummaryData; //LTM summary data
  AliTOFTRMSummaryData *fTRMSummaryData[N_TRM]; //TRM summary data

  ClassDef(AliTOFDRMSummaryData, 1);
};

#endif /* ALITOFDRMSUMMARYDATA_H */
