#ifndef ALITOFDECODERSUMMARYDATA_H
#define ALITOFDECODERSUMMARYDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides a summary for decoder data.         //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliTOFDRMSummaryData.h"

//decoder summary data
class AliTOFDecoderSummaryData : public TObject
{
 public:
  AliTOFDecoderSummaryData(); //default contructor
  AliTOFDecoderSummaryData(const AliTOFDecoderSummaryData &source); //copy contructor
  AliTOFDecoderSummaryData &operator = (const AliTOFDecoderSummaryData &source); //operator =
  virtual ~AliTOFDecoderSummaryData(); //destructor
  /* getters */
  UInt_t   GetRunNumber() const {return fRunNumber;}; //get run number
  UInt_t   GetEventNumber() const {return fEventNumber;}; //get event number
  UInt_t   GetEquipmentID() const {return fEquipmentID;}; //get equipment ID
  UInt_t   GetInputWords() const {return fInputWords;}; //get input words
  UInt_t   GetDecodedWords() const {return fDecodedWords;}; //get decoded words
  UShort_t GetDecoderStatus() const {return fDecoderStatus;}; //get decoder status
  Bool_t   GetErrorDetected() const {return fErrorDetected;}; //get error detected
  UShort_t GetErrorSlotID() const {return fErrorSlotID;}; // get error slot id
  UShort_t GetCurrentDRMID() const {return fCurrentDRMID;}; //get current DRM ID
  UShort_t GetCurrentSlotID() const {return fCurrentSlotID;}; //get current slot ID
  UShort_t GetCurrentChain() const {return fCurrentChain;}; //get current chain
  Bool_t   GetV2718Patch() const {return fV2718Patch;}; //get V2718 patch
  Bool_t   GetRecoverError() const {return fRecoverError;}; //get recover error
  Bool_t   GetRecoveringError() const {return fRecoveringError;}; //get recovering error
  Bool_t   GetSpider() const {return fSpider;}; //get spider flag
  AliTOFDRMSummaryData *GetDRMSummaryData() const {return fDRMSummaryData;}; //get DRM summary data
  /* setters */
  void SetRunNumber(UInt_t RunNumber) {fRunNumber = RunNumber;}; //set run number
  void SetEventNumber(UInt_t EventNumber) {fEventNumber = EventNumber;}; //set event number
  void SetEquipmentID(UInt_t EquipmentID) {fEquipmentID = EquipmentID;}; //set equipment ID
  void SetInputWords(UInt_t InputWords) {fInputWords = InputWords;}; //set input words
  void SetDecodedWords(UInt_t DecodedWords) {fDecodedWords = DecodedWords;}; //set decoded words
  void SetDecoderStatus(UShort_t DecoderStatus) {fDecoderStatus = DecoderStatus;}; //set decoder status
  void SetErrorDetected(Bool_t ErrorDetected) {fErrorDetected = ErrorDetected;}; //set error detected
  void SetErrorSlotID(UShort_t value) {fErrorSlotID |= 1<< value;}; // set error slot id
  void SetCurrentDRMID(UShort_t CurrentDRMID) {fCurrentDRMID = CurrentDRMID;}; //set current DRM ID
  void SetCurrentSlotID(UShort_t CurrentSlotID) {fCurrentSlotID = CurrentSlotID;}; //set current slot ID
  void SetCurrentChain(UShort_t CurrentChain) {fCurrentChain = CurrentChain;}; //set current chain
  void SetV2718Patch(Bool_t V2718Patch) {fV2718Patch = V2718Patch;}; //set V2718 patch
  void SetRecoverError(Bool_t RecoverError) {fRecoverError = RecoverError;}; //set recover error
  void SetRecoveringError(Bool_t RecoveringError) {fRecoveringError = RecoveringError;}; //set recovering error
  void SetSpider(Bool_t Spider) {fSpider = Spider;}; //set spider flag
  /* methods */
  void Reset(); //reset
 private:
  UInt_t   fRunNumber; //run number
  UInt_t   fEventNumber; //event number
  UInt_t   fEquipmentID; //equipment ID
  UInt_t   fInputWords; //input words
  UInt_t   fDecodedWords; //decoded words
  UShort_t fDecoderStatus; //decoder status
  Bool_t   fErrorDetected; //error detected
  UShort_t   fErrorSlotID; //error slot ID
  UShort_t fCurrentDRMID; //current DRM ID
  UShort_t fCurrentSlotID; //current slot ID
  UShort_t fCurrentChain; //current chain
  Bool_t   fV2718Patch; //V2718 patch
  Bool_t   fRecoverError; //recover error
  Bool_t   fRecoveringError; //recoverisng error
  Bool_t   fSpider; //spider flag
  AliTOFDRMSummaryData *fDRMSummaryData; //DRM summary data

  ClassDef(AliTOFDecoderSummaryData, 1);
};

#endif /* ALITOFDECODERSUMMARYDATA */
