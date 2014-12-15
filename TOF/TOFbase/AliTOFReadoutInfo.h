#ifndef ALITOFREADOUTINFO_H
#define ALITOFREADOUTINFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

// *
// *
// *
// * this class defines the TOF object to be stored
// * in Reference data a run-by-run basis in order to have
// * info about readout electronics
// * 
// *
// *
// *

#include "TObject.h"

class TH1F;
class TH2F;

class AliTOFReadoutInfo :
public TObject
{

 public:

  AliTOFReadoutInfo(); // default constructor
  virtual ~AliTOFReadoutInfo(); // default destructor
  AliTOFReadoutInfo(const AliTOFReadoutInfo &source); // copy constructor
  AliTOFReadoutInfo &operator=(const AliTOFReadoutInfo &source); // operator=

  TH1F *GetChainEfficiency() const {return fChainEfficiency;}; // getter
  TH1F *GetTRMData() const {return fTRMData;}; // getter
  TH1F *GetTRMEmptyEvent() const {return fTRMEmptyEvent;}; // getter
  TH1F *GetTRMBadEventCounter() const {return fTRMBadEventCounter;}; // getter
  TH1F *GetTRMBadCRC() const {return fTRMBadCRC;}; // getter
  TH1F *GetChainData() const {return fChainData;}; // getter
  TH1F *GetChainBadStatus() const {return fChainBadStatus;}; // getter
  TH1F *GetChainBadEventCounter() const {return fChainBadEventCounter;}; // getter
  TH1F *GetTDCError() const {return fTDCError;}; // getter
  TH2F *GetTDCErrorFlags() const {return fTDCErrorFlags;}; // getter


  void SetChainEfficiency(TH1F *value) {fChainEfficiency = value;}; // getter
  void SetTRMData(TH1F *value) {fTRMData = value;}; // getter
  void SetTRMEmptyEvent(TH1F *value) {fTRMEmptyEvent = value;}; // getter
  void SetTRMBadEventCounter(TH1F *value) {fTRMBadEventCounter = value;}; // getter
  void SetTRMBadCRC(TH1F *value) {fTRMBadCRC = value;}; // getter
  void SetChainData(TH1F *value) {fChainData = value;}; // getter
  void SetChainBadStatus(TH1F *value) {fChainBadStatus = value;}; // getter
  void SetChainBadEventCounter(TH1F *value) {fChainBadEventCounter = value;}; // getter
  void SetTDCError(TH1F *value) {fTDCError = value;}; // getter
  void SetTDCErrorFlags(TH2F *value) {fTDCErrorFlags = value;}; // getter


 private:

  TH1F *fChainEfficiency; // chain efficiency
  TH1F *fTRMData; // TRM data
  TH1F *fTRMEmptyEvent; // TRM empty event
  TH1F *fTRMBadEventCounter; // TRM bad event counter
  TH1F *fTRMBadCRC; // TRM bad CRC
  TH1F *fChainData; // chain data
  TH1F *fChainBadStatus; // chain bad status
  TH1F *fChainBadEventCounter; // chain bad event counter
  TH1F *fTDCError; // TDC error
  TH2F *fTDCErrorFlags; // TDC error flags

  ClassDef(AliTOFReadoutInfo, 1);
};

#endif /* ALITOFREADOUTINFO_H */
