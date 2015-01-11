#ifndef ALITOFFEELIGHTCONFIG_H
#define ALITOFFEELIGHTCONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the TOFFEE light config objects.   //
//                                                           //
//   authors: Roberto Preghenella (R+)                       //
//   contacts: preghenella@bo.infn.it                        //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFFEEchannelConfig
{

 public:
  enum EStatus_t {
    kStatusEnabled = 0x1
  };
  
 private:
  UChar_t fStatus; // status
  Int_t fMatchingWindow; // matching window [ns]
  Int_t fLatencyWindow; // latency window [ns]

 public:
  AliTOFFEEchannelConfig() : fStatus(0x0), fMatchingWindow(0), fLatencyWindow(0) {}; // default construct
  ~AliTOFFEEchannelConfig() {}; // default destructor

  UChar_t GetStatus() const {return fStatus;}; // get status
  Int_t GetMatchingWindow() const {return fMatchingWindow;}; // get matching window
  Int_t GetLatencyWindow() const {return fLatencyWindow;}; // get latency window

  void SetStatus(UChar_t value) {fStatus = value;}; // set status
  void SetMatchingWindow(Int_t value) {fMatchingWindow = value;}; // set matching window
  void SetLatencyWindow(Int_t value) {fLatencyWindow = value;}; // set latency window

  Bool_t IsEnabled() const {return (GetStatus() & kStatusEnabled);}; // is enabled

};

//_____________________________________________________________________________

class AliTOFFEEtriggerConfig
{

 public:
  
 private:
  UInt_t fStatusMap; // status

 public:
  AliTOFFEEtriggerConfig() : fStatusMap(0x0) {}; // default construct
  ~AliTOFFEEtriggerConfig() {}; // default destructor

  UInt_t GetStatusMap() const {return fStatusMap;}; // get status map

  void SetStatusMap(UInt_t value) {fStatusMap = value;}; // set status map

};

//_____________________________________________________________________________

class AliTOFFEElightConfig
{

 private:
  static const Int_t fgkNumberOfChannels = 172800; // number of channels
  static const Int_t fgkNumberOfTriggerMaps = 72; // number of trigger maps
  Int_t fVersion; // version
  Int_t fRunNumber; // run number
  Int_t fRunType; // run type
  AliTOFFEEchannelConfig fChannelConfig[fgkNumberOfChannels]; // channel config array
  AliTOFFEEtriggerConfig fTriggerConfig[fgkNumberOfTriggerMaps]; // trigger config array

 public:
  AliTOFFEElightConfig() : fVersion(0), fRunNumber(0), fRunType(0), fChannelConfig() {}; // default construct
  ~AliTOFFEElightConfig() {}; // default destructor

  Int_t GetVersion() const {return fVersion;}; // get version
  Int_t GetRunNumber() const {return fRunNumber;}; // get run number
  Int_t GetRunType() const {return fRunType;}; // get run type
  AliTOFFEEchannelConfig *GetChannelConfig(Int_t i) {return (i < fgkNumberOfChannels ? &fChannelConfig[i] : NULL);}; // get channel config
  AliTOFFEEtriggerConfig *GetTriggerConfig(Int_t i) {return (i < fgkNumberOfTriggerMaps ? &fTriggerConfig[i] : NULL);}; // get trigger config

  void SetVersion(Int_t value) {fVersion = value;}; // get version
  void SetRunNumber(Int_t value) {fRunNumber = value;}; // get run number
  void SetRunType(Int_t value) {fRunType = value;}; // get run type

};



#endif /* ALITOFFEELIGHTCONFIG_H */
