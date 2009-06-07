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
  UChar_t GetStatus() const {return fStatus;}; // get status
  Int_t GetMatchingWindow() const {return fMatchingWindow;}; // get matching window
  Int_t GetLatencyWindow() const {return fLatencyWindow;}; // get latency window

  Bool_t IsEnabled() const {return (GetStatus() & kStatusEnabled);}; // is enabled

};

//_____________________________________________________________________________

class AliTOFFEElightConfig
{

 private:
  static const Int_t fgkNumberOfChannels = 157248; // number of channels
  Int_t fVersion; // version
  Int_t fRunNumber; // run number
  Int_t fRunType; // run type
  AliTOFFEEchannelConfig fChannelConfig[fgkNumberOfChannels]; // channel config array

 public:
  Int_t GetVersion() const {return fVersion;}; // get version
  Int_t GetRunNumber() const {return fRunNumber;}; // get run number
  Int_t GetRunType() const {return fRunType;}; // get run type
  AliTOFFEEchannelConfig *GetChannelConfig(Int_t i) {return i < fgkNumberOfChannels ? &fChannelConfig[i] : NULL;}; // get channel config

};



#endif /* ALITOFFEELIGHTCONFIG_H */
