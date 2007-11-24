#ifndef ALITOFCHANNELONLINE_H
#define ALITOFCHANNELONLINE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for TOF Online calibration                 //
////////////////////////////////////////////////

#include "TObject.h"

class AliTOFChannelOnline: public TObject {

public:
  AliTOFChannelOnline();
  AliTOFChannelOnline(UChar_t status, Float_t delay);
  AliTOFChannelOnline(const AliTOFChannelOnline &chan);
  AliTOFChannelOnline& operator= (const AliTOFChannelOnline &chan);
  virtual ~AliTOFChannelOnline(){};
  UChar_t  GetStatus()   	const {return fStatus;}
  Float_t  GetDelay()   	const {return fDelay;}
  void	   SetStatus(UChar_t status) {fStatus=status;}
  void     SetDelay(Float_t val) {fDelay=val;}

  enum{
    kTOFOnlineUnknown=0x0, kTOFOnlineOk=0x2a, kTOFOnlineBad=0x15, 
      kTOFHWOk=0x1, kTOFHWBad=0x2, kTOFHVUnknown=0x0,
      kTOFPulserOk=0x4, kTOFPulserBad=0x8, kTOFPulserUnknown=0x0,
      kTOFNoiseOk=0x10, kTOFNoiseBad=0x20, kTOFNoiseUnknown=0x0
      };
  
  enum{
    kRightShiftHW=0,kRightShiftPulser=2,kRightShiftNoise=4
      };

private:
  UChar_t  fStatus;     // Status of the channel (0:on 1:off)
  Float_t  fDelay;	// Delay
  //
  ClassDef(AliTOFChannelOnline,2)    // TOF Sensor Online Calibration object
};

#endif
