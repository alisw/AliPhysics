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
  AliTOFChannelOnline(Bool_t status, Float_t delay);
  AliTOFChannelOnline(const AliTOFChannelOnline &chan);
  AliTOFChannelOnline& operator= (const AliTOFChannelOnline &chan);
  virtual ~AliTOFChannelOnline(){};
  Bool_t   GetStatus()   	const {return fStatus;}
  Float_t  GetDelay()   	const {return fDelay;}
  void	   SetStatus(Bool_t status) {fStatus=status;}
  void     SetDelay(Float_t val) {fDelay=val;}

private:
  Bool_t   fStatus;     // Status of the channel (0:on 1:off)
  Float_t  fDelay;	// Delay
  //
  ClassDef(AliTOFChannelOnline,1)    // TOF Sensor Online Calibration object
};

#endif
