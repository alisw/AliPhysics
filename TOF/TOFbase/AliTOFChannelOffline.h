#ifndef ALITOFCHANNELOFFLINE_H
#define ALITOFCHANNELOFFLINE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for TOF calibration                 //
////////////////////////////////////////////////

/* $Id$ */

#include "TObject.h"

class AliTOFChannelOffline: public TObject {

public:
  AliTOFChannelOffline();
  AliTOFChannelOffline(Float_t* slewingPar);
  AliTOFChannelOffline(const AliTOFChannelOffline &chan);
  AliTOFChannelOffline& operator= (const AliTOFChannelOffline &chan);
  virtual ~AliTOFChannelOffline();
  Float_t  GetSlewPar(Int_t index)	const {return fSlewPar[index];}
  Float_t* GetSlewPar()   const {return (float*)fSlewPar;}
  void     SetSlewPar(Int_t index, Float_t val) {fSlewPar[index]=val;}
  void 	   SetSlewPar(Float_t* SlewPar);

private:
  Float_t  fSlewPar[6];	// Time slewing parameters
  //
  ClassDef(AliTOFChannelOffline,1)    // TOF Sensor Calibration object
};

#endif
