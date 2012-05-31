#ifndef ALITOFCHANNELTASK_H
#define ALITOFCHANNELTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for TOF calibration                 //
////////////////////////////////////////////////

/* $Id$ */

#include "TObject.h"

class AliTOFChannelTask: public TObject {

public:
  AliTOFChannelTask();
  AliTOFChannelTask(Float_t* slewingPar);
  AliTOFChannelTask(const AliTOFChannelTask &chan);
  AliTOFChannelTask& operator= (const AliTOFChannelTask &chan);
  virtual ~AliTOFChannelTask(){};
  Float_t  GetSlewPar(Int_t index)	const {return fSlewPar[index];}
  Float_t* GetSlewPar()   const {return (float*)fSlewPar;}
  void     SetSlewPar(Int_t index, Float_t val) {fSlewPar[index]=val;}
  void 	   SetSlewPar(Float_t* SlewPar);

private:
  Float_t  fSlewPar[6];	// Time slewing parameters
  //
  ClassDef(AliTOFChannelTask,1)    // TOF Sensor Calibration object
};

#endif
