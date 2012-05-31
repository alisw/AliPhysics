#ifndef ALITOFChannel_H
#define ALITOFChannel_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for TOF calibration                 //
////////////////////////////////////////////////

#include "TObject.h"

class AliTOFChannel: public TObject {

public:
  AliTOFChannel();
  AliTOFChannel(Bool_t status, Float_t delay, Float_t* slewingPar);
  AliTOFChannel(const AliTOFChannel &chan);
  AliTOFChannel& operator= (const AliTOFChannel &chan);
  virtual ~AliTOFChannel(){};
  Bool_t   GetStatus()   	const {return fStatus;}
  Float_t  GetDelay()   	const {return fDelay;}
  Float_t  GetSlewPar(Int_t index)	const {return fSlewPar[index];}
  Float_t* GetSlewPar()   const {return (float*)fSlewPar;}
  void	   SetStatus(Bool_t status) {fStatus=status;}
  void     SetDelay(Float_t val) {fDelay=val;}
  void     SetSlewPar(Int_t index, Float_t val) {fSlewPar[index]=val;}
  void 	   SetSlewPar(Float_t* SlewPar);

private:
  Bool_t   fStatus;     // Status of the channel (0:on 1:off)
  Float_t  fDelay;	// Delay
  Float_t  fSlewPar[6];	// Time slewing parameters
  //
  ClassDef(AliTOFChannel,1)    // TOF Sensor Calibration object
};

#endif
