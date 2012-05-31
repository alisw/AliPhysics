#ifndef ALITOFCHSIM_H
#define ALITOFCHSIM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for TOF calibration                 //
////////////////////////////////////////////////

#include "TObject.h"

class AliTOFChSim: public TObject {

public:
  AliTOFChSim();
  AliTOFChSim(const AliTOFChSim &chan);
  AliTOFChSim& operator= (const AliTOFChSim &chan);
  virtual ~AliTOFChSim(){};
  //
  Bool_t  IsSlewed() const {return fSlewedStatus;}
  void    SetSlewedStatus(Bool_t status) { fSlewedStatus = status;}
  void 	  SetSpectrum(Int_t ispectrum) {fSpectrum = ispectrum;}
  Int_t   GetSpectrum() const {return fSpectrum;}

private:
  Bool_t   fSlewedStatus;  // flag for decalibration status
  Int_t    fSpectrum;      // index of the spectrum used during decalibration

  ClassDef(AliTOFChSim,1)    // TOF Sensor Calibration object
};

#endif
