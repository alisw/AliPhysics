#ifndef ALITOFPIDPARAMS_H
#define ALITOFPIDPARAMS_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliTODPIDparams
// class to store PID parameters for TOF in OADB
// Author: P. Antonioli, pietro.antonioli@to.infn.it
//***********************************************************

#include <TObject.h>
#include <TNamed.h>
#include "AliPIDResponse.h"

class AliTOFPIDParams : public TNamed {

 public:
  AliTOFPIDParams();
  AliTOFPIDParams(Char_t * name);
  virtual ~AliTOFPIDParams();

  enum {kSigPparams = 4};

  Float_t GetTOFresolution(void) const {return fTOFresolution;}
  AliPIDResponse::EStartTimeType_t GetStartTimeMethod(void) const {return fStartTime;}
  Float_t GetSigParams(Int_t i) const {
    return ((i >= 0)  && (i<kSigPparams)) ? fSigPparams[i] : 0;}    
  void SetTOFresolution(Float_t res){fTOFresolution = res;}
  void SetStartTimeMethod(AliPIDResponse::EStartTimeType_t method){fStartTime=method;}
  void SetSigPparams(Float_t *params);

 private:
  AliPIDResponse::EStartTimeType_t fStartTime;      // startTime method
  Float_t fTOFresolution;                           // TOF MRPC intrinsic resolution
  Float_t fSigPparams[kSigPparams];                // parameterisation of sigma(p) dependency 

  ClassDef(AliTOFPIDParams,1);

};

#endif

