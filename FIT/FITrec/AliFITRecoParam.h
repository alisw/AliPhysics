#ifndef ALIFITRECOPARAM_H
#define ALIFITRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with FIT reconstruction parameters                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliDetectorRecoParam.h"

class AliFITRecoParam : public AliDetectorRecoParam
{
 public: 
  AliFITRecoParam();
  AliFITRecoParam(const AliFITRecoParam &p); //copy constructor
  AliFITRecoParam& operator=(const AliFITRecoParam &p);
  virtual ~AliFITRecoParam();
  
  static   AliFITRecoParam *GetLowFluxParam();        // make reco parameters for low  flux env
  static   AliFITRecoParam *GetHighFluxParam();       // make reco parameters for high flux env 
   
  //new staff
   Int_t  GetBadChannels(Int_t i) const {return fBadChannels[i];}
  void SetBadChannels(Int_t i,Int_t v) {fBadChannels[i] = v;}
 
  void PrintParameters() const;
  
 protected:
  Int_t  fBadChannels[288];   // bad channels map
  
  ClassDef(AliFITRecoParam, 1);
 
};
#endif
