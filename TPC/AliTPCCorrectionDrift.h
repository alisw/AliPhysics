#ifndef ALITPCCORRECTIONDRIFT_H
#define ALITPCCORRECTIONDRIFT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTPCCorrectionDrift class                                                   //
// date: 02/05/2010                                                       //
// Authors: Maarian Ivanov, Jim Thomas, Magnus Mager, Stefan Rossegger                    //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"

class AliTPCCorrectionDrift : public AliTPCCorrection {
public:
  AliTPCCorrectionDrift();
  virtual ~AliTPCCorrectionDrift();
  // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);
  void Print(const Option_t* option) const;

public:
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);
  Double_t fZ0Aside;     // z- t0*vdrift shift A side
  Double_t fZ0Cside;     // z- t0*vdrift shift C side
  Double_t fVScale0;     // drift velocity scaling - constant
  Double_t fVScaleR;     // drift velocity scaling - radial
  Double_t fVScaleX;     // drift velocity scaling - global x
  Double_t fVScaleY;     // drift velocity scaling - global y
  //
  Double_t fIROCZ0;      // IROC to OROC shift due unknown reason (clusterer shift Ampl. dependents?)
  Double_t fOROCDZ;      // IROC to OROC slope shift due unknown reason (clusterer shift amplitude dependent?)
private:
  AliTPCCorrectionDrift(const AliTPCCorrectionDrift&);
  AliTPCCorrectionDrift &operator=(const AliTPCCorrectionDrift&);
  ClassDef(AliTPCCorrectionDrift,1);
};

#endif
