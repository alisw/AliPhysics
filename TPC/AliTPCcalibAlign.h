#ifndef ALITPCCALIBALIGN_H
#define ALITPCCALIBALIGN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

#include "TObject.h"
#include "TObjArray.h"
#include "TLinearFitter.h"

class AliExternalTrackParam;

class AliTPCcalibAlign:public TObject {
public:
  AliTPCcalibAlign();

  virtual ~AliTPCcalibAlign();

  void Process(const AliExternalTrackParam &t1,
	       const AliExternalTrackParam &t2,
	       Int_t s1,Int_t s2);
  void Eval();
  TLinearFitter* GetFitter12(Int_t s1,Int_t s2) {
    return static_cast<TLinearFitter*>(fFitterArray12[s1*72+s2]);
  }
  TLinearFitter* GetFitter9(Int_t s1,Int_t s2) {
    return static_cast<TLinearFitter*>(fFitterArray9[s1*72+s2]);
  }
  TLinearFitter* GetFitter6(Int_t s1,Int_t s2) {
    return static_cast<TLinearFitter*>(fFitterArray6[s1*72+s2]);
  }
  Bool_t GetTransformation12(Int_t s1,Int_t s2,TMatrixD &a);
  Bool_t GetTransformation9(Int_t s1,Int_t s2,TMatrixD &a);
  Bool_t GetTransformation6(Int_t s1,Int_t s2,TMatrixD &a);
private:
  void Process12(Double_t *t1,
		 Double_t *t2,
		 TLinearFitter *fitter);
  void Process9(Double_t *t1,
		Double_t *t2,
		TLinearFitter *fitter);
  void Process6(Double_t *t1,
		Double_t *t2,
		TLinearFitter *fitter);
  TLinearFitter* GetOrMakeFitter12(Int_t s1,Int_t s2) {
    //get or make fitter
    if (!fFitterArray12[s1*72+s2])
      fFitterArray12[s1*72+s2]=new TLinearFitter(12,"x0++x1++x2++x3++x4++x5++x6++x7++x8++x9++x10++x11");
    return GetFitter12(s1,s2);
  }
  TLinearFitter* GetOrMakeFitter9(Int_t s1,Int_t s2) {
    //get or make fitter
    if (!fFitterArray9[s1*72+s2])
      fFitterArray9[s1*72+s2]=new TLinearFitter(9,"x0++x1++x2++x3++x4++x5++x6++x7++x8");
    return GetFitter9(s1,s2);
  }
  TLinearFitter* GetOrMakeFitter6(Int_t s1,Int_t s2) {
    //get or make fitter
    if (!fFitterArray6[s1*72+s2])
      fFitterArray6[s1*72+s2]=new TLinearFitter(6,"x0++x1++x2++x3++x4++x5");
    return GetFitter6(s1,s2);
  }
  TObjArray fFitterArray12;
  TObjArray fFitterArray9;
  TObjArray fFitterArray6;
  Int_t fPoints[72*72];

  ClassDef(AliTPCcalibAlign,1)
};

#endif
