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
#include "AliTPCcalibBase.h"
#include "TH1.h"

class AliExternalTrackParam;
class AliTPCseed;
class TGraphErrors;

class AliTPCcalibAlign:public AliTPCcalibBase {
public:
  enum HistoType {kY=0, kZ =1, kPhi=2, kTheta=3};
  enum FitType{ k6=0, k9=1, k12};
  AliTPCcalibAlign();
  AliTPCcalibAlign(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibAlign();
  virtual void Process(AliTPCseed *track);
  virtual void Analyze();
  virtual void Terminate();  
  virtual Long64_t Merge(TCollection* list);
  //
  virtual void EvalFitters();
  TH1 * GetHisto(HistoType type, Int_t s1, Int_t s2, Bool_t force=kFALSE);
  void  MakeTree(const char *fname="alignTree.root");
  TGraphErrors * MakeGraph(Int_t sec0, Int_t sec1, Int_t dsec, 
			   Int_t i0, Int_t i1, FitType type); 
  void ProcessTracklets(const AliExternalTrackParam &t1,
			const AliExternalTrackParam &t2,
			const AliTPCseed * seed,
			Int_t s1,Int_t s2);
  inline Int_t GetIndex(Int_t s1,Int_t s2){return 72*s1+s2;}
  //
  inline TLinearFitter* GetFitter12(Int_t s1,Int_t s2);
  inline TLinearFitter* GetFitter9(Int_t s1,Int_t s2);
  inline TLinearFitter* GetFitter6(Int_t s1,Int_t s2);
  //
  Bool_t GetTransformation12(Int_t s1,Int_t s2,TMatrixD &a);
  Bool_t GetTransformation9(Int_t s1,Int_t s2,TMatrixD &a);
  Bool_t GetTransformation6(Int_t s1,Int_t s2,TMatrixD &a);
  
  void ProcessDiff(const AliExternalTrackParam &t1,
		   const AliExternalTrackParam &t2,
		   const AliTPCseed *seed,
		   Int_t s1,Int_t s2);

//   Bool_t GetTransformationCovar12(Int_t s1,Int_t s2,TMatrixD &a, Bool_t norm=kFALSE);
//   Bool_t GetTransformationCovar9(Int_t s1,Int_t s2,TMatrixD &a, Bool_t norm=kFALSE);
//   Bool_t GetTransformationCovar6(Int_t s1,Int_t s2,TMatrixD &a, Bool_t norm=kFALSE);
  void Add(AliTPCcalibAlign * align);
  Int_t *GetPoints() {return fPoints;}
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliESDEvent *event){AliTPCcalibBase::Process(event);}
private:
  void FillHisto(const AliExternalTrackParam &t1,
			const AliExternalTrackParam &t2,
			Int_t s1,Int_t s2);

  void Process12(const Double_t *t1, const Double_t *t2,
		 TLinearFitter *fitter);
  void Process9(Double_t *t1, Double_t *t2, TLinearFitter *fitter);
  void Process6(Double_t *t1, Double_t *t2, TLinearFitter *fitter);
  TLinearFitter* GetOrMakeFitter12(Int_t s1,Int_t s2);
  TLinearFitter* GetOrMakeFitter9(Int_t s1,Int_t s2);
  TLinearFitter* GetOrMakeFitter6(Int_t s1,Int_t s2);
  TObjArray fDphiHistArray;    // array of residual histograms  phi
  TObjArray fDthetaHistArray;  // array of residual histograms  theta
  TObjArray fDyHistArray;      // array of residual histograms  y
  TObjArray fDzHistArray;      // array of residual histograms  z
  TObjArray fFitterArray12;    // array of fitters
  TObjArray fFitterArray9;     // array of fitters
  TObjArray fFitterArray6;     // array of fitters
  Int_t fPoints[72*72];        // number of points in the fitter
  ClassDef(AliTPCcalibAlign,1)
};


TLinearFitter* AliTPCcalibAlign::GetFitter12(Int_t s1,Int_t s2) {
  return static_cast<TLinearFitter*>(fFitterArray12[GetIndex(s1,s2)]);
}
TLinearFitter* AliTPCcalibAlign::GetFitter9(Int_t s1,Int_t s2) {
  return static_cast<TLinearFitter*>(fFitterArray9[GetIndex(s1,s2)]);
}
TLinearFitter* AliTPCcalibAlign::GetFitter6(Int_t s1,Int_t s2) {
  return static_cast<TLinearFitter*>(fFitterArray6[GetIndex(s1,s2)]);
}





#endif
