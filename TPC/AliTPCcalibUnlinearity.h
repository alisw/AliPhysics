#ifndef ALITPCUNLINEARITY_H
#define ALITPCUNLINEARITY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
#include "TH2F.h"
#include "TF1.h"
#include "TArrayD.h"
#include "TObjArray.h"
#include "TTreeStream.h"
#include "TVectorD.h"

class TH1F;
class TH3F;
class TH2F;
class THnSparse;
class TList;
class AliESDEvent;
class AliESDtrack;
class TLinearFitter;
class AliTPCClusterParam;

 
class AliTPCcalibUnlinearity:public AliTPCcalibBase {
public:
  AliTPCcalibUnlinearity(); 
  AliTPCcalibUnlinearity(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibUnlinearity();
  //
  virtual void Process(AliTPCseed *track);
  virtual void Analyze(){return;}
  virtual void Terminate();
  virtual Long64_t Merge(TCollection* list);
  void    Add(AliTPCcalibUnlinearity * calib);
  //
  void ProcessTree(TTree * tree, Long64_t nmaxPoints);
  void AddPoint(Int_t sector, Double_t cx, Double_t cy, Double_t cz, Double_t ty, Double_t tz,  Double_t ky, Double_t kz, Int_t npoints=1);
  void AddPointRPHI(Int_t sector, Double_t cx, Double_t cy, Double_t cz, Double_t ty, Double_t tz,  Double_t ky, Double_t kz, Int_t npoints=1);
  //
  void MakeHisto();
  void ProcessDiff(AliTPCseed *track, Int_t isec);
  void AlignOROC(AliTPCseed *track, Int_t isec);
  //
  void DumpTree(const char *fname="unlinResidual.root");
  void Init();
  void EvalFitters();
  TLinearFitter * GetFitterOutR(Int_t sector) {return (TLinearFitter*)fFittersOutR.At(sector);}
  TLinearFitter * GetFitterOutZ(Int_t sector) {return (TLinearFitter*)fFittersOutZ.At(sector);}
  TVectorD * GetParamOutR(Int_t sector) {return (TVectorD*)fParamsOutR.At(sector);}
  TVectorD * GetParamOutZ(Int_t sector) {return (TVectorD*)fParamsOutZ.At(sector);}
  //
  //TMatrixD * GetNormCovariance(Int_t sector, Int_t type);
  //TMatrixD * GetNormCovariance(Int_t sector, Int_t type);
  static void MakeQPosNormAll(TTree * chainres, AliTPCClusterParam * param, Int_t maxPoints);
  void     Process(AliESDEvent *event) {AliTPCcalibBase::Process(event);};
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
public:
  Bool_t      fInit;              // initialization flag
  THnSparse * fDiffHistoLineY;    // matrix with cluster residuals - linear fit
  THnSparse * fDiffHistoParY;     // matrix with cluster residuals - parabolic fit
  THnSparse * fDiffHistoLineZ;    // matrix with cluster residuals - linear fit
  THnSparse * fDiffHistoParZ;     // matrix with cluster residuals - parabolic fit
  //
  // Outer residula fit
  //
  TObjArray   fFittersOutR;      // Unlinearity fitters for radial distortion  - outer field cage
  TObjArray   fFittersOutZ;      // Unlinearity fitters for z      distortion  - outer field cage
  TObjArray   fParamsOutR;       // Parameters  for radial distortion  - outer field cage
  TObjArray   fParamsOutZ;      // Parameters  for z      distortion  - outer field cage
  TObjArray   fErrorsOutR;       // Parameters  for radial distortion  - outer field cage
  TObjArray   fErrorsOutZ;      // Parameters  for z      distortion  - outer field cage
  //
  // R-phi residual histogram
  //
  TObjArray   fDistRPHIPlus;     // edge effect histograms  - plus  direction 
  TObjArray   fDistRPHIMinus;    // edge effect histograms  - minus direction
  //
  // Quadrant fitters
  //
   TObjArray   fFitterQuadrantY;        //qudrant misalignemnt fit Y
   TObjArray   fFitterQuadrantPhi;      //qudrant misalignemnt fit Phi
private:
  AliTPCcalibUnlinearity(const AliTPCcalibUnlinearity&); 
  AliTPCcalibUnlinearity& operator=(const AliTPCcalibUnlinearity&); 
 
  ClassDef(AliTPCcalibUnlinearity, 4); 
};

#endif


