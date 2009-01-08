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

 
class AliTPCcalibUnlinearity:public AliTPCcalibBase {
public:
  AliTPCcalibUnlinearity(); 
  AliTPCcalibUnlinearity(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibUnlinearity();
  //
  virtual void Process(AliTPCseed *track);
  virtual void Analyze(){return;}
  virtual void Terminate();
  virtual Long64_t Merge(TCollection* list){return 0;}
  //
  void ProcessTree(TTree * tree, Long64_t nmaxPoints);
  void AddPoint(Int_t sector, Int_t row, Double_t dz, Double_t dy, Double_t p2, Double_t p3, Double_t dr, Int_t npoints=1);
  //
  void MakeHisto();
  void ProcessDiff(AliTPCseed *track, Int_t isec);
  void DumpTree();
  void MakeFitters();
  void EvalFitters();
  TLinearFitter * GetFitterOutR(Int_t sector) {return (TLinearFitter*)fFittersOutR.At(sector);}
  TLinearFitter * GetFitterOutZ(Int_t sector) {return (TLinearFitter*)fFittersOutZ.At(sector);}
  TVectorD * GetParamOutR(Int_t sector) {return (TVectorD*)fParamsOutR.At(sector);}
  TVectorD * GetParamOutZ(Int_t sector) {return (TVectorD*)fParamsOutZ.At(sector);}
  //
  Double_t      GetDr(Int_t sector, Double_t dout, Double_t dr);
  Double_t      GetDz(Int_t sector, Double_t dout, Double_t dr);
  Double_t      GetGDr(Int_t stype,Float_t gx, Float_t gy,Float_t gz);
  //
  static Double_t      SGetDr(Int_t sector, Double_t dout, Double_t dr);
  static Double_t      SGetDz(Int_t sector, Double_t dout, Double_t dr);
  static Double_t      SGetGDr(Int_t stype,Float_t gx, Float_t gy,Float_t gz);

  static AliTPCcalibUnlinearity* Instance();
  void SetInstance(AliTPCcalibUnlinearity*param){fgInstance = param;}
  //TMatrixD * GetNormCovariance(Int_t sector, Int_t type);
  //TMatrixD * GetNormCovariance(Int_t sector, Int_t type);
  static void MakeQPosNormAll(TTree * chainres, AliTPCClusterParam * param, Int_t maxPoints);
public:
  THnSparse * fDiffHistoLine;    // matrix with cluster residuals - linear fit
  THnSparse * fDiffHistoPar;     // matrix with cluster residuals - parabolic fit
  //
  //
  TObjArray   fFittersOutR;      // Unlinearity fitters for radial distortion  - outer field cage
  TObjArray   fFittersOutZ;      // Unlinearity fitters for z      distortion  - outer field cage
  TObjArray   fParamsOutR;       // Parameters  for radial distortion  - outer field cage
  TObjArray   fParamsOutZ;      // Parameters  for z      distortion  - outer field cage
  TObjArray   fErrorsOutR;       // Parameters  for radial distortion  - outer field cage
  TObjArray   fErrorsOutZ;      // Parameters  for z      distortion  - outer field cage
  //

private:
  AliTPCcalibUnlinearity(const AliTPCcalibUnlinearity&); 
  AliTPCcalibUnlinearity& operator=(const AliTPCcalibUnlinearity&); 
 static AliTPCcalibUnlinearity*   fgInstance; //! Instance of this class (singleton implementation)
 
  ClassDef(AliTPCcalibUnlinearity, 1); 
};

#endif


