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
  void ProcessTree(TTree * tree, Int_t nmaxPoints);
  void AddPoint(Int_t sector, Int_t row, Float_t dz, Float_t dy, Float_t p2, Float_t p3, Float_t dr, Int_t npoints=1);
  //
  void MakeHisto();
  void ProcessDiff(AliTPCseed *track, Int_t isec);
  void DumpTree();
  void MakeFitters();
  void EvalFitters();
  //
public:
  THnSparse * fDiffHistoLine;    // matrix with cluster residuals - linear fit
  THnSparse * fDiffHistoPar;     // matrix with cluster residuals - parabolic fit
  //
  //
  TObjArray   fFittersOutR;      // Unlinearity fitters for radial distortion  - outer field cage
  TObjArray   fFittersOutZ;      // Unlinearity fitters for z      distortion  - outer field cage
  //
  
private:
  AliTPCcalibUnlinearity(const AliTPCcalibUnlinearity&); 
  AliTPCcalibUnlinearity& operator=(const AliTPCcalibUnlinearity&); 

  ClassDef(AliTPCcalibUnlinearity, 1); 
};

#endif


