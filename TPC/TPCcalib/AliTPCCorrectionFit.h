#ifndef ALITPCCORRECTIONFIT_H
#define ALITPCCORRECTIONFIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//
#include "TNamed.h"
class TObjArray;
class THnBase;
class THnSparse;
class TTreeSRedirector;
class TTree;

class AliTPCCorrectionFit:public TNamed { 
public:
  AliTPCCorrectionFit();
  virtual ~AliTPCCorrectionFit();
  //
  //
  //
  //static Double_t EvalAt(Double_t phi, Double_t refX, Double_t theta,Double_t evalX,  Int_t corr, Int_t ptype, Float_t wt=100, Float_t t1=1, Float_t t2=1);
  static Double_t EvalAtPar(Double_t phi, Double_t snp, Double_t refX, Double_t evalX, Double_t theta, Int_t corr, Int_t ptype, Int_t nstep, Float_t wt=100,  Float_t t1=1, Float_t t2=1, Bool_t *pstatus=0);
  static Double_t EvalAtHelix(Double_t phi, Double_t snp, Double_t refX, Double_t evalX, Double_t theta, Int_t corr, Int_t ptype, Int_t nstep,  Float_t wt=100, Float_t t1=1, Float_t t2=1, Bool_t *pstatus=0);
  //
  // Make distortion maps
  //
  static void  CreateAlignMaps(Double_t  bz, Int_t run);
  static void  MakeClusterDistortionMap(THnBase * hisInput, TTreeSRedirector *pcstream, Int_t ptype, Int_t dtype=0);
  static void  MakeDistortionMap(THnSparse * his0, TTreeSRedirector *pcstream, const char* hname, Int_t run,  Float_t refX, Int_t type, Int_t integ, Double_t bz);
  static void  MakeDistortionMapCosmic(THnSparse * his0, TTreeSRedirector *pcstream, const char* hname, Int_t run,  Float_t refX, Int_t type);
  static void  MakeDistortionMapSector(THnSparse * his0, TTreeSRedirector *pcstream, const char* hname, Int_t run, Int_t type, Double_t bz);
  //
  // Create a distortion trees with the numerical derivatives
  // 
  //static void MakeLaserDistortionTree(TTree* tree, TObjArray *corrArray, Int_t itype);
  static void MakeTrackDistortionTree(TTree *tinput, Int_t dtype, Int_t ptype, const TObjArray * corrArray, Int_t step=1, Int_t offset=0, Bool_t debug=0);
  static void MakeSectorDistortionTree(TTree *tinput, Int_t dtype, Int_t ptype, const TObjArray * corrArray, Int_t step=1, Int_t offset=0, Bool_t debug=0);
  void SetMaxChi2(Double_t maxChi2) {fgMaxChi2HelixAt=maxChi2;}
private:

private:
  static Double_t fgMaxChi2HelixAt;  // max chi2 per point in helix fit 
  AliTPCCorrectionFit& operator=(const AliTPCCorrectionFit&); // not implemented
  AliTPCCorrectionFit(const AliTPCCorrectionFit&); // not implemented
  ClassDef(AliTPCCorrectionFit,1)
};

#endif
