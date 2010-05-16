#ifndef ALITPCKALMANALIGN_H
#define ALITPCKALMANALIGN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TNamed.h"
#include "TMatrixD.h"
class TTreeSRedirector;
class TObjArray;
class AliTPCcalibAlign;
class TTreeSRedirector;
class TTree;
class AliTPCCalPad;


class AliTPCkalmanAlign: public TNamed{
public:
  AliTPCkalmanAlign();
  AliTPCkalmanAlign(const char* name, const char* title);
  void ReadAlign(const char *fname="CalibObjects.root");
  void MakeGlobalAlign();
  void DrawDeltaAlign();
  void UpdateOCDBTime0( AliTPCCalPad  *pad, Int_t ustartRun, Int_t uendRun,  const char* storagePath );
  static void UpdateAlign1D(Double_t delta, Double_t sigma, Int_t s1, Int_t s2, TMatrixD &param, TMatrixD &covar);
  static void BookAlign1D(TMatrixD &param, TMatrixD &covar, Double_t sigma, Double_t mean);
  void DumpOldAlignment(TTreeSRedirector *pcstream);
  void MakeNewAlignment(Bool_t add,TTreeSRedirector *pcstream=0);
  void DrawAlignmentTrends();
  void FitCE();
  static void MakeAliasCE(TTree * chain);
public:
  AliTPCcalibAlign * fCalibAlign;    // kalman alignemnt
  TClonesArray     *fOriginalAlign;  // original alignment objects
  TClonesArray     *fNewAlign;       // new alignment objects
  //
  AliTPCCalPad     *fPadTime0;       // pad time0 - for z alignment
  //                                 // time offset parameterization
  TObjArray       *fFitCEGlobal;    // vector of parameter of the CE fits
  TObjArray       *fFitCELocal;     // vector of parameter delta to global
  //
  TMatrixD * fDelta1D[4];            // deltas
  TMatrixD * fCovar1D[4];            // covariance
private:
  AliTPCkalmanAlign(const AliTPCkalmanAlign&);
  AliTPCkalmanAlign &operator=(const AliTPCkalmanAlign&);
  ClassDef(AliTPCkalmanAlign,1);
};

#endif

