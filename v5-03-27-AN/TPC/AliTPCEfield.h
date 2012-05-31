#ifndef ALITPCEFIELD_H
#define ALITPCEFIELD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCEfield.h 35613 2009-10-16 03:24:40Z marian $ */


#include "TNamed.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TNamed.h"
class TTreeSRedirector;

class AliTPCEfield:public TNamed {
public:
  AliTPCEfield();  
  AliTPCEfield(const char * name, Int_t maxFreq, Bool_t is2D, Bool_t useLinear=kTRUE);  
  virtual ~AliTPCEfield(); 
  void SetRange(Double_t x0, Double_t x1, Double_t y0, Double_t y1, Double_t z00,Double_t z1=0);
  void AddBoundaryLine(Double_t x0,Double_t y0, Double_t z0, Double_t v0, Double_t x1, Double_t y1, Double_t z1, Double_t v1, Int_t id=0, Int_t npoints=100);
  TTree * GetTree(const char * tname="Boundary");
  //
  void MakeFitFunctions(Int_t maxFreq);
  void FitField();
  void DumpField(Double_t gridSize=5, Double_t step=0.5);
  //
  Double_t EvalField(Int_t ifun, Double_t x, Double_t y, Double_t z,  Int_t type=0); 
  Double_t Eval(Double_t x, Double_t y, Double_t z,  Int_t type=0); 
  //
  static Double_t EvalS(Double_t x, Double_t y, Double_t z,  Int_t type=0); 
  //
  
  Double_t Field(Int_t ftype, Double_t ifx, Double_t ify, Double_t ifz, Double_t x, Double_t y, Double_t z);
  Double_t FieldDn(Int_t ftype, Double_t ifx, Double_t ify, Double_t ifz, Int_t dn, Double_t x, Double_t y, Double_t z);
  TMatrixD* MakeCorrelation(TMatrixD &matrix);

  // get rid of numerical instabilities
  Double_t SinHNorm(Double_t x, Double_t norm){ return 0.5*(TMath::Exp(x-norm)-TMath::Exp(-x-norm));}
  Double_t CosHNorm(Double_t x, Double_t norm){ return 0.5*(TMath::Exp(x-norm)+TMath::Exp(-x-norm));}
 public:
  Double_t fMin[3];      // range of coordinates from Min to Max
  Double_t fMax[3];      //  
  Double_t fScale;       // scaling factor
  Int_t    fMaxFreq;     // maximal frequency of expansion
  Bool_t   fIs2D;        // flag for 2D field
  Bool_t   fUseLinear;   // flag to use also linear term of the field 
  //
  TTreeSRedirector * fWorkspace;    //! workspace
  TMatrixD  *fFitFunctions;         // fit function description
  TVectorD  *fFitParam;             // fit parameters - coeficients
  TMatrixD  *fFitCovar;             // fit covariance
  TLinearFitter *fFitter;           // linear fitter - temporary solution - integrals to be calculated
  static AliTPCEfield* fgInstance;  // instance of fied  - for visualization
  ClassDef(AliTPCEfield,1)
};

#endif
