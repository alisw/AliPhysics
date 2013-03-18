#ifndef ALIDIELECTRONBTOJPSITOELECDFFITFCNFITTER_H
#define ALIDIELECTRONBTOJPSITOELECDFFITFCNFITTER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//                      Class AliDielectronBtoJPSItoEleCDFfitFCNfitter
//                          Origin: A. Mastroserio
//     Interface class to fit invariant mass and pseudoproperdecay length separately
//     Contact: Annalisa.Mastroserio@ba.infn.it
//_________________________________________________________________________


#include "AliDielectronBtoJPSItoEleCDFfitFCN.h"
class TF1;
class TF2;
class TH1F;
class TH2F;
class TArrayD;

class AliDielectronBtoJPSItoEleCDFfitFCNfitter {

 public:		
  AliDielectronBtoJPSItoEleCDFfitFCNfitter();		
  virtual ~AliDielectronBtoJPSItoEleCDFfitFCNfitter();


  // INVARIANT MASS
  enum {kInvMassSignal=5, kInvMassBkg=4, kInvMassTotal=9};

  void SetInvMassParameters(const Double_t massPar[kInvMassTotal]);
  void GetInvMassParameters(TArrayD &massPar);
  void SetInvMassSignalParameters(const Double_t massPar[kInvMassSignal]);
  void GetInvMassSignalParameters(TArrayD &massPar);
  void SetInvMassBkgParameters(const Double_t massPar[kInvMassBkg]);
  void GetInvMassBkgParameters(TArrayD &massPar);		

  Double_t CDFInvMassSignal(const Double_t *x, const Double_t *par);
  Double_t CDFInvMassBkg(const Double_t *x, const Double_t *par);
  Double_t CDFInvMassTotal(const Double_t *x, const Double_t *par);
  void SetInvMass(TH1F *mass) {fInvMass=mass;}	
  void SetParameterToFixInInvMass(Bool_t fixed[kInvMassTotal]) {for(Int_t i=0; i<kInvMassTotal; i++) 
   fParameterInvMassToFix[i]=fixed[i];}
  TH1F * FitInvMass(Double_t norm[2], Double_t mMin=0, Double_t mMax=0);
  TH1F * FitInvMassSignal( Double_t norm=1, Double_t mMin=0, Double_t mMax=0);


  // RESOLUTION FUNC
  // use SetResolutionConstants to set the parameters for FF, FS, SS		
  enum {kPseudo=9, kPseudoBkg=7};
  void SetPseudoProper(TH1F *hX) {if(fX)delete fX; fX=hX;}
  void GetPseudoProperParameters(TArrayD &xPar, Int_t type);
  Double_t CDFResolutionFunction(const Double_t *x, const Double_t *par);
  Double_t PsProperBackFitFunc(const Double_t* x, const Double_t* par);


  void SetParameterToFixInX(Bool_t fixed[kPseudo]) {for(Int_t i=0; i<kPseudo; i++) 
   fParameterXToFix[i]=fixed[i];}	
  void SetPseudoProperBkg(TH2F *hX2D) {fX2D=hX2D;}
  void SetPseudoProperBkgParameters(Double_t xBkgPar[kPseudoBkg]);
  void GetPseudoProperBkgParameters(TArrayD& bkgParam);
  void SetParameterToFixInXbkg(Bool_t fixed[kPseudoBkg]) {for(Int_t i=0; i<kPseudoBkg; i++) fParameterXbkgToFix[i]=fixed[i];}
  TH1F *FitResolutionFunction(Double_t xmin, Double_t xmax, Int_t type, Double_t norm=1);	
  TH2F *FitBkgPsudoProper(Double_t xMin, Double_t xMax, Double_t norm);


  // GENERAL
  void SetFitOption(char *opt) {fFitOpt=opt;}
  void SetCrystalBallFunction(Bool_t isCB) {fFCN->SetCrystalBallFunction(isCB);} 
  void PrintParamStatus() {fFCN->PrintStatus();}
  AliDielectronBtoJPSItoEleCDFfitFCN * GetFCN() {return fFCN;}

 protected :

  TH1F *fInvMass; // 1D inv mass distribution
  TH1F *fX;       // 1D X distribution
  TH2F *fX2D;     // 2D X distribution
  AliDielectronBtoJPSItoEleCDFfitFCN *fFCN;  // fitter

  TString fFitOpt; // fit options
  Bool_t fParameterInvMassToFix[kInvMassTotal]; // booleans to fix single parameters
  Bool_t fParameterXToFix[kPseudo];              // booleans to fix single parameters
  Bool_t fParameterXbkgToFix[kPseudoBkg];        // booleans to fix single parameters

  AliDielectronBtoJPSItoEleCDFfitFCNfitter(const AliDielectronBtoJPSItoEleCDFfitFCNfitter& source); 
  AliDielectronBtoJPSItoEleCDFfitFCNfitter& operator=(const AliDielectronBtoJPSItoEleCDFfitFCNfitter& source);  

  ClassDef (AliDielectronBtoJPSItoEleCDFfitFCNfitter,1);     
};

#endif
