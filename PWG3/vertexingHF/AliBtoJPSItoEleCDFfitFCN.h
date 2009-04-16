#ifndef ALIBTOJPSITOELECDFFITFCN_H
#define ALIBTOJPSITOELECDFFITFCN_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//                      Class AliBtoJPSItoEleCDFfitFCN
//                    Definition of main function used in 
//                     unbinned log-likelihood fit for
//                 the channel B -> JPsi + X -> e+e- + X
//      
//                          Origin: C.Di Giglio
//     Contact: Carmelo.Digiglio@ba.infn.it , Giuseppe.Bruno@ba.infn.it
//_________________________________________________________________________

#include <TNamed.h>
#include <TDatabasePDG.h>
#include "TH1F.h"
#include "TMath.h"

enum IntegralType {kBkg, 
                   kBkgNorm, 
                   kSig, 
                   kSigNorm};

enum PtBins       {kallpt, 
                   kptbin1,kptbin2,kptbin3,
                   kptbin4,kptbin5,kptbin6,
                   kptbin7,kptbin8,kptbin9};
//_________________________________________________________________________________________________
class AliBtoJPSItoEleCDFfitFCN : public TNamed {
 public:
  //
  AliBtoJPSItoEleCDFfitFCN();
  AliBtoJPSItoEleCDFfitFCN(const AliBtoJPSItoEleCDFfitFCN& source); 
  AliBtoJPSItoEleCDFfitFCN& operator=(const AliBtoJPSItoEleCDFfitFCN& source);  
  virtual ~AliBtoJPSItoEleCDFfitFCN();

  Double_t EvaluateLikelihood(const Double_t* pseudoproperdecaytime,
                              const Double_t* invariantmass, const Int_t ncand);
 
  Double_t GetFPlus() const {return fFPlus;}
  Double_t GetFMinus() const {return fFMinus;}
  Double_t GetFSym() const {return fFSym;}
  Double_t GetRadius() const {return fParameters[0];}
  Double_t GetTheta() const {return fParameters[1];}
  Double_t GetPhi() const {return fParameters[2];}
  Double_t GetLamPlus() const {return fParameters[3];}
  Double_t GetLamMinus() const {return fParameters[4];}
  Double_t GetLamSym() const {return fParameters[5];}
  Double_t GetMassSlope() const {return fParameters[6];}
  Double_t GetFractionJpsiFromBeauty() const {return fParameters[7];}
  Double_t GetFsig() const {return fParameters[8];}
  Double_t GetCrystalBallMmean() const {return fParameters[9];}
  Double_t GetCrystalBallNexp() const {return fParameters[10];}
  Double_t GetCrystalBallSigma() const {return fParameters[11];}
  Double_t GetCrystalBallAlpha() const {return fParameters[12];}
  Double_t GetIntegral() const {return fIntegral;}
  Bool_t GetCrystalBallParam() const {return fCrystalBallParam;}

  void SetFPlus(Double_t plus) {fFPlus = plus;}
  void SetFMinus(Double_t minus) {fFMinus = minus;}
  void SetFSym(Double_t sym) {fFSym = sym;}
  void SetRadius(Double_t radius) {fParameters[0] = radius;}
  void SetTheta(Double_t theta) {fParameters[1] = theta;}
  void SetPhi(Double_t phi) {fParameters[2] = phi;}
  void SetLamPlus(Double_t lamplus) {fParameters[3] = lamplus;}
  void SetLamMinus(Double_t lamminus) {fParameters[4] = lamminus;}
  void SetLamSym(Double_t lamsym) {fParameters[5] = lamsym;}
  void SetMassSlope(Double_t slope) {fParameters[6] = slope;}
  void SetFractionJpsiFromBeauty(Double_t B) {fParameters[7] = B;}
  void SetFsig(Double_t Fsig) {fParameters[8] = Fsig;}
  void SetCrystalBallMmean(Double_t CrystalBallMmean) {fParameters[9] = CrystalBallMmean;}
  void SetCrystalBallNexp(Double_t CrystalBallNexp) {fParameters[10] = CrystalBallNexp;}
  void SetCrystalBallSigma(Double_t CrystalBallSigma) {fParameters[11] = CrystalBallSigma;}
  void SetCrystalBallAlpha(Double_t CrystalBallAlpha) {fParameters[12] = CrystalBallAlpha;}

  void SetAllParameters(const Double_t* parameters);
  void SetIntegral(Double_t integral) {fIntegral = integral;}
  void SetCsiMC(const TH1F* MCtemplate) {fhCsiMC = (TH1F*)MCtemplate->Clone("fhCsiMC");}
  void SetResolutionConstants(Int_t BinNum);
  void SetMassWndHigh(Double_t limit) { fMassWndHigh = TDatabasePDG::Instance()->GetParticle(443)->Mass() + limit ;}//here use pdg code instead
  void SetMassWndLow(Double_t limit) { fMassWndLow = TDatabasePDG::Instance()->GetParticle(443)->Mass() - limit ;}//here use pdg code instead
  void SetCrystalBallParam(Bool_t okCB) {fCrystalBallParam = okCB;}

  void ConvertFromSpherical() { fFPlus  = TMath::Power((fParameters[0]*TMath::Cos(fParameters[1])),2.);
                                fFMinus = TMath::Power((fParameters[0]*TMath::Sin(fParameters[1])*TMath::Sin(fParameters[2])),2.);
                                fFSym   = TMath::Power((fParameters[0]*TMath::Sin(fParameters[1])*TMath::Cos(fParameters[2])),2.);} 

  void ComputeIntegral(); 

  void ReadMCtemplates(Int_t BinNum);

  void PrintStatus();

 private:  
  //
  Double_t fParameters[13];        /*  par[0]  = fRadius;                
                                       par[1]  = fTheta;
                                       par[2]  = fPhi;
                                       par[3]  = fOneOvLamPlus;
                                       par[4]  = fOneOvLamMinus;
                                       par[5]  = fOneOvLamSym;
                                       par[6]  = fMassBkgSlope;
                                       par[7]  = fFractionJpsiFromBeauty;
                                       par[8]  = fFsig;
                                       par[9]  = fCrystalBallMmean;
                                       par[10] = fCrystalBallNexp;
                                       par[11] = fCrystalBallSigma;
                                       par[12] = fCrystalBallAlpha;*/

  Double_t fFPlus;                  // parameters of the log-likelihood function
  Double_t fFMinus;                 // Slopes of the x distributions of the background 
  Double_t fFSym;                   // functions 

  Double_t fIntegral;               // integral values of log-likelihood function terms
  TH1F *fhCsiMC;                    // X distribution used as MC template for JPSI from B
  Double_t fResolutionConstants[7]; // constants for the parametrized resolution function R(X)
  Double_t fMassWndHigh;            // JPSI Mass window higher limit
  Double_t fMassWndLow;             // JPSI Mass window lower limit
  Bool_t fCrystalBallParam;         // Boolean to switch to Crystall Ball parameterisation

  ////

  Double_t EvaluateCDFfunc(Double_t x, Double_t m) const ;

  Double_t EvaluateCDFfuncNorm(Double_t x, Double_t m) const ;

  ////

  Double_t EvaluateCDFfuncSignalPart(Double_t x, Double_t m) const ;      // Signal part 

  Double_t EvaluateCDFDecayTimeSigDistr(Double_t x) const ;

  Double_t EvaluateCDFInvMassSigDistr(Double_t m) const ;

  Double_t EvaluateCDFfuncBkgPart(Double_t x,Double_t m) const ;          // Background part

  Double_t EvaluateCDFDecayTimeBkgDistr(Double_t x) const ;
 
  Double_t EvaluateCDFInvMassBkgDistr(Double_t m) const ;

  ////

  Double_t FunB(Double_t x) const ;

  Double_t FunP(Double_t x) const ;

  Double_t CsiMC(Double_t x) const ;

  Double_t FunBkgPos(Double_t x) const ;

  Double_t FunBkgNeg(Double_t x) const ;

  Double_t FunBkgSym(Double_t x) const ;

  Double_t ResolutionFunc(Double_t x) const ;
  //
  ClassDef (AliBtoJPSItoEleCDFfitFCN,0);         // Unbinned log-likelihood fit 

};

#endif
