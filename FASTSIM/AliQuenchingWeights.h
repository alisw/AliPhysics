#ifndef ALIQUENCHINGWEIGHTS_H
#define ALIQUENCHINGWEIGHTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//----------------------------------------------------------------------------
//     Implementation of the class to calculate the parton energy loss
//  Based on the "BDMPS" quenching weights by C.A.Salgado and U.A.Wiedemann
//
//  References:
//   C.A.Salgado and U.A.Wiedemann, Phys.Rev.D68 (2003) 014008 [hep-ph/0302184]
//   A.Dainese, Eur.Phys.J.C, in press, [nucl-ex/0312005]             
//
//            Origin:  C. Loizides   constantin.loizides@cern.ch
//                     A. Dainese    andrea.dainese@pd.infn.it            
//----------------------------------------------------------------------------

#include <TObject.h>
class TH1F;

class AliQuenchingWeights : public TObject {
 public:
  enum kECMethod {kDefault=0,kReweight=1};

  AliQuenchingWeights();
  AliQuenchingWeights(const AliQuenchingWeights& a);
  virtual ~AliQuenchingWeights();

  void Reset();
  Int_t SampleEnergyLoss();
  Double_t GetELossRandom(Int_t ipart, Double_t length, Double_t e=1.e6) const;
  Double_t CalcQuenchedEnergy(Int_t ipart, Double_t length, Double_t e)  const;
  Double_t GetELossRandom(Int_t ipart, TH1F *hell, Double_t e=1.e6) const;
  Double_t CalcQuenchedEnergy(Int_t ipart, TH1F *hell, Double_t e)  const;

  //multiple soft scattering approximation
  Int_t InitMult(const Char_t *contall="$(ALICE_ROOT)/FASTSIM/data/cont_mult.all",
                 const Char_t *discall="$(ALICE_ROOT)/FASTSIM/data/disc_mult.all"); 

  //single hard scattering approximation
  Int_t InitSingleHard(const Char_t *contall="$(ALICE_ROOT)/FASTSIM/data/cont_lin.all",
                       const Char_t *discall="$(ALICE_ROOT)/FASTSIM/data/disc_lin.all"); 

  Int_t CalcMult(Int_t ipart, Double_t rrrr,Double_t xxxx,
                 Double_t &continuous,Double_t &discrete) const;
  Int_t CalcMult(Int_t ipart, 
		 Double_t w, Double_t qtransp, Double_t length,
                 Double_t &continuous,Double_t &discrete) const;
  Int_t CalcSingleHard(Int_t ipart, Double_t rrrr,Double_t xxxx,
		       Double_t &continuous,Double_t &discrete) const;
  Int_t CalcSingleHard(Int_t ipart, 
  		       Double_t w, Double_t mu, Double_t length,
                       Double_t &continuous,Double_t &discrete) const;

  Double_t CalcWC(Double_t q, Double_t l) const 
    {return 0.5*q*l*l*gkConvFmToInvGeV;}

  Double_t CalcWCbar(Double_t mu, Double_t l) const 
    {return 0.5*mu*mu*l*gkConvFmToInvGeV;}

  Double_t CalcWC(Double_t l) const 
    {if(fMultSoft) return CalcWC(fQTransport,l);
     else return CalcWCbar(fMu,l);}

  Double_t CalcR(Double_t wc, Double_t l) const; 

  Int_t CalcLengthMax(Double_t q) const
    {Double_t l3max=gkRMax/.5/q/gkConvFmToInvGeV/gkConvFmToInvGeV;
     return (Int_t)TMath::Power(l3max,1./3.);} 

  const TH1F* GetHisto(Int_t ipart,Int_t l) const;

  void SetMu(Double_t m=1.) {fMu=m;}
  void SetQTransport(Double_t q=1.) {fQTransport=q;}
  void SetECMethod(kECMethod type=kDefault);
  void SetLengthMax(Int_t l=20) {fLengthMax=l;}

  Float_t GetMu()           const {return fMu;}
  Float_t GetQTransport()   const {return fQTransport;}
  Bool_t  GetECMethod()     const {return fECMethod;}
  Bool_t  GetTablesLoaded() const {return fTablesLoaded;}
  Bool_t  GetMultSoft()     const {return fMultSoft;}
  Int_t   GetLengthMax()    const {return fLengthMax;}

  TH1F* ComputeQWHisto (Int_t ipart,Double_t medval,Double_t length)  const; 
  TH1F* ComputeQWHistoX(Int_t ipart,Double_t medval,Double_t length)  const; 
  TH1F* ComputeELossHisto (Int_t ipart,Double_t medval,Double_t l,Double_t e=1.e6) const; 
  TH1F* ComputeELossHisto (Int_t ipart,Double_t medval,TH1F *hEll,Double_t e=1.e6) const; 
  
  void PlotDiscreteWeights(Int_t len=4) const;
  void PlotContWeights(Int_t itype,Int_t len) const;
  void PlotContWeights(Int_t itype,Double_t medval) const;
  void PlotAvgELoss(Int_t len ,Double_t e=1.e6) const;
  void PlotAvgELoss(TH1F *hEll,Double_t e=1.e6) const;
  void PlotAvgELossVsPt(Double_t medval,Int_t len) const;
  void PlotAvgELossVsPt(Double_t medval,TH1F *hEll) const;

 protected:
  static const Double_t gkConvFmToInvGeV; 
  static const Double_t gkRMax; 
  static Int_t gCounter;//static instance counter
  Int_t fInstanceNumber;//instance number of class

  Bool_t fMultSoft;     //approximation type
  Bool_t fECMethod;     //energy constraint method
  Double_t fQTransport; //transport coefficient
  Double_t fMu;         //Debye screening mass
  Int_t fLengthMax;     //maximum length
  Int_t fLengthMaxOld;  //maximum length used for histos

  //discrete and cont part of quenching for
  //both parton type and different lengths
  TH1F ***fHistos; //!

  // data strucs for tables
  Double_t fxx[400];
  Double_t fxxg[400];
  Double_t fdaq[30];
  Double_t fdag[30];
  Double_t fcaq[30][261];
  Double_t fcag[30][261];  
  Double_t frrr[30];
  Double_t frrrg[30];
  Bool_t fTablesLoaded;

  ClassDef(AliQuenchingWeights,1)    // Base class for Quenching Weights
};

#endif
