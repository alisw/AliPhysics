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

#include <TMath.h>
#include <TObject.h>
class TH1F;

class AliQuenchingWeights : public TObject {
 public:
  enum kECMethod {kDefault=0,kReweight=1,kReweightCont=2};

  AliQuenchingWeights();
  AliQuenchingWeights(const AliQuenchingWeights& a);
  AliQuenchingWeights& operator=(const AliQuenchingWeights& a)
      {a.Copy(*this); return(*this);}
  virtual ~AliQuenchingWeights();

  void Reset();
  Int_t SampleEnergyLoss();
  Int_t SampleEnergyLoss(Int_t ipart, Double_t r);

  Double_t GetELossRandom(Int_t ipart, Double_t length, Double_t e=1.e10) const;
  Double_t CalcQuenchedEnergy(Int_t ipart, Double_t length, Double_t e)  const;
  Double_t GetELossRandom(Int_t ipart, TH1F *hell, Double_t e=1.e10) const;
  Double_t CalcQuenchedEnergy(Int_t ipart, TH1F *hell, Double_t e)  const;
  Double_t GetELossRandomK(Int_t ipart, Double_t I0, Double_t I1, Double_t e=1.e10);
  Double_t CalcQuenchedEnergyK(Int_t ipart, Double_t I0, Double_t I1, Double_t e);
  Double_t GetELossRandomKFast(Int_t ipart, Double_t I0, Double_t I1, Double_t e=1.e10);
  Double_t GetELossRandomKFastR(Int_t ipart, Double_t r, Double_t wc, Double_t e=1.e10);
  Double_t CalcQuenchedEnergyKFast(Int_t ipart, Double_t I0, Double_t I1, Double_t e);

  Double_t GetDiscreteWeight(Int_t ipart, Double_t I0, Double_t I1);
  Double_t GetDiscreteWeightR(Int_t ipart, Double_t r);
  void GetZeroLossProb(Double_t &p,Double_t &prw,Double_t &prwcont,
		       Int_t ipart,Double_t I0,Double_t I1,Double_t e=1.e10);
  void GetZeroLossProbR(Double_t &p,Double_t &prw, Double_t &prwcont,
			Int_t ipart,Double_t r,Double_t wc,Double_t e=1.e10);

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
    {return 0.5*q*l*l*fgkConvFmToInvGeV;}

  Double_t CalcWCbar(Double_t mu, Double_t l) const 
    {return 0.5*mu*mu*l*fgkConvFmToInvGeV;}

  Double_t CalcWC(Double_t l) const 
    {if(fMultSoft) return CalcWC(fQTransport,l);
     else return CalcWCbar(fMu,l);}

  Double_t CalcWCk(Double_t I1) const 
    {if(fMultSoft) return CalcWCk(fK,I1);
     else return -1;} //not implemented!

  Double_t CalcWCk(Double_t k, Double_t I1) const 
    {if(fMultSoft) return k*I1/fgkConvFmToInvGeV;
     else return -1;} //not implemented!

  Double_t CalcR(Double_t wc, Double_t l) const; 

  Double_t CalcRk(Double_t I0, Double_t I1) const
    {return CalcRk(fK,I0,I1);} 

  Double_t CalcRk(Double_t k, Double_t I0, Double_t I1) const; 

  Double_t CalcQk(Double_t I0, Double_t I1) const
    {return CalcQk(fK,I0,I1);} 

  Double_t CalcQk(Double_t k, Double_t I0, Double_t I1) const
    {return I0*I0/2./I1/fgkConvFmToInvGeV/fgkConvFmToInvGeV*k;}

  Double_t CalcLk(Double_t i0, Double_t i1) const
    {return 2.*i1/i0;}

  Int_t CalcLengthMax(Double_t q) const
    {Double_t l3max=fgkRMax/.5/q/fgkConvFmToInvGeV/fgkConvFmToInvGeV;
     return (Int_t)TMath::Power(l3max,1./3.);} 

  const TH1F* GetHisto(Int_t ipart,Double_t length) const;

  void SetMu(Double_t m=1.) {fMu=m;}
  void SetQTransport(Double_t q=1.) {fQTransport=q;}
  void SetK(Double_t k=4.e5) {fK=k;} //about 1 GeV^2/fm
  void SetECMethod(kECMethod type=kDefault);
  void SetLengthMax(Int_t l=20) {fLengthMax=l;}

  Float_t GetMu()           const {return fMu;}
  Float_t GetQTransport()   const {return fQTransport;}
  Float_t GetK()            const {return fK;}
  Bool_t  GetECMethod()     const {return fECMethod;}
  Bool_t  GetTablesLoaded() const {return fTablesLoaded;}
  Bool_t  GetMultSoft()     const {return fMultSoft;}
  Int_t   GetLengthMax()    const {return fLengthMax;}

  TH1F* ComputeQWHisto (Int_t ipart,Double_t medval,Double_t length)  const; 
  TH1F* ComputeQWHistoX(Int_t ipart,Double_t medval,Double_t length)  const; 
  TH1F* ComputeQWHistoX(Int_t ipart,Double_t r)                       const; 
  TH1F* ComputeELossHisto(Int_t ipart,Double_t medval,Double_t l,Double_t e=1.e10) const; 
  TH1F* ComputeELossHisto(Int_t ipart,Double_t medval,TH1F *hEll,Double_t e=1.e10) const; 
  TH1F* ComputeELossHisto(Int_t ipart,Double_t r)                                  const; 

  Double_t GetMeanELoss(Int_t ipart,Double_t medval,Double_t l) const;
  Double_t GetMeanELoss(Int_t ipart,Double_t medval,TH1F *hEll) const; 
  Double_t GetMeanELoss(Int_t ipart,Double_t r)                 const; 
  
  void PlotDiscreteWeights(Double_t len=4,Double_t qm=5)         const; 
  void PlotContWeights(Int_t itype,Double_t len)                 const;
  void PlotContWeightsVsL(Int_t itype,Double_t medval)           const;
  void PlotAvgELoss(Double_t len,Double_t qm=5,Double_t e=1.e10) const;
  void PlotAvgELoss(TH1F *hEll,Double_t e=1.e10)                 const;
  void PlotAvgELossVsL(Double_t e=1.e10)                         const;
  void PlotAvgELossVsPt(Double_t medval,Double_t len)            const;
  void PlotAvgELossVsPt(Double_t medval,TH1F *hEll)              const;

 protected:
  Int_t GetIndex(Double_t len) const;

  static const Double_t fgkConvFmToInvGeV; //conversion factor
  static const Int_t    fgkBins;           //number of bins for hists
  static const Double_t fgkMaxBin;         //max. value of wc
  static const Double_t fgkRMax;           //max. tabled value of R

  static Int_t fgCounter;//static instance counter
  Int_t fInstanceNumber; //instance number of class

  Bool_t fMultSoft;     //approximation type
  kECMethod fECMethod;     //energy constraint method
  Double_t fQTransport; //transport coefficient [GeV^2/fm]]
  Double_t fMu;         //Debye screening mass
  Double_t fK;          //proportional constant [fm]
  Int_t fLengthMax;     //maximum length
  Int_t fLengthMaxOld;  //maximum length used for histos

  //discrete and cont part of quenching for
  //both parton type and different lengths
  TH1F ***fHistos; //!
  TH1F *fHisto; //!

  // data strucs for tables
  Double_t fxx[400];      //sampled energy quark
  Double_t fxxg[400];     //sampled energy gluon
  Double_t fdaq[34];      //discrete weight quark
  Double_t fdag[34];      //discrete weight gluon
  Double_t fcaq[34][261]; //continuous weights quarks
  Double_t fcag[34][261]; //continuous weights gluons  
  Double_t frrr[34];      //r value quark
  Double_t frrrg[34];     //r value gluon
  Bool_t fTablesLoaded;   //tables loaded

  ClassDef(AliQuenchingWeights,1)    // Base class for Quenching Weights
};

#endif
