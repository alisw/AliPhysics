#ifndef ALIEVENTPLANERESOLUTION_H
#define ALIEVENTPLANERESOLUTION_H


/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to compute event plane resolution for flow analyses     //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TH1F.h"

class AliEventPlaneResolution : public TObject{
 public:
  AliEventPlaneResolution();
  AliEventPlaneResolution(Int_t k);
  virtual ~AliEventPlaneResolution() {};

  void SetK(Int_t k){fK=k;}
  void SetSubEvResol(Double_t res){fSubRes=res;}
  void SetSubEventHisto(const TH1F* hSub){
    fSubRes=GetSubEvResol(hSub);
  }

  Int_t GetK() const  {return fK;}
  Double_t GetSubEvResol() const  {return fSubRes;}

  Double_t Pol(Double_t x) const {return Pol(x,fK);}
  Double_t FindChi() const {return FindChi(fSubRes,fK);}
  Double_t GetFullEvResol() const {return GetFullEvResol(fSubRes,fK);}

  static Double_t FindChi(Double_t res,  Int_t k=1);
  static Double_t Pol(Double_t x, Int_t k);
  static Double_t ResolK1(Double_t x);
  static Double_t GetSubEvResol(const TH1F* hSubEvCorr){
    if(hSubEvCorr) return TMath::Sqrt(hSubEvCorr->GetMean());
    else return 1.;
  }
  static Double_t GetFullEvResol(Double_t resSub, Int_t k=1);
  static Double_t GetFullEvResol(const TH1F* hSubEvCorr, Int_t k=1);

 private:

  Int_t fK;             // ratio of measured harmonic to event plane harmonic
  Double_t fSubRes;     // sub-event resolution = sqrt(<cos[n(phiA-phiB)] >)

  ClassDef(AliEventPlaneResolution,0) 
};
#endif
