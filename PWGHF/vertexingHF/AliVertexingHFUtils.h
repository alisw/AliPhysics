#ifndef ALIVERTEXINGHFUTILS_H
#define ALIVERTEXINGHFUTILS_H


/* $Id:  $ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class with functions useful for different D2H analyses        //
// - event plane resolution                                      //
// - <pt> calculation with side band subtraction                 //
// - tracklet multiplcity calculation                            //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TH1F.h"

class AliVertexingHFUtils : public TObject{
 public:
  AliVertexingHFUtils();
  AliVertexingHFUtils(Int_t k);
  virtual ~AliVertexingHFUtils() {};

  // Functions for Event plane resolution
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
  static Double_t GetSubEvResolLowLim(const TH1F* hSubEvCorr){
    if(hSubEvCorr) return TMath::Sqrt(hSubEvCorr->GetMean()-hSubEvCorr->GetMeanError());
    else return 1.;
  }
  static Double_t GetSubEvResolHighLim(const TH1F* hSubEvCorr){
    if(hSubEvCorr) return TMath::Sqrt(hSubEvCorr->GetMean()+hSubEvCorr->GetMeanError());
    else return 1.;
  }
  static Double_t GetFullEvResol(Double_t resSub, Int_t k=1);
  static Double_t GetFullEvResol(const TH1F* hSubEvCorr, Int_t k=1);
  static Double_t GetFullEvResolLowLim(const TH1F* hSubEvCorr, Int_t k=1);
  static Double_t GetFullEvResolHighLim(const TH1F* hSubEvCorr, Int_t k=1);

  // Functions for tracklet multiplcity calculation
  void SetEtaRangeForTracklets(Double_t mineta, Double_t maxeta){
    fMinEtaForTracklets=mineta; 
    fMaxEtaForTracklets=maxeta;
  }
  static Int_t GetNumberOfTrackletsInEtaRange(AliAODEvent* ev, Double_t mineta, Double_t maxeta);
  Int_t GetNumberOfTrackletsInEtaRange(AliAODEvent* ev) const {
    return GetNumberOfTrackletsInEtaRange(ev,fMinEtaForTracklets,fMaxEtaForTracklets);
  }

 private:

  Int_t fK;             // ratio of measured harmonic to event plane harmonic
  Double_t fSubRes;     // sub-event resolution = sqrt(<cos[n(phiA-phiB)] >)
  Double_t fMinEtaForTracklets; // min eta for counting tracklets
  Double_t fMaxEtaForTracklets; // min eta for counting tracklets

  ClassDef(AliVertexingHFUtils,0) 
};
#endif
