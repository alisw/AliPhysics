#ifndef ALIVERTEXINGHFUTILS_H
#define ALIVERTEXINGHFUTILS_H


/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class with functions useful for different D2H analyses        //
// - event plane resolution                                      //
// - <pt> calculation with side band subtraction                 //
// - tracklet multiplicity calculation                            //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliAODMCParticle;
class AliAODMCHeader;
class AliGenEventHeader;
class AliAODEvent;
class TProfile;
class TClonesArray;
class TH1F;
class TH2F;
class TF1;

class AliVertexingHFUtils : public TObject{
 public:
  AliVertexingHFUtils();
  AliVertexingHFUtils(Int_t k);
  virtual ~AliVertexingHFUtils() {};

  // Significance calculator
  static void ComputeSignificance(Double_t signal, Double_t  errsignal, Double_t  background, Double_t  errbackground, Double_t &significance,Double_t &errsignificance);

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
  static TString  GetGenerator(Int_t label, AliAODMCHeader* header); 
  // Functions for tracklet multiplcity calculation
  void SetEtaRangeForTracklets(Double_t mineta, Double_t maxeta){
    fMinEtaForTracklets=mineta; 
    fMaxEtaForTracklets=maxeta;
  }
  static Int_t GetNumberOfTrackletsInEtaRange(AliAODEvent* ev, Double_t mineta, Double_t maxeta);
  Int_t GetNumberOfTrackletsInEtaRange(AliAODEvent* ev) const {
    return GetNumberOfTrackletsInEtaRange(ev,fMinEtaForTracklets,fMaxEtaForTracklets);
  }
  static Int_t GetGeneratedMultiplicityInEtaRange(TClonesArray* arrayMC, Double_t mineta, Double_t maxeta);
  static Int_t GetGeneratedPrimariesInEtaRange(TClonesArray* arrayMC, Double_t mineta, Double_t maxeta);
  static Int_t GetGeneratedPhysicalPrimariesInEtaRange(TClonesArray* arrayMC, Double_t mineta, Double_t maxeta);

  // Functions for computing average pt 
  static void AveragePt(Float_t& averagePt, Float_t& errorPt, Float_t ptmin, Float_t ptmax, TH2F* hMassD, Float_t massFromFit, Float_t sigmaFromFit, TF1* funcB2, Float_t sigmaRangeForSig=2.5, Float_t sigmaRangeForBkg=4.5, Int_t rebin=4);

  // Functions for computing true impact parameter of D meson
  static Double_t GetTrueImpactParameterDzero(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partDp);
  static Double_t GetTrueImpactParameterDplus(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partDp);

  static Double_t GetCorrectedNtracklets(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult);

 private:

  Int_t fK;             // ratio of measured harmonic to event plane harmonic
  Double_t fSubRes;     // sub-event resolution = sqrt(<cos[n(phiA-phiB)] >)
  Double_t fMinEtaForTracklets; // min eta for counting tracklets
  Double_t fMaxEtaForTracklets; // min eta for counting tracklets

  ClassDef(AliVertexingHFUtils,0) 
};
#endif
