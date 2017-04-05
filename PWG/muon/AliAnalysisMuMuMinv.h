#ifndef ALIANALYSISMUMUMINV_H
#define ALIANALYSISMUMUMINV_H

/**
 *
 * \class AliAnalysisMuMuNch
 * \brief Invariant mass dimuon analysis
 * \author L. Aphecetche, J. Martin Blanco and B. Audurier (Subatech)
 */

#include "AliAnalysisMuMuBase.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliMergeableCollection.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TH2.h"

class TH2F;
class AliVParticle;
class TLorentzVector;
class AliMergeableCollectionProxy;

class AliAnalysisMuMuMinv : public AliAnalysisMuMuBase
{
public:

  AliAnalysisMuMuMinv(TH2* AccEffHisto=0x0, Int_t systLevel=0);
  virtual ~AliAnalysisMuMuMinv();

  Bool_t IsPtInRange(const AliVParticle& t1, const AliVParticle& t2,
                           Double_t& ptmin, Double_t& ptmax) const;

  void NameOfIsPtInRange(TString& name, Double_t& ymin, Double_t& ymax) const;

  Bool_t IsRapidityInRange(const AliVParticle& t1, const AliVParticle& t2) const;
  void NameOfIsRapidityInRange(TString& name) const { name = "PAIRY"; }

  Bool_t ShouldCorrectDimuonForAccEff() { return (fAccEffHisto != 0x0); }

  void FillMeanPtHisto() { fComputeMeanPt=kTRUE; }

  void SetMCptCut(Double_t mcptmin, Double_t mcptmax) { fmcptcutmin=mcptmin;fmcptcutmax=mcptmax; }

  void SetMuonWeight() { fWeightMuon=kTRUE; }

  void SetLegacyBinNaming() { fMinvBinSeparator = ""; }

  void SetBinsToFill(const char* particle, const char* bins);

  // create the original function with the parameters used in simulation to generate the pT distribution
  void SetOriginPtFunc(TString formula, const Double_t *param, Double_t xMin, Double_t xMax);
  // create the new function with its initial parameters to fit the generated/weighted pT distribution
  void SetNewPtFunc(TString formula, const Double_t *param, Double_t xMin, Double_t xMax);
  // create the original function with the parameters used in simulation to generate the y distribution
  void SetOriginYFunc(TString formula, const Double_t *param, Double_t xMin, Double_t xMax);
  // create the new function with its initial parameters to fit the generated/weighted y distribution
  void SetNewYFunc(TString formula, const Double_t *param, Double_t xMin, Double_t xMax);

  void DefineMinvRange(Double_t minvMin, Double_t minvMax, Double_t minvBinSize);

protected:

  void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                 const char* centrality,Bool_t =kFALSE);

  virtual void FillHistosForPair(const char* eventSelection,const char* triggerClassName,
                                 const char* centrality,
                                 const char* pairCutName,
                                 const AliVParticle& part,
                                 const AliVParticle& part2,
                                 const Bool_t IsMixedHisto);

  void FillHistosForMCEvent(const char* eventSelection,const char* triggerClassName,const char* centrality);

  void FillMinvHisto(TString* minvName,TProfile* hprof,TProfile* hprof2,AliMergeableCollectionProxy* proxy, TLorentzVector* pair4Momentum, Double_t inputWeight);

private:

  void CreateMinvHistograms(const char* eventSelection, const char* triggerClassName, const char* centrality);

  // normalize the function to its integral in the given range
  void NormFunc(TF1 *f, Double_t min, Double_t max);

  TString GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, Bool_t accEffCorrected, Double_t PairCharge=0, Bool_t mix =kFALSE) const;

  Double_t GetAccxEff(Double_t pt,Double_t rapidity);

  Double_t WeightMuonDistribution(Double_t pt);

  Double_t WeightPairDistribution(Double_t pt,Double_t rapidity);

  Double_t TriggerLptApt(Double_t *x, Double_t *par);

  Bool_t  CheckBinRangeCut(AliAnalysisMuMuBinning::Range* r, TLorentzVector* pair4Momentum, AliMergeableCollectionProxy* proxy);

  Bool_t CheckMCTracksMatchingStackAndMother(Int_t labeli, Int_t labelj, AliVParticle* mcTracki, AliVParticle* mcTrackj, Double_t inputWeightMC);

private:
  Bool_t fComputeMeanPt;
  Bool_t fWeightMuon;
  TH2F     * fAccEffHisto;
  TString fMinvBinSeparator;
  Int_t fsystLevel;
  TF1      *fPtFuncOld;              ///< original generated pT function with original parameters
  TF1      *fPtFuncNew;              ///< new generated pT fit function with new parameters
  TF1      *fYFuncOld;               ///< original generated y function with original parameters
  TF1      *fYFuncNew;               ///< new generated y fit function with new parameters
  TObjArray* fBinsToFill;
  Double_t fMinvBinSize;
  Double_t fMinvMin;
  Double_t fMinvMax;
  Double_t fmcptcutmin;
  Double_t fmcptcutmax;

  ClassDef(AliAnalysisMuMuMinv,8) // implementation of AliAnalysisMuMuBase for muon pairs
};

#endif
