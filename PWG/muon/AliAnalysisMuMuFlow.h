#ifndef ALIANALYSISMUMUFLOW_H
#define ALIANALYSISMUMUFLOW_H

/**
 *
 * \class AliAnalysisMuMuEP
 * \brief Elliptic dimuon analysis with the event plane method
 * \author A. Francisco (Subatech)
 */

#include "AliAnalysisMuMuBase.h"
#include "AliAnalysisMuMuBinning.h"
#include "TString.h"
#include "TH2.h"
#include "TVector2.h"

class TH2F;
class AliVParticle;

class AliAnalysisMuMuFlow : public AliAnalysisMuMuBase
{
public:

  AliAnalysisMuMuFlow(TH2* AccEffHisto=0x0, Int_t systLevel=0);
  virtual ~AliAnalysisMuMuFlow();

  Bool_t IsQnInRange(const AliVEvent& event, Double_t& qnMin, Double_t& qnMax) const;
  void NameOfIsQnInRange(TString& name, Double_t& qnmin, Double_t& qnmax) const;

  Bool_t IsDPhiInPlane(const AliVParticle& t1, const AliVParticle& t2) const;
  Bool_t IsDPhiOutOfPlane(const AliVParticle& t1, const AliVParticle& t2) const;
  void NameOfIsDPhiInPlane(TString& name) const;
  void NameOfIsDPhiOutOfPlane(TString& name) const;

  Bool_t ShouldCorrectDimuonForAccEff() { return (fAccEffHisto != 0x0); }

  void SetMuonWeight() { fWeightMuon=kTRUE; }

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

private:

  void CreateMinvHistograms(const char* eventSelection, const char* triggerClassName, const char* centrality);

  // normalize the function to its integral in the given range
  void NormFunc(TF1 *f, Double_t min, Double_t max);

  TString GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, Bool_t accEffCorrected) const;

  Double_t GetAccxEff(Double_t pt,Double_t rapidity);

  Double_t WeightMuonDistribution(Double_t pt);

  Double_t WeightPairDistribution(Double_t pt,Double_t rapidity);

  Double_t TriggerLptApt(Double_t *x, Double_t *par);

  Double_t GetEventPlane(const char* detector, Int_t step = 3) const;

  TVector2 GetQn(const char* detector, Int_t step = 3) const;

private:
  Bool_t fcomputeMeanV2; //mv2 with EP method
  Bool_t fcomputeEP; //dist in dphi bins
  Bool_t fcomputeSP; //scalar product
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
  Int_t fNDetectors;
  Int_t fHar;
  TString fEqSteps  [5] = {"raw", "plain", "rec", "align","twist"};
  TString fDetectors[3] = {"SPD","VZEROA", "VZEROC"};

  ClassDef(AliAnalysisMuMuFlow,1) // implementation of AliAnalysisMuMuBase for muon pairs
};

#endif
