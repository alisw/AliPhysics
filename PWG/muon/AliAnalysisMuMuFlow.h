#ifndef ALIANALYSISMUMUFLOW_H
#define ALIANALYSISMUMUFLOW_H

/**
 *
 * \class AliAnalysisMuMuEP
 * \brief Elliptic dimuon analysis with the event plane or scalar product method
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

  AliAnalysisMuMuFlow(TH2* AccEffHisto=0x0, TList *q2Map=0x0, Int_t systLevel=0);
  virtual ~AliAnalysisMuMuFlow();

  Bool_t Isq2InSmallRange(const AliVEvent& event) const;
  Bool_t Isq2InLargeRange(const AliVEvent& event) const;
  void NameOfIsq2InSmallRange(TString& name) const;
  void NameOfIsq2InLargeRange(TString& name) const;

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

  void ComputeMeanV2(Bool_t b){fcomputeMeanV2=b;}
  void UseEventPlaneMethod(Bool_t b){fcomputeEP=b;}
  void UseScalarProductMethod(Bool_t b){fcomputeSP=b;}
  void UsSESEMethod(Bool_t b){fESE=b;}


protected:

  void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                 const char* centrality,Bool_t =kFALSE);

  virtual void FillHistosForPair(const char* eventSelection,const char* triggerClassName,
                                 const char* centrality,
                                 const char* pairCutName,
                                 const AliVParticle& part,
                                 const AliVParticle& part2,
                                 const Bool_t IsMixedHisto);

  void FillHistosForEvent(const char* eventSelection,const char* triggerClassName,const char* centrality);
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

  // TVector2 GetQn(Int_t detector) const { return Q2[detector]; }
  // void SetQn(Int_t detector, TVector2 Qn) { Q2[detector].Set(Qn); }
  // Double_t GetEventPlane(Int_t detector) const { return phiEP[detector]; }
  // void SetEventPlane(Int_t detector, Double_t psi) { phiEP[detector]=psi; }

  Double_t GetCentrality() const;

private:
  Bool_t fcomputeMeanV2; //mv2 with EP method
  Bool_t fcomputeEP; //dist in dphi bins
  Bool_t fcomputeSP; //scalar product
  Bool_t fESE; //event shape engineering
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
  Double_t EP[3];
  Double_t Q2[3][2];
  TString fEqSteps  [5] = {"raw", "plain", "rec", "align","twist"};
  TString fDetectors[3] = {"SPD","VZEROA", "VZEROC"};
  TH1F *fq2Map[2];

  ClassDef(AliAnalysisMuMuFlow,3) // implementation of AliAnalysisMuMuFlow for muon pairs
};

#endif
