#ifndef ALIANALYSISMUMUMINV_H
#define ALIANALYSISMUMUMINV_H

/**
 *
 * \class AliAnalysisMuMuNch
 * \brief Invariant mass dimuon analysis
 * \author L. Aphecetche and J. Martin Blanco (Subatech)
 */

#include "AliAnalysisMuMuBase.h"
#include "AliAnalysisMuMuBinning.h"
#include "TString.h"
#include "TH2.h"

class TH2F;
class AliVParticle;

class AliAnalysisMuMuMinv : public AliAnalysisMuMuBase
{
public:

  AliAnalysisMuMuMinv(TH2* AccEffHisto=0x0, Bool_t computeMeanPt=kFALSE, Int_t systLevel=0);
  virtual ~AliAnalysisMuMuMinv();
  
  Bool_t IsPtInRange(const AliVParticle& t1, const AliVParticle& t2,
                           Double_t& ptmin, Double_t& ptmax) const;
  
  void NameOfIsPtInRange(TString& name, Double_t& ymin, Double_t& ymax) const;
  
  Bool_t IsRapidityInRange(const AliVParticle& t1, const AliVParticle& t2) const;
  void NameOfIsRapidityInRange(TString& name) const { name = "PAIRY"; }
  
  Bool_t ShouldCorrectDimuonForAccEff() { return (fAccEffHisto != 0x0); }
  
  void SetLegacyBinNaming() { fMinvBinSeparator = ""; }
  
protected:
  
  void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                 const char* centrality);

  virtual void FillHistosForPair(const char* eventSelection,const char* triggerClassName,
                                 const char* centrality,
                                 const char* pairCutName,
                                 const AliVParticle& part,
                                 const AliVParticle& part2);
  
  void FillHistosForMCEvent(const char* eventSelection,const char* triggerClassName,const char* centrality);
  
private:
  
  void CreateMinvHistograms(const char* eventSelection, const char* triggerClassName, const char* centrality);
  
  TString GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, Bool_t accEffCorrected) const;
  
  Double_t GetAccxEff(Double_t pt,Double_t rapidity);
  
  Double_t WeightDistribution(Double_t pt,Double_t rapidity);
  
  Double_t powerLaw3Par(Double_t *x, Double_t *par);
  
  Double_t normPol12Par(Double_t *x, Double_t *par);

private:
  Bool_t fcomputeMeanPt;
  TH2F* fAccEffHisto;
  TString fMinvBinSeparator;
  
  Int_t fsystLevel;
  
  ClassDef(AliAnalysisMuMuMinv,3) // implementation of AliAnalysisMuMuBase for muon pairs
};

#endif
