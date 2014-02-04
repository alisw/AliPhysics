#ifndef ALIANALYSISMUMUMINV_H
#define ALIANALYSISMUMUMINV_H

/**
 * \class AliAnalysisMuMuMinv
 * \brief Invariant mass analysis of muon pairs
 * \author L. Aphecetche (Subatech)
 */

#include "AliAnalysisMuMuBase.h"
#include "AliAnalysisMuMuBinning.h"
#include "TString.h"

class TH2F;
class AliVParticle;

class AliAnalysisMuMuMinv : public AliAnalysisMuMuBase
{
public:

  AliAnalysisMuMuMinv();
  virtual ~AliAnalysisMuMuMinv();
  
  Bool_t IsRapidityInRange(const AliVParticle& t1, const AliVParticle& t2,
                           Double_t& ymin, Double_t& ymax) const;

  virtual void FillHistosForMCEvent(const char* /*eventSelection*/,const char* /*triggerClassName*/,const char* /*centrality*/);

  void NameOfIsRapidityInRange(TString& name, Double_t& ymin, Double_t& ymax) const;

  Bool_t ShouldCorrectDimuonForAccEff() const { return (fAccEffHisto != 0x0); }
  
  void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                 const char* centrality);

  virtual void FillHistosForPair(const char* eventSelection,const char* triggerClassName,
                                 const char* centrality,
                                 const char* pairCutName,
                                 const AliVParticle& part,
                                 const AliVParticle& part2);
  
private:
  
  /// not implemented on purpose
  AliAnalysisMuMuMinv(const AliAnalysisMuMuMinv& rhs);
  /// not implemented on purpose
  AliAnalysisMuMuMinv& operator=(const AliAnalysisMuMuMinv& rhs);
  
  void CreateMinvHistograms(const char* eventSelection, const char* triggerClassName, const char* centrality);
  
  TString GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, Bool_t accEffCorrected) const;
  

private:
  TH2F* fAccEffHisto; // input Acc x Eff (optional)
  
  ClassDef(AliAnalysisMuMuMinv,1) // implementation of AliAnalysisMuMuBase for muon pairs
};

#endif

