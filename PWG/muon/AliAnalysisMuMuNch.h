#ifndef ALIANALYSISMUMUNCH_H
#define ALIANALYSISMUMUNCH_H

/**
 *
 * \class AliAnalysisMuMuNch
 * \brief Charged particle multiplicity analysis (for plotting muon variables against it)
 * \author L. Aphecetche and J. Martin-Bianco (Subatech)
 */

#include "AliAnalysisMuMuBase.h"

class TH2F;
class TH2;
class AliVEvent;
class TAxis;
class TF1;
class AliAODVertex;

class AliAnalysisMuMuNch : public AliAnalysisMuMuBase
{
public:
  
  AliAnalysisMuMuNch(TH2* spdCorrection=0x0, Double_t etaMin=-0.5, Double_t etaMax=0.5
                     , Double_t zmin=-40, Double_t zmax=40,Bool_t disableHistos=kFALSE,
                     Bool_t computeResolution=kFALSE);
  virtual ~AliAnalysisMuMuNch();
  
  Bool_t HasAtLeastNTrackletsInEtaRange(const AliVEvent& event, Int_t n,
                                        Double_t& etaMin, Double_t& etaMax) const;

  void NameOfHasAtLeastNTrackletsInEtaRange(TString& name, Int_t n,
                                            Double_t& etaMin, Double_t& etaMax) const;
  
  virtual void SetEvent(AliVEvent* event, AliMCEvent* mcEvent=0x0);

  virtual void Terminate(Option_t* /*opt*/="");

protected:

  void AddHisto(const char* eventSelection,
                const char* triggerClassName,
                const char* centrality,
                const char* histoname,
                Double_t z,
                TH1* h,
                Bool_t isMC=kFALSE);
  
  void AttachSPDAcceptance(UInt_t dataType,
                           const char* eventSelection,
                           const char* triggerClassName,
                           const char* centrality,const char* histoname);


  void FillHistosForMCEvent(const char* eventSelection, const char* triggerClassName,
              const char* centrality);
  

  void FillHistosForEvent(const char* eventSelection, const char* triggerClassName,
                          const char* centrality);
  
  virtual void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                         const char* centrality);

  void DefineSPDAcceptance();
  
  void GetEtaRangeSPD(Double_t spdZVertex, Double_t etaRange[]);

  Double_t GetSPDCorrection(Double_t zvert, Double_t eta) const;

private:
  TH2F* fSPDCorrection; // Nch/Tracklets_SPD (eta vs z)
  TAxis* fEtaAxis; // Eta axis used for the histos
  TAxis* fZAxis;  // Z vertex axis used for the histos
  AliVEvent* fCurrentEvent; //! cache of the current event
  Double_t fEtaMin; // Minimum tracklet eta value 
  Double_t fEtaMax; // Maximum tracklet eta value 
  Double_t fZMin; // Minimum z vertex value 
  Double_t fZMax; // Maximum z vertex value 
  Bool_t fResolution; // Flag to set the resolution computation
  
  TF1* fSPD1LR; // SPD acceptance shape
  TF1* fSPD1LL; // SPD acceptance shape
  TF1* fSPD2LR; // SPD acceptance shape
  TF1* fSPD2LL; // SPD acceptance shape
  
  ClassDef(AliAnalysisMuMuNch,1) // implementation of AliAnalysisMuMuBase for Nch analysis
};

#endif
