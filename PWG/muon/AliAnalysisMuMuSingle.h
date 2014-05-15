#ifndef ALIANALYSISMUMUSINGLE_H
#define ALIANALYSISMUMUSINGLE_H

/**
 *
 * \class AliAnalysisMuMuSingle
 * 
 * \brief Histogramming of single muon tracks.
 *
 * \author L. Aphecetche (Subatech)
 *
 */

#include "AliAnalysisMuMuBase.h"

class AliMuonTrackCuts;
class TH2F;
class TObjArray;

class AliAnalysisMuMuSingle : public AliAnalysisMuMuBase
{
public:
  
  AliAnalysisMuMuSingle();
  virtual ~AliAnalysisMuMuSingle();
  
  virtual void ShouldSeparatePlusAndMinus(Bool_t value) { fShouldSeparatePlusAndMinus = value; }
  
  virtual Bool_t ShouldSeparatePlusAndMinus() const { return fShouldSeparatePlusAndMinus; }

  AliMuonTrackCuts* MuonTrackCuts();
  
  void SetMuonTrackCuts(const AliMuonTrackCuts& trackCuts);
  
  Bool_t IsRabsOK(const AliVParticle& part) const;
  void NameOfIsRabsOK(TString& name) const { name = "RABS"; }

  Bool_t IsEtaInRange(const AliVParticle& part, Double_t& etamin, Double_t& etamax) const;
  void NameOfIsEtaInRange(TString& name, Double_t& etamin, Double_t& etamax) const
  { name.Form("ETA%3.1f-%3.1f",etamin,etamax); }

  void SetRun(const AliInputEventHandler* eventHandler);

protected:
  
  void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                 const char* centrality);

  virtual void FillHistosForTrack(const char* eventSelection, const char* triggerClassName,
                                  const char* centrality,
                                  const char* trackCutName,
                                  const AliVParticle& part);

private:
  
  void CreateTrackHisto(const char* eventSelection,
                        const char* triggerClassName,
                        const char* centrality,
                        const char* hname, const char* htitle,
                        Int_t nbinsx, Double_t xmin, Double_t xmax,
                        Int_t nbinsy=-1, Double_t ymin=0.0, Double_t ymax=0.0,
                        Bool_t separatePlusAndMinus=kFALSE) const;

  Double_t GetTrackTheta(const AliVParticle& particle) const;
  
  /* methods prefixed with EA should really not exist at all. They are there
   only because the some of our base interfaces are shamelessly incomplete or
   inadequate...
   */
  
  Int_t EAGetNumberOfMuonTracks() const;
  
//  Int_t EAGetNumberOfSelectMuonTracks() const;
  
  Double_t EAGetTrackDCA(const AliVParticle& particle) const;
     
private:
  
  /// not implemented on purpose
  AliAnalysisMuMuSingle& operator=(const AliAnalysisMuMuSingle& rhs);
  /// not implemented on purpose
  AliAnalysisMuMuSingle(const AliAnalysisMuMuSingle& rhs);
  
  AliMuonTrackCuts* fMuonTrackCuts; //! common cuts for muon tracks (from Diego)
  Bool_t fShouldSeparatePlusAndMinus; // whether or not to histogram mu+ and mu- separately
  TH2F* fAccEffHisto; // dimuon acc x eff (y vs pt)
  
  ClassDef(AliAnalysisMuMuSingle,1) // implementation of AliAnalysisMuMuBase for single mu analysis
};

#endif
