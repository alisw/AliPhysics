#ifndef ALIANALYSISMUMUSINGLE_H
#define ALIANALYSISMUMUSINGLE_H

#include "AliAnalysisMuMuBase.h"

/**
 *
 * \class AliAnalysisMuMuSingle
 *
 * \brief Histogramming of single muon tracks.
 *
 * \author L. Aphecetche (Subatech)
 *
 */

#include "AliAnalysisMuonUtility.h"

class AliMergeableCollectionProxy;
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

  Bool_t IsPDCAOK(const AliVParticle& part);
  void NameOfIsPDCAOK(TString& name) const { name = "PDCA";}

  Bool_t IsMatchingTriggerAnyPt(const AliVParticle& part) const { return ( AliAnalysisMuonUtility::GetMatchTrigger(&part) >= 1 ); }
  void NameOfIsMatchingTriggerAnyPt(TString& name) const { name = "MATCHANY";}

  Bool_t IsMatchingTriggerLowPt(const AliVParticle& part) const { return ( AliAnalysisMuonUtility::GetMatchTrigger(&part) >= 2 ); }
  void NameOfIsMatchingTriggerLowPt(TString& name) const { name = "MATCHLOW";}

  Bool_t IsMatchingTriggerHighPt(const AliVParticle& part) const { return ( AliAnalysisMuonUtility::GetMatchTrigger(&part) >= 3 ); }
  void NameOfIsMatchingTriggerHighPt(TString& name) const { name = "MATCHHIGH";}

  Bool_t IsRabsOK(const AliVParticle& part) const;
  void NameOfIsRabsOK(TString& name) const { name = "RABS"; }

  Bool_t IsEtaInRange(const AliVParticle& part) const;
  void NameOfIsEtaInRange(TString& name) const
  { name = "ETA"; }

  void SetRun(const AliInputEventHandler* eventHandler);

  void MakePtEtaSpectraPerBunchCrossing() { fPtEtaSpectraPerBCX = kTRUE; }

  void MakeDCAHistos() { fDCAHistos = kTRUE; }

protected:

  void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                 const char* centrality,Bool_t mix =kFALSE);

  virtual void FillHistosForTrack(const char* eventSelection, const char* triggerClassName,
                                  const char* centrality,
                                  const char* trackCutName,
                                  const AliVParticle& part);

  void FillHistosForMuonTrack(AliMergeableCollectionProxy& proxy, const AliVParticle& track);


private:

  void CreateTrackHisto(const char* eventSelection,
                        const char* triggerClassName,
                        const char* centrality,
                        const char* hname, const char* htitle,
                        Int_t nbinsx, Double_t xmin, Double_t xmax,
                        Int_t nbinsy=-1, Double_t ymin=0.0, Double_t ymax=0.0,
                        Bool_t separatePlusAndMinus=kFALSE) const;

  void CreateTrackTHnSparse(const char* eventSelection,
                            const char* triggerClassName,
                            const char* centrality,
                            const char* hname, const char* htitle,
                            Int_t nDim, Int_t* nbinsx, Double_t* xmin, Double_t* xmax,
                            Bool_t separatePlusAndMinus) const;

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

  Bool_t fPtEtaSpectraPerBCX; // make pt vs eta spectra bunch by bunch (caution : much slower !)
  Bool_t fDCAHistos; // make DCA histograms

  ClassDef(AliAnalysisMuMuSingle,3) // implementation of AliAnalysisMuMuBase for single mu analysis
};

#endif
