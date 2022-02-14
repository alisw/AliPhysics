#ifndef ALIANALYSISTASKPTEMCALTRIGGER_H_
#define ALIANALYSISTASKPTEMCALTRIGGER_H_
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliAnalysisTaskEmcalJet.h"
#include "AliESDtrackCuts.h"
#include <TClonesArray.h>
#include <TList.h>
#include "AliCutValueRange.h"

class TArrayD;
class Axis;
class THistManager;
class AliClusterContainer;
class AliEmcalJet;
class AliESDtrack;
class AliJetContainer;
class AliParticleContainer;
class AliVTrack;
class AliVParticle;

namespace PWGJE {
  
namespace EMCALJetTasks {

/**
 * \class AliAnalysisTaskPtEMCalTrigger
 * \brief Old charged hadron analysis in EMCAL-triggered events
 *
 * Analysis task of the pt analysis on EMCal-triggered events
 */
class AliAnalysisTaskPtEMCalTrigger : public AliAnalysisTaskEmcalJet {
public:
  enum EEMCalTriggerType_t{
    kEMCalJetLow = 0,
    kEMCalJetHigh = 1,
    kEMCalGammaLow = 2,
    kEMCalGammaHigh = 3
  };
  static const Int_t kNJetRadii;
  static const Double_t *kJetRadii;

  AliAnalysisTaskPtEMCalTrigger();
  AliAnalysisTaskPtEMCalTrigger(const char *name);
  ~AliAnalysisTaskPtEMCalTrigger();

  virtual void UserCreateOutputObjects();
  virtual Bool_t Run();

  void AddESDTrackCuts(AliESDtrackCuts *trackCuts);
  void AddCutsForAOD(AliESDtrackCuts *trackCuts, UInt_t filterBits);
  void SetEtaRange(double etamin, double etamax) { fEtaRange.SetLimits(etamin, etamax); }
  void SetPtRange(double ptmin, double ptmax) { fPtRange.SetLimits(ptmin, ptmax); }
  void SetVertexRange(double vmin, double vmax) { fVertexRange.SetLimits(vmin, vmax); }
  void SetClusterEnergyRange(double emin, double emax) { fEnergyRange.SetLimits(emin,emax); }
  void SetSwapEta() { fSwapEta = kTRUE; }
  void UseTriggersFromTriggerMaker() { fUseTriggersFromTriggerMaker = kTRUE; }
  void AddJetContainerName(const Char_t * contname, Bool_t isMC = kFALSE);
  void SelectAllTracks(Bool_t doAll) { fSelectAllTracks = doAll; }

private:
  AliAnalysisTaskPtEMCalTrigger(const AliAnalysisTaskPtEMCalTrigger &);
  AliAnalysisTaskPtEMCalTrigger &operator=(const AliAnalysisTaskPtEMCalTrigger &);
  void CreateDefaultPtBinning(TArrayD &binning) const;
  void CreateDefaultZVertexBinning(TArrayD &binning) const;
  void CreateDefaultEtaBinning(TArrayD &binning) const;
  void DefineAxis(TAxis &axis, const char *name, const char *title, const TArrayD &binning, const char **labels = NULL);
  void DefineAxis(TAxis &axis, const char *name, const char *title, int nbins, double min, double max, const char **labels = NULL);
  void FillEventHist(const char *trigger, double vz, bool isPileup);
  void FillTrackHist(const char *trigger, const AliVTrack *track, double vz, bool isPileup, int cut, bool isMinBias, double jetradius = -1.);
  void FillClusterHist(const char *trigger, const AliVCluster *clust, double vz, bool isPileup, bool isMinBias);
  void FillMCParticleHist(const char *histname, const AliVParticle * const part, double vz, bool isPileup);
  bool IsTrueTrack(const AliVTrack *const) const;
  TString BuildTriggerString();
  const AliVVertex *GetSPDVertex() const;
  const AliEmcalJet *FoundTrackInJet(const AliVParticle * const track, AliJetContainer *const jets) const;
  const AliEmcalJet *FoundClusterInJet(const AliVCluster * const clust, AliJetContainer *const jets) const;
  bool TrackInJet(const AliVParticle *const track, const AliEmcalJet *reconstructedJet, const AliParticleContainer *const particles) const;
  bool ClusterInJet(const AliVCluster *const clust, const AliEmcalJet *reconstructedJet, const AliClusterContainer *const particles) const;
  bool IsInRadius(const AliVParticle *const track, const AliEmcalJet *reconstructedJet, Double_t radius) const;
  bool IsInRadius(const AliVCluster *const clust, const AliEmcalJet *reconstructedJet, Double_t radius) const;

  /// Histogram container for the task
  THistManager                  *fHistos;                 //!
  TList 						            *fListTrackCuts;		      ///< List of track cuts

  // Cuts
  AliCutValueRange<double>      fEtaRange;                ///< Eta Selection Range
  AliCutValueRange<double>	    fPtRange;				          ///< Pt Selection Range
  AliCutValueRange<double>      fEnergyRange;             ///< Cluster energy selection range
  AliCutValueRange<double>      fVertexRange;             ///< Vertex cut

  // Jet containers
  TList                         fJetContainersMC;          ///< List of jet containers for MC
  TList                         fJetContainersData;        ///< List of jet containers for Data

  // Settings
  Bool_t                        fSelectAllTracks;         ///< Loop over all tracks
  Bool_t						            fSwapEta;				          ///< Allow swapping of the eta sign in asymmetric collision systems
  Bool_t 						            fUseTriggersFromTriggerMaker; ///< Use trigger classes from trigger maker

  ClassDef(AliAnalysisTaskPtEMCalTrigger, 1);           // Analysis of EMCal triggered events
};

}

}
#endif
