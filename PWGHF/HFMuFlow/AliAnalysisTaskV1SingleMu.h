#ifndef ALIANALYSISTASKV1SINGLEMU_H
#define ALIANALYSISTASKV1SINGLEMU_H

//
// AliAnalysisTaskV1SingleMu
// Analysis task for v1 of single muons in the spectrometer with the scalar product method
//
//  Author: Audrey Francisco and Jacopo Margutti from AliAnalysisDimu
//


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisMuonUtility.h"
// #include "AliAODEvent.h"
// #include "AliESDEvent.h"
#include "AliMuonTrackCuts.h"
#include "AliMuonEventCuts.h"
// #include "AliMuonTrackCuts.h"

class TObjArray;
class THnSparse;
class AliMergeableCollection;
class AliUtilityMuonAncestor;


class AliAnalysisTaskV1SingleMu : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskV1SingleMu();
  AliAnalysisTaskV1SingleMu(const char *name);//, const AliMuonTrackCuts& cuts);
  virtual ~AliAnalysisTaskV1SingleMu();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void NotifyRun();
  void Terminate(Option_t *option);

  /// Get muon event cuts
  AliMuonEventCuts* GetMuonEventCuts() { return &fMuonEventCuts; }
  /// Get muon track cuts
  AliMuonTrackCuts* GetMuonTrackCuts() { return &fMuonTrackCuts; }

  /// Set muon event cuts
  void SetMuonEventCuts ( AliMuonEventCuts* muonEventCuts ) { fMuonEventCuts = *muonEventCuts; }
  /// Set muon track cuts
  void SetMuonTrackCuts ( AliMuonTrackCuts* muonTrackCuts ) { fMuonTrackCuts = *muonTrackCuts; }

  void UseZDCCalibRuns ( Bool_t b ) { fZdcCalibOnly = b; }

  void SetNPtBins  ( Int_t nBins ) { fNPtBins=nBins; }
  void SetNEtaBins ( Int_t nBins ) { fNEtaBins=nBins; }
  void SetNSPBins  ( Int_t nBins ) { fNSPBins=nBins; }
  void SetNCentBins( Int_t nBins ) { fNCentBins=nBins; }
  void SetNPhiBins ( Int_t nBins ) { fNPhiBins=nBins; }


  enum {
    kStepReconstructed,  ///< Reconstructed tracks
    kStepGeneratedMC,    ///< Generated tracks (MC)
    kNsteps              ///< Number of steps
  };
  enum {
    kHCentrality,    ///< event centrality
    kHvarPt,         ///< Pt at vertex
    kHvarEta,        ///< Pseudo-Rapidity
    kHvarCharge,     ///< Particle charge
    kHvarPhi,        ///< Phi
    kHvarV1QA,
    kHvarV1QB,
    kNvars           ///< THnSparse dimensions
  };

 private:
  TObject* GetMergeableObject ( TString identifier, TString objectName );

  AliAnalysisTaskV1SingleMu(const AliAnalysisTaskV1SingleMu&);
  AliAnalysisTaskV1SingleMu& operator=(const AliAnalysisTaskV1SingleMu&);

  Int_t GetParticleType ( AliVParticle* track );

  Bool_t IsZDCCalibrated(Int_t run);

  //add enum for MC (AliUtilityMuonAncestor)
    enum {
    kCharmMu,       ///< Mu from charm
    kBeautyMu,      ///< Mu from beauty
    kQuarkoniumMu,  ///< Mu from resonance
    kWbosonMu,      ///< Mu from W
    kZbosonMu,      ///< Mu from Z
    kDecayMu,       ///< Decay mu
    kSecondaryMu,   ///< Secondary mu
    kRecoHadron,    ///< Reconstructed hadron
    kUnidentified,  ///< Particle that fails matching kine
    kNtrackSources  ///< Total number of track sources
  };

  AliMuonEventCuts fMuonEventCuts;  ///< Muon event cuts
  AliMuonTrackCuts fMuonTrackCuts;  ///< Muon event cuts
  AliUtilityMuonAncestor* fUtilityMuonAncestor; ///< Utility to get the muon ancestor for MC
  Int_t fNPtBins;
  Int_t fNEtaBins;
  Int_t fNSPBins;
  Int_t fNCentBins;
  Int_t fNPhiBins;
  Int_t fHarmonic;
  AliMergeableCollection* fMergeableCollection; //!<! collection of mergeable objects
  Bool_t fZdcCalibOnly;
  THnSparse* fSparse; ///< CF container

  ClassDef(AliAnalysisTaskV1SingleMu, 3); // Single muon v1 analysis
};

#endif
