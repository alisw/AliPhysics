#ifndef ALIANALYSISTASKSINGLEMU_H
#define ALIANALYSISTASKSINGLEMU_H

/* $Id$ */ 

//
// AliAnalysisTaskSingleMu
// Analysis task for single muons in the spectrometer
//
//  Author: Diego Stocco
//

#include "AliVAnalysisMuon.h"

class TObjArray;
class TString;
class TAxis;
class AliMuonTrackCuts;

class AliAnalysisTaskSingleMu : public AliVAnalysisMuon {
 public:
  AliAnalysisTaskSingleMu();
  AliAnalysisTaskSingleMu(const char *name, const AliMuonTrackCuts& cuts);
  virtual ~AliAnalysisTaskSingleMu();

  virtual void   Terminate(Option_t *option);

  void MyUserCreateOutputObjects();
  void ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality);
  
  /// Apply cut on dimuon invariant mass (to reject Z contribution)
  void SetCutDimu ( Bool_t cutOnDimu = kTRUE ) { fCutOnDimu = cutOnDimu; }

  /// Use associated MC kinematics for reconstructed tracks
  void SetUseMCKineForRecoTracks ( Bool_t useMCKineForRecoTracks = kTRUE ) { fUseMCKineForRecoTracks = useMCKineForRecoTracks; }
  
  enum {
    kIPVz,           ///< Interaction point vertex distribution
    kTrackContainer, ///< CF container for tracks
    kNobjectTypes    ///< Number of objects
  };
  
  enum {
    kThetaAbs23,  ///< Theta abs 2-3 deg
    kThetaAbs310, ///< Theta abs 3-10 deg
    kNthetaAbs    ///< Number of theta abs bins
  };
  
  enum {
    kStepReconstructed,  ///< Reconstructed tracks
    kStepGeneratedMC,    ///< Generated tracks (MC)
    kNsteps              ///< Number of steps
  };
  
  enum {
    kHvarPt,         ///< Pt at vertex
    kHvarEta,        ///< Pseudo-Rapidity
    kHvarPhi,        ///< Phi
    kHvarVz,         ///< Z vertex position
    kHvarCharge,     ///< Particle charge
    kHvarThetaAbs,   ///< Theta abs bin
    kHvarMotherType, ///< Mother type (MC only)
    kNvars           ///< THnSparse dimensions
  };

 private:

  AliAnalysisTaskSingleMu(const AliAnalysisTaskSingleMu&);
  AliAnalysisTaskSingleMu& operator=(const AliAnalysisTaskSingleMu&);

  TObjArray* fThetaAbsKeys;    ///< Name of theta at absorber end
  Bool_t fCutOnDimu;           ///< Cut on dimuons
  Bool_t fUseMCKineForRecoTracks; ///< Use MC kinematics for reconstructed tracks

  ClassDef(AliAnalysisTaskSingleMu, 5); // Single muon analysis
};

#endif
