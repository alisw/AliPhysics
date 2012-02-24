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
#include "TVector3.h"

class TObjArray;
class AliHistogramCollection;
class TString;
class TAxis;
class AliVParticle;
class AliAODEvent;
class AliMuonTrackCuts;

class AliAnalysisTaskSingleMu : public AliVAnalysisMuon {
 public:
  AliAnalysisTaskSingleMu();
  AliAnalysisTaskSingleMu(const char *name, const AliMuonTrackCuts& cuts);
  virtual ~AliAnalysisTaskSingleMu();

  virtual void   Terminate(Option_t *option);

  void MyUserCreateOutputObjects();
  void ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality);
  
  /// Set minimum number of vertex contributors
  void SetMinNvtxContributors(Int_t minNvtxContributors) { fMinNvtxContirbutors = minNvtxContributors; }

 private:

  AliAnalysisTaskSingleMu(const AliAnalysisTaskSingleMu&);
  AliAnalysisTaskSingleMu& operator=(const AliAnalysisTaskSingleMu&);

  // Histograms to extract average DCA position
  enum {
    kIPVz,           ///< Interaction point vertex distribution
    kTrackContainer, ///< CF container for tracks
    kNhistoTypes     ///< Number of histograms
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
  
  Int_t fMinNvtxContirbutors;  ///< Minimum number of vertex contributors

  TObjArray* fThetaAbsKeys;    ///< Name of theta at absorber end

  ClassDef(AliAnalysisTaskSingleMu, 2); // Single muon analysis
};

#endif
