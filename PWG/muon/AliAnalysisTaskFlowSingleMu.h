#ifndef ALIANALYSISTASKFLOWSINGLEMU_H
#define ALIANALYSISTASKFLOWSINGLEMU_H

/* $Id: AliAnalysisTaskFlowSingleMu.h 55545 2012-04-04 07:16:39Z pcrochet $ */ 

//
// AliAnalysisTaskFlowSingleMu
// Analysis task for flow of single muons in the spectrometer
//
//  Author: Diego Stocco
//

#include "AliVAnalysisMuon.h"
#include "TRandom3.h"

class TObjArray;
class TString;
class TArrayD;
class TAxis;
class AliMuonTrackCuts;

class AliAnalysisTaskFlowSingleMu : public AliVAnalysisMuon {
 public:
  AliAnalysisTaskFlowSingleMu();
  AliAnalysisTaskFlowSingleMu(const char *name, const AliMuonTrackCuts& cuts);
  virtual ~AliAnalysisTaskFlowSingleMu();

  virtual void   Terminate(Option_t *option);

  void MyUserCreateOutputObjects();
  void ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality);

 private:

  AliAnalysisTaskFlowSingleMu(const AliAnalysisTaskFlowSingleMu&);
  AliAnalysisTaskFlowSingleMu& operator=(const AliAnalysisTaskFlowSingleMu&);

  TArrayD GetCentralityRange(TString sRange);
  
  enum {
    kStepReconstructed,  ///< Reconstructed tracks
    kStepGeneratedMC,    ///< Generated tracks (MC)
    kNsteps              ///< Number of steps
  };  
  
  enum {
    kHvarPt,         ///< Pt at vertex
    kHvarEta,        ///< Pseudo-Rapidity
    kHvarPhi,        ///< Phi
    kHvarDeltaPhi,   ///< Phi_mu - Psi_plane
    kHvarCharge,     ///< Particle charge
    kHvarMotherType, ///< Mother type (MC only)
    kNvars           ///< THnSparse dimensions
  };

  TObjArray* fEPKeys; ///< EP keys
  TRandom3* fRandom; //!< Random number generator

  ClassDef(AliAnalysisTaskFlowSingleMu, 1); // Single muon analysis
};

#endif
