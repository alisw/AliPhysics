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

  /*
  enum {
    kEPV0A, ///< EP form V0A
    kEPTPC, ///< EP form TPC
    kEPrandom ///< Random EP
  };

  void SetEPtype ( Int_t epType = kEPV0A ) { fEPtype = epType; }
  */

 private:

  AliAnalysisTaskFlowSingleMu(const AliAnalysisTaskFlowSingleMu&);
  AliAnalysisTaskFlowSingleMu& operator=(const AliAnalysisTaskFlowSingleMu&);

  TArrayD GetCentralityRange(TString sRange);

  /*
  enum {
    kTrackContainer, ///< CF container for tracks
    kHistoEP,      ///< Event plane distribution
    kNobjectTypes    ///< Number of objects
  };
  */
  
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
