#ifndef ALIANALYSISTASKTRACKINGEFFPID
#define ALIANALYSISTASKTRACKINGEFFPID

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskTrackingEffPID
// AliAnalysisTaskSE to compute tracking and PID efficiencies for 
//  different particle species
//
// Authors:
//          M. Puccio
//          F. Prino
//          
//*************************************************************************


#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliPID.h"
#include "AliEventCuts.h"
#include <THnSparse.h>

class TList;

class AliAnalysisTaskTrackingEffPID : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskTrackingEffPID();
  virtual ~AliAnalysisTaskTrackingEffPID();

  static bool  HasTOF(AliVTrack *t);
  bool ConvertAndSelectAODTrack(AliAODTrack* aTrack, const AliESDVertex vESD, Double_t magField);

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  void SetTrackCuts(AliESDtrackCuts* cuts){
    if(fTrackCuts) delete fTrackCuts;
    fTrackCuts = new AliESDtrackCuts(*cuts);
  }
  AliESDtrackCuts* GetTrackCuts() const{
    return fTrackCuts;
  }
  void SetPrimarySelectionOption(Int_t opt){
    fPrimarySelectionOpt=opt;
  }
  void SetTrackletMultiplicityEstimator(){ fMultEstimator=0;}
  void SetVertexContribMultiplicityEstimator(){ fMultEstimator=1;}
  void SetTracksMultiplicityEstimator(){ fMultEstimator=2;}

  void SetCollisionSystem(TString collsy){
    collsy.ToLower();
    collsy.ReplaceAll("-","");
    if(collsy=="pbpb" || collsy=="xexe") fIsAA=kTRUE;
    else fIsAA=kFALSE;
  }
  void SetFilterBitCutForAODTracks(int fb){
    fFilterBit=fb;
  }
  void UseOnlyFilterBitCutForAODTracks(){
    fUseTrackCutsForAOD=kFALSE;    
  }
  void UseTrackCutObjectForAODTracks(){
    fUseTrackCutsForAOD=kTRUE;    
  }

  void SetUseGeneratedKine(bool flag) {fUseGeneratedKine=flag;}

  void SetTriggerMask(ULong64_t mask){
    fEventCut.OverrideAutomaticTriggerSelection(mask);
  }
  void SetUseSPDPileup(bool multDep, int nContrCut=5, double dzCut=0.8){
    fEventCut.OverridePileUpCuts(nContrCut,dzCut,3.,2.,5.,kTRUE);
    fEventCut.fUseMultiplicityDependentPileUpCuts=multDep;
  }
  void SetUseMVPileup(bool flag) {fEventCut.fPileUpCutMV=flag;}

  void SetUseImpactParameter(bool flag) {fUseImpPar=flag;}
  void SetPtHardRange(double pmin, double pmax){
    fSelectPtHardRange=kTRUE; fMinPtHard=pmin; fMaxPtHard=pmax;
  }
  void SelectedGeneratorName(TString name){
    fGenerToKeep=name.Data(); fSelectOnGenerator=kTRUE;}
  void ExcludedGeneratorName(TString name){
    fGenerToExclude=name.Data(); fSelectOnGenerator=kTRUE;}
  void KeepOnlyInjectedParticles(bool opt){fKeepOnlyInjected=opt;}
  void KeepOnlyUnderlyingEventParticles(bool opt){fKeepOnlyUE=opt;}
  TString GetGenerator(int label, TList *lh);
  bool IsInjectedParticle(int lab, TList *lh);
  AliEventCuts  fEventCut;


private:
  AliAnalysisTaskTrackingEffPID (const AliAnalysisTaskTrackingEffPID &source);
  AliAnalysisTaskTrackingEffPID &operator=(const AliAnalysisTaskTrackingEffPID &source);

  bool fUseTrackCutsForAOD;       /// flag to switch off/on fTrackCuts for AOD
  bool fUseGeneratedKine;         /// flag to use the generated pt, eta phi
  int  fPrimarySelectionOpt;      /// 0=no selection, 1=IsPhysicalPrimary, 2= cut on the origin
  int  fMultEstimator;            /// multiplicity estimator: 0=trackelts, 1=ITS+TPCtracks, 2=primary vertex contributors
  bool fIsAA;                     /// flag to control collision system
  int  fFilterBit;                /// filter-bit selection for AOD tracks
  AliESDtrackCuts* fTrackCuts;                            /// cut object
  bool fSelectOnGenerator;       /// flag to select events with generator name
  TString fGenerToKeep;          /// generator name to analyse
  TString fGenerToExclude;       /// generator name to exclude
  bool fKeepOnlyInjected;        /// flag to keep only injected particles
  bool fKeepOnlyUE;              /// flag to keep only underlying event
  bool fUseImpPar;               /// flag to enable plots vs. impact parameter
  bool fSelectPtHardRange;       /// flag to enable the cut on pthard
  double fMinPtHard;             /// min pthard
  double fMaxPtHard;             /// max pthard
  TList* fOutputList;                                     //!<! Output list
  TList* fListCuts;                                       //!<! Output with cuts
  TH1F*  fHistNEvents;                                    //!<!  histo with N of events  
  TH1D*  fHistNParticles;                                    //!<!  histo with N of particles
  TH1D*  fHistNTracks;                                    //!<!  histo with N of tracks
  TH1D*  hHistXsecVsPtHard;                               //!<!  control plot
  THnSparseF* fGenerated[AliPID::kSPECIESC][2];           //!<! Generated particles (pt, eta, phi, mult, zvert)
  THnSparseF* fGeneratedEvSel[AliPID::kSPECIESC][2];      //!<! Generated particles after event selection
  THnSparseF* fReconstructed[AliPID::kSPECIESC][2];       //!<! Reconstructed particles (pt, eta, phi, mult, zvert)
  THnSparseF* fReconstructedTOF[AliPID::kSPECIESC][2];    //!<! Reconstructed particles after PID (pt, eta, phi, mult, zvert)
  THnSparseF* fReconstructedPID[AliPID::kSPECIESC][2];    //!<! Reconstructed particles after PID (pt, eta, phi, mult, zvert)


  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskTrackingEffPID, 8);
  /// \endcond
};


#endif 
