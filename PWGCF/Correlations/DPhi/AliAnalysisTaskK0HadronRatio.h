/* Copyright(c) 1998-2024, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full copyright notice */
/*
                AliAnalysisTaskK0HadronRatio class
                This task is for determining the k0/hadron (pion) ratio
                in different kinematic regions using a two-particle azimuthal
                correlation method
                Origin: Ryan Hannigan, August 2024, ryan.hannigan@cern.ch
*/

#ifndef AliAnalysisTaskK0HadronRatio_H
#define AliAnalysisTaskK0HadronRatio_H

// Includes order: standard, ROOT, AliRoot (for objects in this file only)
#include "TString.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskSE.h"

// Forward declarations order: standard, ROOT, AliRoot (for pointers in this file only)
class TList;
class TH1D;
class TH2D;
class TH1F;
class THnSparse;

class AliAODEvent;
class AliEventPoolManager;
class AliPIDResponse;
class AliMultSelection;
class AliEventPool;
class AliAODTrack;
class AliAODv0;


class AliAnalysisTaskK0HadronRatio : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskK0HadronRatio();
  AliAnalysisTaskK0HadronRatio(const char *name);
  virtual ~AliAnalysisTaskK0HadronRatio();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void LoadEfficiencies(TString filePath);

  void SetMultBounds(float multLow, float multHigh);
  void SetTriggerBit(float trigBit);
  void SetAssociatedBit(float assocBit);
  void SetCentEstimator(TString estimator);
  void SetPIDCuts(float nSigmaTPC_pion, float nSigmaTOF_pion, bool tofVeto);

  struct AliMotherContainer {
    AliAODv0* vzero;
    int daughter1ID;
    int daughter2ID;
  };

private:
  float fMultLow; // lower bound for multiplicity
  float fMultHigh; // upper bound for multiplicity
  float fDaughterBit; // filter bit for daughter particle
  float fAssociatedBit; // filter bit for associated particle
  float fTriggerBit; // filter bit for trigger particle
  float fTPCnSigmaPionCut; // TPC n sigma cut for pions
  float fTOFnSigmaPionCut; // TOF n sigma cut for pions (veto if fTOFVeto is true, else just cut)
  bool fTOFVeto; // if true, use TOF n sigma cut as a veto, else just cut

  TString fCentEstimator;

  AliAODEvent* fAOD; //!>! input event
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fCorPoolMgr; //!>! correlation pool manager
  AliEventPoolManager *fCorPoolMgr_highestPt; //!>! correlation pool manager for highest pt trigger
  
  // for efficiency: x axis is pt, y axis is eta

  TH2D* fTriggerEff_0_20; ///> trigger efficiency (0-20% mult)
  TH2D* fAssociatedEff_0_20; ///> associated efficiency (0-20% mult)
  TH2D* fK0Eff_0_20; ///> k0 efficiency (0-20% mult)
  
  TH2D* fTriggerEff_20_50; ///> trigger efficiency (20-50% mult)
  TH2D* fAssociatedEff_20_50; ///> associated efficiency (20-50% mult)
  TH2D* fK0Eff_20_50; ///> k0 efficiency (20-50% mult)
  
  TH2D* fTriggerEff_50_80; ///> trigger efficiency (50-80% mult)
  TH2D* fAssociatedEff_50_80; ///> associated efficiency (50-80% mult)
  TH2D* fK0Eff_50_80; ///> k0 efficiency (50-80% mult)

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection

  TH1D* fMultDist; //!>! mult dist
  TH2D* fEventSelection; //!>! event selection hist (x axis bins are: minbias, pass nCont, pass zVertex, have trigger, have trigger + k0_candidate; y axis bins are: 0-80 mult. in increments of 10)

  TH2D* fTPCnSigmaPion; //!>! TPC n sigma for pions (vs pt)
  TH2D* fTOFnSigmaPion; //!>! TOF n sigma for pions (vs pt)

  TH2D* fTOFvTPCnSigmaPion; //!>! TPC vs TOF n sigma for pions (TPC is X axis, TOF is Y axis)

  THnSparse* fTriggerDist;  //!>! single particle trigger dist (corrected for efficiency)
  THnSparse* fTriggerDist_highestPt;  //!>! single particle trigger dist (corrected for efficiency, highest pt between 4 and 8)
  THnSparse* fAssociatedHDist;  //!>! single particle associated hadron dist

  THnSparse* fK0Dist;  //!>! single particle k0 dist (eff corrected)
  THnSparse* fTriggeredK0Dist;  //!>! single particle k0 dist within a triggered event (eff corrected)

  THnSparse* fDphiHK0;  //!>! hadron-k0 correlation hist (efficiency corrected)
  THnSparse* fDphiHK0_highestPt;  //!>! hadron-k0 correlation hist (efficiency corrected, highest pt trigger between 4 and 8)
  THnSparse* fDphiHH;   //!>! hadron-hadron correlation hist (efficiency corrected)
  THnSparse* fDphiHH_highestPt;   //!>! hadron-hadron correlation hist (efficiency corrected, highest pt trigger between 4 and 8)
  THnSparse* fDphiHK0Mixed; //!>! hadron-K0 mixed correlation hist
  THnSparse* fDphiHK0Mixed_highestPt; //!>! hadron-K0 mixed correlation hist (highest pt trigger between 4 and 8)
  THnSparse* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparse* fDphiHHMixed_highestPt; //!>! hadron-hadron mixed correlation hist (highest pt trigger between 4 and 8)

  THnSparse* fV0TopDist; //!>! dist containing v0 top variables to see which cuts are applied by default


  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, double multPercentile, THnSparse* fDist, bool trig_eff=false);
  void FillMotherDist(std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist, bool k0_eff=true);
  void MakeSameHK0Correlations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> k0_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff=true);
  void MakeMixedHK0Correlations(AliEventPool *fPool, std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> k0_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx, double multPercentile, bool eff=true);
  bool PassDaughterCuts(AliAODTrack *track);
  bool PassTriggerCuts(AliAODTrack *track);
  bool PassAssociatedCuts(AliAODTrack *track);

  ClassDef(AliAnalysisTaskK0HadronRatio, 3);

};
#endif
