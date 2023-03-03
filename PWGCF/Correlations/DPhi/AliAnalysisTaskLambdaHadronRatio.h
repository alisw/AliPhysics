/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full copyright notice */
/*
                AliAnalysisTaskLambdaHadronRatio class
                This task is for determining the lambda/hadron (pion) ratio
                in different kinematic regions using a two-particle azimuthal
                correlation method
                Origin: Ryan Hannigan, January 2021, ryan.hannigan@cern.ch
*/

#ifndef AliAnalysisTaskLambdaHadronRatio_H
#define AliAnalysisTaskLambdaHadronRatio_H

// Includes order: standard, ROOT, AliRoot (for objects in this file only)
#include "TString.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskSE.h"

// Forward declarations order: standard, ROOT, AliRoot (for pointers in this file only)
class TList;
class TH1D;
class TH1F;
class THnSparse;

class AliAODEvent;
class AliEventPoolManager;
class AliPIDResponse;
class AliMultSelection;
class AliEventPool;
class AliAODTrack;
class AliAODv0;


class AliAnalysisTaskLambdaHadronRatio : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskLambdaHadronRatio();
  AliAnalysisTaskLambdaHadronRatio(const char *name);
  virtual ~AliAnalysisTaskLambdaHadronRatio();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void LoadEfficiencies(TString filePath);

  void SetMultBounds(float multLow, float multHigh);
  void SetTriggerBit(float trigBit);
  void SetAssociatedBit(float assocBit);
  void SetCentEstimator(TString estimator);
  void SetPIDCuts(float nSigmaTPC_proton, float nSigmaTOF_proton, float nSigmaTPC_pion, float nSigmaTOF_pion, bool tofVeto);

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
  float fTPCnSigmaProtonCut; // TPC n sigma cut for protons
  float fTOFnSigmaProtonCut; // TOF n sigma cut for protons (veto if fTOFVeto is true, else just cut)
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
  TH2D* fLambdaEff_0_20; ///> lambda efficiency (0-20% mult)
  
  TH2D* fTriggerEff_20_50; ///> trigger efficiency (20-50% mult)
  TH2D* fAssociatedEff_20_50; ///> associated efficiency (20-50% mult)
  TH2D* fLambdaEff_20_50; ///> lambda efficiency (20-50% mult)
  
  TH2D* fTriggerEff_50_80; ///> trigger efficiency (50-80% mult)
  TH2D* fAssociatedEff_50_80; ///> associated efficiency (50-80% mult)
  TH2D* fLambdaEff_50_80; ///> lambda efficiency (50-80% mult)

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection

  TH2D* fTPCnSigmaProton; //!>! TPC n sigma for protons (vs pt)
  TH2D* fTPCnSigmaPion; //!>! TPC n sigma for pions (vs pt)
  TH2D* fTOFnSigmaProton; //!>! TOF n sigma for protons (vs pt)
  TH2D* fTOFnSigmaPion; //!>! TOF n sigma for pions (vs pt)

  TH2D* fTOFvTPCnSigmaProton; //!>! TPC vs TOF n sigma for protons (TPC is X axis, TOF is Y axis)
  TH2D* fTOFvTPCnSigmaPion; //!>! TPC vs TOF n sigma for pions (TPC is X axis, TOF is Y axis)

  THnSparse* fTriggerDist;  //!>! single particle trigger dist (corrected for efficiency)
  THnSparse* fTriggerDist_highestPt;  //!>! single particle trigger dist (corrected for efficiency, highest pt between 4 and 8)
  THnSparse* fAssociatedHDist;  //!>! single particle associated hadron dist

  THnSparse* fLambdaDist;  //!>! single particle lambda dist (eff corrected)
  THnSparse* fTriggeredLambdaDist;  //!>! single particle lambda dist within a triggered event (eff corrected)

  THnSparse* fDphiHLambda;  //!>! hadron-lambda correlation hist (efficiency corrected)
  THnSparse* fDphiHLambda_highestPt;  //!>! hadron-lambda correlation hist (efficiency corrected, highest pt trigger between 4 and 8)
  THnSparse* fDphiHH;   //!>! hadron-hadron correlation hist (efficiency corrected)
  THnSparse* fDphiHH_highestPt;   //!>! hadron-hadron correlation hist (efficiency corrected, highest pt trigger between 4 and 8)
  THnSparse* fDphiHLambdaMixed_highestPt; //!>! hadron-lambda mixed correlation hist (highest pt trigger between 4 and 8)
  THnSparse* fDphiHLambdaMixed; //!>! hadron-lambda mixed correlation hist (highest pt trigger between 4 and 8)
  THnSparse* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparse* fDphiHHMixed_highestPt; //!>! hadron-hadron mixed correlation hist (highest pt trigger between 4 and 8)


  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, double multPercentile, THnSparse* fDist, bool trig_eff=false);
  void FillMotherDist(std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist, bool isAntiLambda, bool lambda_eff=true);
  void MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff, bool isAntiLambda);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff=true);
  void MakeMixedHLambdaCorrelations(AliEventPool *fPool, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff, bool isAntiLambda);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx, double multPercentile, bool eff=true);
  bool PassDaughterCuts(AliAODTrack *track);
  bool PassTriggerCuts(AliAODTrack *track);
  bool PassAssociatedCuts(AliAODTrack *track);

  ClassDef(AliAnalysisTaskLambdaHadronRatio, 3);

};
#endif