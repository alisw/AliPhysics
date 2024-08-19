/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2_H
#define AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TTree.h"
#include "AliEventCuts.h"
#include "TExMap.h"

class TList;
class TTree;

class AliESDEvent;
//class AliAODEvent;
class AliESDtrackCuts;
class AliMCEvent;
class AliStack;
class AliVTrack;
class AliAODTrack;
class AliAODv0;

class TH1D;
class TH2D;
class TH3D;
class TProfile;
class AliPIDResponse;
class AliPIDCombined;
class AliMultSelection;


class AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2();
  AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2(const char *name);
  virtual ~AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t GetEvent();
  Bool_t PassedTrackQualityCuts (AliAODTrack *track);
  Bool_t KaonSelector (AliVTrack *track, Double_t nSigmaCut); 
  Bool_t ProtonSelector (AliVTrack *track, Double_t nSigmaCut); 
  Bool_t PionSelector (AliVTrack *track, Double_t nSigmaCut); 
  Bool_t ElectronRejectionCut(AliVTrack *track, Int_t fCut);
  Bool_t PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type, Double_t PIDcut);
  Bool_t PassedSingleParticlePileUpCuts(AliAODTrack *track);
  void GetMCEffCorrectionHist();
  void FilltrackQAplots_beforeCut(Double_t fDcaXY, Double_t fDcaZ, Double_t fEta, Double_t fITSchi2perNcls, Double_t fTPCchi2perNcls, Double_t fTPCcrossedrows);
  void FilltrackQAplots_afterCut(Double_t fDcaXY, Double_t fDcaZ, Double_t fEta, Double_t fITSchi2perNcls, Double_t fTPCchi2perNcls, Double_t fTPCcrossedrows);
  void FillPIDQAplots_beforeCut(AliVTrack *track);
  void FillPIDQAplots_afterCut(AliVTrack *track, Bool_t Pionflag, Bool_t Kaon_flag, Bool_t Proton_flag);
  void GlobalTracksAOD(AliAODEvent *aAOD);
  Bool_t PassedDaughterTrackDCAtoVertexSelectionCutsV0 (AliAODv0 *v0);
  Bool_t PassedV0SelectionTopologicalCutsForLambda (AliAODv0 *v0);
  Bool_t PassedV0SelectionTopologicalCutsForK0s (AliAODv0 *v0);
  Int_t CheckFlagITSHitOrTOFhit(AliAODTrack *track, Double_t lMagField);
  Bool_t IsLambdaCandidate (AliAODv0 *V0, AliAODTrack *pos, AliAODTrack *neg, Double_t PIDcut, Double_t centrality, Double_t LambdaMassCut);
  Bool_t IsAntiLambdaCandidate (AliAODv0 *V0, AliAODTrack *pos, AliAODTrack *neg, Double_t PIDcut, Double_t centrality, Double_t LambdaMassCut);
  Bool_t IsK0sCandidate (AliAODv0 *V0, AliAODTrack *pos, AliAODTrack *neg, Double_t PIDcut, Double_t centrality,Double_t K0sMassCut);
  Bool_t HasTrackPIDTPC(AliVTrack *track);
  Bool_t HasTrackPIDTOF(AliVTrack *track);
  Int_t IdentifyTrackBayesian(AliVTrack *track);

  void SetListForTrkCorr (TList *fList)
  {
    this->fListTRKCorr = (TList*) fList->Clone();
  }
  void SetVzRangeMax(Double_t VzMax)
  {
    this->fVertexZMax = VzMax;
  }
  void SetTrackFilterBit(Int_t FBno)
  {
    this->fFBNo = FBno;
  }
  void SetMaxChi2PerTPCClusterRange(Double_t chi2tpc)
  {
    this->fChi2TPC = chi2tpc;
  }
  void SetMaxChi2PerITSClusterRange(Double_t chi2its)
  {
    this->fChi2ITS = chi2its;
  }
  void SetPIDnSigmaCut(Double_t PIDnSigmaCut)
  {
    this->fPIDnSigmaCut = PIDnSigmaCut;
  }
  void SetEtaLeftCut(Double_t etaleft)
  {
    this->fEtaLeftCut = etaleft;
  }
  void SetEtaCut(Double_t eta_min)
  {
    this->fEtaMin = eta_min;
  }
  void SetMinNCrossedRowsTPC(Double_t nCrossedRow)
  {
    this->fTPCcrossedrows = nCrossedRow;
  }
  void SetTreeName(TString TreeName)
  {
    this->fTreeName = TreeName;
  }
  void SetCentralityEstimator (Int_t val)
  {
    this->fCentralityEstimator_flag = val;
  }
  void SetPileupCutValue (Int_t cutval)
  {
    this->fPileupCutVal = cutval;
  }
  void SetEfficiencyEffectImposeFlag (Int_t flagg)
  {
    this->fEffFlag = flagg;
  }
  void SetEfficiencyCorrectionFlag (Int_t flagval)
  {
    this->fEffCorrectionFlag = flagval;
  }
  void SetRejectElectronFlag (Int_t elrej_flag)
  {
    this->fRejectElectron_cut = elrej_flag;
  }
  void SetExclusivePIDCutFlag (Int_t flag1)
  {
    this->fExclusivePIDCut_flag = flag1;
  }
  void SetFillTrackQAHistogramsFlag (Int_t flag2)
  {
    this->fFillTrackQAhists_flag = flag2;
  }
  void SetFillTrackPIDQAHistogramsFlag (Int_t flag3)
  {
    this->fFillPIDhists_flag = flag3;
  }
  void SetSelectPiKaPrByBayesianPIDFlag (Int_t flagg)
  {
    this->fBayesianPID_flag = flagg;
  }
  void SetMinV0TracksTpcClustersNo (Int_t nClusters)
  {
    this->fMinV0TracksTpcClustersNo = nClusters; //Default: 80
  }
  void SetMinV0TracksTpcCrossedRowsNo (Int_t nCrossedRows)
  {
    this->fMinV0TracksTpcCrossedRowsNo = nCrossedRows; //Default: 80
  }
  void SetMaxV0TracksTpcSharedClustersNo (Int_t nSharedClusters)
  {
    this->fMaxV0TracksTpcSharedClustersNo = nSharedClusters; //Default: 5
  }
  void SetMaxV0TracksChi2TPCperClstr (Double_t chi2)
  {
    this->fMaxV0TracksChi2TPCperClstr = chi2; //Default: 2.5
  }
  void SetMaxV0TracksChi2ITSperClstr (Double_t chi2)
  {
    this->fMaxV0TracksChi2ITSperClstr = chi2; //Default: 36
  }
  void SetRatio_TPCcrossedrowsByFindableclusters (Double_t ratio)
  {
    this->fRatioTPCcrossedrowsByFindableclusters = ratio;  //Default: 0.8
  }
  void SetDcaV0DaughterTracksToPV (Double_t DcaDauTracksToPV)
  {
    this->fDcaV0DaughterTracksToPV = DcaDauTracksToPV; //Default: 0.1
  }
  void SetLambdaPropLifetime (Double_t propLifetime)
  {
    this->fLambdaPropLifetime = propLifetime; //Default: 40 cm
  }
  void SetLambdaTransDecayRadius_Min(Double_t min)
  {
    this->fMinLambdaTransDecayRadius = min;  //Default: 5cm
  }
  void SetLambdaTransDecayRadius_Max(Double_t max)
  {
    this->fMaxLambdaTransDecayRadius = max;  //Default: 100cm
  }
  void SetLambdaDCAv0daughters (Double_t DcaV0daughters)
  {
    this->fLambdaDcaV0daughters = DcaV0daughters;  //Default: 0.2
  }
  void SetLambdaDcaV0toPV (Double_t DcaV0toPV)
  {
    this->fLambdaDcaV0toPV = DcaV0toPV; //Default: 0.5
  }
  void SetLambdaCosPA (Double_t cosPAval)
  {
    this->fLambdaCosPAval = cosPAval; //Default: 0.997
  }
  void SetK0sPropLifetime (Double_t propLifetime)
  {
    this->fK0sPropLifetime = propLifetime; //Default: 10 cm
  }
  void SetK0sTransDecayRadius_Min(Double_t min)
  {
    this->fMinK0sTransDecayRadius = min;  //Default: 5cm
  }
  void SetK0sTransDecayRadius_Max(Double_t max)
  {
    this->fMaxK0sTransDecayRadius = max;  //Default: 100cm
  }
  void SetK0sDCAv0daughters (Double_t DcaV0daughters)
  {
    this->fK0sDcaV0daughters = DcaV0daughters;  //Default: 0.4
  }
  void SetK0sCosPA (Double_t cosPAval)
  {
    this->fK0sCosPAval = cosPAval; //Default: 0.999
  }
  void SetK0sArmentousCut (Double_t cutval)
  {
    this->fArmentousCutVal = cutval; //Default: 0.2
  }
  void SetLambdaDaughtersPIDcut (Double_t PIDcut)
  {
    this->fLambdaDaughtersPIDcut = PIDcut;
  }
  void SetK0sDaughtersPIDcut (Double_t PIDcut)
  {
    this->fK0sDaughtersPIDcut = PIDcut;
  }
  void SetLambdaFinalMassCut (Double_t cutval)
  {
    this->fLambdaMassCut = cutval; //Default val = 0.003
  }
  void SetK0sFinalMassCut (Double_t cutval)
  {
    this->fK0sMassCut = cutval; //Default val = 0.015
  }
  void SetBayesPIDPionVal (Double_t bayesPidPi)
  {
    this->fPIDbayesPion = bayesPidPi;
  }
  void SetBayesPIDKaonVal (Double_t bayesPidKa)
  {
    this->fPIDbayesKaon = bayesPidKa;
  }
  void SetBayesPIDProtonVal (Double_t bayesPidPr)
  {
    this->fPIDbayesProton = bayesPidPr;
  }
  
 private:
  
  
  AliESDEvent *fESDevent;
  AliAODEvent *fAODevent;
  AliVEvent *fInputEvent;
  AliPIDResponse   *fPIDResponse;
  AliPIDCombined   *fPIDCombined;
  AliESDtrackCuts  *fESDtrackCuts;
  AliESDtrackCuts  *fESDtrackCuts_primary;
  AliAnalysisUtils *fUtils;
  AliEventCuts fAODeventCuts;
  TList *fOutputList;
  TList *fQAList;
  Double_t fMultLow;
  Double_t fMultHigh;
  TH1D *hNumberOfEvents;
  TH1D *hNumberOfKaonEtaLess0;
  TH1D *hNumberOfPionEtaLess0;
  TH1D *hNumberOfProtonEtaLess0;
  TH1D *hNumberOfLambdaEtaLess0;
  TH1D *hNumberOfK0sEtaLess0;

  //ftreeEvent object variables
  TTree *fTreeEvent;
  Float_t fTreeVariableCentrality;
  Float_t fPtsum_hadrons_less0;
  Float_t fPtsum_hadrons_greaterEtaMin;
  Float_t fNsum_hadrons_less0;
  Float_t fNsum_hadrons_greaterEtaMin;
  Float_t fNsum_pions_less0;
  Float_t fNsum_kaons_less0;
  Float_t fNsum_protons_less0;
  Float_t fNsum_lambdas_less0;
  Float_t fNsum_K0s_less0;
  Float_t fPt_no_hadron[20];
  Float_t fPt_no_pion[20];
  Float_t fPt_no_kaon[20];
  Float_t fPt_no_proton[20];
  Float_t fPt_no_lambda[20];
  Float_t fPt_no_K0s[20];
  
  //Efficiency list of histograms
  TList *fListTRKCorr; 
  TH1D *fHistMCEffKaonPlus;
  TH1D *fHistMCEffKaonMinus;
  TH1D *fHistMCEffPionPlus;
  TH1D *fHistMCEffPionMinus;
  TH1D *fHistMCEffProtonPlus;
  TH1D *fHistMCEffProtonMinus;
  TH1D *fHistMCEffHadronPlus;
  TH1D *fHistMCEffHadronMinus;
  TH1D *fEffPionPlus[9];
  TH1D *fEffKaonPlus[9];
  TH1D *fEffProtonPlus[9];
  TH1D *fEffPionMinus[9];
  TH1D *fEffKaonMinus[9];
  TH1D *fEffProtonMinus[9];
  
  //Track QA histograms
  TH1D *hist_beforeCut_DCAxy;
  TH1D *hist_beforeCut_DCAz;
  TH1D *hist_beforeCut_eta;
  TH1D *hist_beforeCut_chi2perTPCclstr;
  TH1D *hist_beforeCut_chi2perITSclstr;
  TH1D *hist_beforeCut_TPCncrossedrows;
  TH1D *hist_afterCut_DCAxy;
  TH1D *hist_afterCut_DCAz;
  TH1D *hist_afterCut_eta;
  TH1D *hist_afterCut_chi2perTPCclstr;
  TH1D *hist_afterCut_chi2perITSclstr;
  TH1D *hist_afterCut_TPCncrossedrows;
  
  
  //TPC, TOF, TPC+TOF nSigma histograms
  //before cuts
  TH2D *f2Dhist_beforeCut_nSigmaTPC_pion;
  TH2D *f2Dhist_beforeCut_nSigmaTPC_kaon;
  TH2D *f2Dhist_beforeCut_nSigmaTPC_proton;
  TH2D *f2Dhist_beforeCut_nSigmaTOF_pion;
  TH2D *f2Dhist_beforeCut_nSigmaTOF_kaon;
  TH2D *f2Dhist_beforeCut_nSigmaTOF_proton;
  TH2D *f2Dhist_beforeCut_nSigmaTPCplusTOF_pion;
  TH2D *f2Dhist_beforeCut_nSigmaTPCplusTOF_kaon;
  TH2D *f2Dhist_beforeCut_nSigmaTPCplusTOF_proton;
  TH2D *f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_pion;
  TH2D *f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_kaon;
  TH2D *f2Dhist_beforeCut_nSigmaTPC_vs_nSigmaTOF_proton;
  TH2D *f2Dhist_beforeCut_TPCdEdx_all;
  TH2D *f2Dhist_beforeCut_TOFtime_all;
  
  //after cuts
  TH2D *f2Dhist_afterCut_nSigmaTPC_pion;
  TH2D *f2Dhist_afterCut_nSigmaTPC_kaon;
  TH2D *f2Dhist_afterCut_nSigmaTPC_proton;
  TH2D *f2Dhist_afterCut_nSigmaTOF_pion;
  TH2D *f2Dhist_afterCut_nSigmaTOF_kaon;
  TH2D *f2Dhist_afterCut_nSigmaTOF_proton;
  TH2D *f2Dhist_afterCut_nSigmaTPCplusTOF_pion;
  TH2D *f2Dhist_afterCut_nSigmaTPCplusTOF_kaon;
  TH2D *f2Dhist_afterCut_nSigmaTPCplusTOF_proton;
  TH2D *f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_pion;
  TH2D *f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_kaon;
  TH2D *f2Dhist_afterCut_nSigmaTPC_vs_nSigmaTOF_proton;
  TH2D *f2Dhist_afterCut_TPCdEdx_all;
  TH2D *f2Dhist_afterCut_TOFtime_all;
  TH2D *f2Dhist_afterCut_TPCdEdx_pion;
  TH2D *f2Dhist_afterCut_TOFtime_pion;
  TH2D *f2Dhist_afterCut_TPCdEdx_kaon;
  TH2D *f2Dhist_afterCut_TOFtime_kaon;
  TH2D *f2Dhist_afterCut_TPCdEdx_proton;
  TH2D *f2Dhist_afterCut_TOFtime_proton;
  TH3D *f3DhistMassLambdaAll_vs_Pt_beforeMasscut_Cent;
  TH3D *f3DhistMassK0s_vs_Pt_beforeMasscut_Cent;
  TH3D *f3DhistMassLambdaAll_vs_Pt_afterMasscut_Cent;
  TH3D *f3DhistMassK0s_vs_Pt_afterMasscut_Cent;
  
 
  Double_t fVertexZMax;
  Int_t fFBNo;
  Double_t fChi2TPC;
  Double_t fChi2ITS;
  Double_t fPIDnSigmaCut;
  Double_t fTPCcrossedrows;
  Int_t fCentralityEstimator_flag;
  Int_t fPileupCutVal;
  Double_t fEtaLeftCut;
  Double_t fEtaMin;
  Int_t fEffFlag;
  TString fTreeName;
  Int_t fEffCorrectionFlag;
  Int_t fExclusivePIDCut_flag;
  Int_t fRejectElectron_cut;
  Int_t fFillTrackQAhists_flag;
  Int_t fFillPIDhists_flag;
  Int_t fBayesianPID_flag;
  Int_t fMinV0TracksTpcClustersNo;
  Int_t fMinV0TracksTpcCrossedRowsNo;
  Int_t fMaxV0TracksTpcSharedClustersNo;
  Double_t fMaxV0TracksChi2TPCperClstr;
  Double_t fMaxV0TracksChi2ITSperClstr;
  Double_t fRatioTPCcrossedrowsByFindableclusters;
  Double_t fDcaV0DaughterTracksToPV;
  Double_t fLambdaPropLifetime;
  Double_t fMinLambdaTransDecayRadius;
  Double_t fMaxLambdaTransDecayRadius;
  Double_t fLambdaDcaV0daughters;
  Double_t fLambdaDcaV0toPV;
  Double_t fLambdaCosPAval;
  Double_t fK0sPropLifetime;
  Double_t fMinK0sTransDecayRadius;
  Double_t fMaxK0sTransDecayRadius;
  Double_t fK0sDcaV0daughters;
  Double_t fK0sCosPAval;
  Double_t fArmentousCutVal;
  Double_t fLambdaDaughtersPIDcut;
  Double_t fK0sDaughtersPIDcut;
  Double_t fLambdaMassCut;
  Double_t fK0sMassCut;
  Double_t fPIDbayesPion;
  Double_t fPIDbayesKaon;
  Double_t fPIDbayesProton;
  
  TExMap *fGlobalTracksAOD; //! global tracks in AOD for FB128 **Ante**


  
  AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2(const AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2&);
  AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2& operator=(const AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2&);  
  ClassDef(AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2, 1);
};

#endif
