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
#ifndef AliAnalysisTaskCorrPbPbMC_H
#define AliAnalysisTaskCorrPbPbMC_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TTree.h"
#include "AliEventCuts.h"

class TList;
class TTree;

class AliESDEvent;
//class AliAODEvent;
class AliESDtrackCuts;
//class AliESDcascade;
class AliMCEvent;
class AliStack;

class TH1D;
class TH2D;
class TH3D;
class TProfile;
class AliPIDResponse;
class AliMultSelection;


class AliAnalysisTaskCorrPbPbMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCorrPbPbMC();
  AliAnalysisTaskCorrPbPbMC(const char *name);
  virtual ~AliAnalysisTaskCorrPbPbMC();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t GetEvent();
  Bool_t PassedTrackQualityCuts (AliAODTrack *track);
  Bool_t KaonSelector (AliVTrack *track, Double_t nSigmaCut); 
  Bool_t ProtonSelector (AliVTrack *track, Double_t nSigmaCut); 
  Bool_t PionSelector (AliVTrack *track, Double_t nSigmaCut); 
  Bool_t PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type);
  Bool_t PassedSingleParticlePileUpCuts(AliAODTrack *track);
  void GetMCEffCorrectionHist();

  void SetListForTrkCorr (TList *fList)
  {
 	this->fListTRKCorr = (TList*) fList->Clone();
  }
  
  void SetVzRangeMax(Double_t VzMax)
  {
    this->fVertexZMax = VzMax;
  }
  void SetPileupCutValue (Int_t cutval)
  {
    this->fPileupCutVal = cutval;
  }
  void SetCentralityEstimator (Int_t val)
  {
    this->fCentralityEstimator_flag = val;
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
  void SetPIDnSigmaCut(Double_t PIDnSigmaCut_pion, Double_t PIDnSigmaCut_kaon, Double_t PIDnSigmaCut_proton)
  {
    this->fPIDnSigmaPionCut = PIDnSigmaCut_pion;
    this->fPIDnSigmaKaonCut = PIDnSigmaCut_kaon;
    this->fPIDnSigmaProtonCut = PIDnSigmaCut_proton;
  }
  void SetMinNoTPCCrossedRows(Double_t tpccrossedrows)
  {
    this->fTPCcrossedrows = tpccrossedrows;
  }
  
  
 private:
  
  
  AliESDEvent *fESDevent;
  AliAODEvent *fAODevent;
  AliMCEvent *fMCevent;
  AliStack *fMCstack;
  AliVEvent *fInputEvent;
  AliPIDResponse   *fPIDResponse;
  AliESDtrackCuts  *fESDtrackCuts;
  AliESDtrackCuts  *fESDtrackCuts_primary;
  AliAnalysisUtils *fUtils;
  AliEventCuts fAODeventCuts;
  TList *fOutputList;
  TList *fQAList;
  TList *fTreeList;
  Double_t fMultLow;
  Double_t fMultHigh;
  TH1D *hNumberOfEvents;
  TH1D *hNumberOfKaonPlus;
  TH1D *hNumberOfKaonMinus;
  TH1D *hNumberOfPionPlus;
  TH1D *hNumberOfPionMinus;
  TH1D *hNumberOfProtonPlus;
  TH1D *hNumberOfProtonMinus;

  //ftreeEvent object variables
  TTree *fTreeEvent;
  Float_t fTreeVariableCentrality;
  //generated
  Float_t fNoGenKaonPlus_ptmax2;
  Float_t fNoGenKaonMinus_ptmax2;
  Float_t fNoGenKaonPlus_ptmax3;
  Float_t fNoGenKaonMinus_ptmax3;
  Float_t fNoGenPionPlus_ptmax2;
  Float_t fNoGenPionMinus_ptmax2;
  Float_t fNoGenPionPlus_ptmax3;
  Float_t fNoGenPionMinus_ptmax3;
  Float_t fNoGenProtonPlus_ptmax2;
  Float_t fNoGenProtonMinus_ptmax2;
  Float_t fNoGenProtonPlus_ptmax3;
  Float_t fNoGenProtonMinus_ptmax3;
  //reconstruced
  Float_t fNoKaonPlus_ptmax2;
  Float_t fNoKaonMinus_ptmax2;
  Float_t fNoKaonPlus_ptmax3;
  Float_t fNoKaonMinus_ptmax3;
  Float_t fNoPionPlus_ptmax2;
  Float_t fNoPionMinus_ptmax2;
  Float_t fNoPionPlus_ptmax3;
  Float_t fNoPionMinus_ptmax3;
  Float_t fNoProtonPlus_ptmax2;
  Float_t fNoProtonMinus_ptmax2;
  Float_t fNoProtonPlus_ptmax3;
  Float_t fNoProtonMinus_ptmax3;
  //CORRECTED
  Float_t fCorrectedNoKaonPlus_ptmax2;
  Float_t fCorrectedNoKaonMinus_ptmax2;
  Float_t fCorrectedNoKaonPlus_ptmax3;
  Float_t fCorrectedNoKaonMinus_ptmax3;
  Float_t fCorrectedNoPionPlus_ptmax2;
  Float_t fCorrectedNoPionMinus_ptmax2;
  Float_t fCorrectedNoPionPlus_ptmax3;
  Float_t fCorrectedNoPionMinus_ptmax3;
  Float_t fCorrectedNoProtonPlus_ptmax2;
  Float_t fCorrectedNoProtonMinus_ptmax2;
  Float_t fCorrectedNoProtonPlus_ptmax3;
  Float_t fCorrectedNoProtonMinus_ptmax3;
  //Eff.square factors
  Float_t fEffSqrFactrPionMinus_ptmax2;
  Float_t fEffSqrFactrPionPlus_ptmax2;
  Float_t fEffSqrFactrProtonMinus_ptmax2;
  Float_t fEffSqrFactrProtonPlus_ptmax2;
  Float_t fEffSqrFactrKaonMinus_ptmax2;
  Float_t fEffSqrFactrKaonPlus_ptmax2;
  Float_t fEffSqrFactrPionMinus_ptmax3;
  Float_t fEffSqrFactrPionPlus_ptmax3;
  Float_t fEffSqrFactrProtonMinus_ptmax3;
  Float_t fEffSqrFactrProtonPlus_ptmax3;
  Float_t fEffSqrFactrKaonMinus_ptmax3;
  Float_t fEffSqrFactrKaonPlus_ptmax3;
  

  //Histograms for purity check
  TH1F *hist_KaonPlusWithoutPdg;
  TH1F *hist_KaonPlusWithPdg;
  TH1F *hist_KaonMinusWithoutPdg;
  TH1F *hist_KaonMinusWithPdg;
  TH1F *hist_PionPlusWithoutPdg;
  TH1F *hist_PionPlusWithPdg;
  TH1F *hist_PionMinusWithoutPdg;
  TH1F *hist_PionMinusWithPdg;
  TH1F *hist_ProtonPlusWithoutPdg;
  TH1F *hist_ProtonPlusWithPdg;
  TH1F *hist_ProtonMinusWithoutPdg;
  TH1F *hist_ProtonMinusWithPdg;

  //Generated histograms
  TH1F *hist_GenPionPlus;
  TH1F *hist_GenPionMinus;
  TH1F *hist_GenKaonPlus;
  TH1F *hist_GenKaonMinus;
  TH1F *hist_GenProtonPlus;
  TH1F *hist_GenProtonMinus;

  //Histograms for calculating centrality wise efficiency
  TH2D *f2Dhist_GenPionPlus;
  TH2D *f2Dhist_GenPionMinus;
  TH2D *f2Dhist_GenKaonPlus;
  TH2D *f2Dhist_GenKaonMinus;
  TH2D *f2Dhist_GenProtonPlus;
  TH2D *f2Dhist_GenProtonMinus;

  TH2D *f2Dhist_RecPionPlus;
  TH2D *f2Dhist_RecPionMinus;
  TH2D *f2Dhist_RecKaonPlus;
  TH2D *f2Dhist_RecKaonMinus;
  TH2D *f2Dhist_RecProtonPlus;
  TH2D *f2Dhist_RecProtonMinus;

  Double_t fVertexZMax;
  Int_t fFBNo;
  Double_t fChi2TPC;
  Double_t fChi2ITS;
  Double_t fPIDnSigmaPionCut;
  Double_t fPIDnSigmaKaonCut;
  Double_t fPIDnSigmaProtonCut;
  Double_t fTPCcrossedrows;
  
  //Efficiency list of histograms
  TList *fListTRKCorr; 
  TH1D *fHistMCEffKaonPlus;
  TH1D *fHistMCEffKaonMinus;
  TH1D *fHistMCEffPionPlus;
  TH1D *fHistMCEffPionMinus;
  TH1D *fHistMCEffProtonPlus;
  TH1D *fHistMCEffProtonMinus;
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
  TH2D *f2Dhist_nSigmaTPC_pion;
  TH2D *f2Dhist_nSigmaTPC_kaon;
  TH2D *f2Dhist_nSigmaTPC_proton;
  TH2D *f2Dhist_nSigmaTOF_pion;
  TH2D *f2Dhist_nSigmaTOF_kaon;
  TH2D *f2Dhist_nSigmaTOF_proton;
  TH2D *f2Dhist_nSigmaTPCplusTOF_pion;
  TH2D *f2Dhist_nSigmaTPCplusTOF_kaon;
  TH2D *f2Dhist_nSigmaTPCplusTOF_proton;

  //Pileup cut val
  Int_t fPileupCutVal;
  
  //Flag to select which centrality estimator
  Int_t fCentralityEstimator_flag;
 
  //Efficienc hist per centralities
  TH1D *fEffPionPlus[9];
  TH1D *fEffKaonPlus[9];
  TH1D *fEffProtonPlus[9];
  TH1D *fEffPionMinus[9];
  TH1D *fEffKaonMinus[9];
  TH1D *fEffProtonMinus[9];
  

  /*
  //Custom Functions:
  TF1 *fCenCutLowPU;
  TF1 *fCenCutHighPU;
  TF1 *fSPDCutPU;
  TF1 *fV0CutPU;
  TF1 *fMultCutPU;

  //Argument variables
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TString fTreeName;
  Double_t fEtaMax;
  Double_t fVertexZMax;
  Double_t fDCAxyMax;
  Double_t fDCAzMax;
  Double_t fChi2TPC;
  Double_t fChi2ITS;
  Double_t fNCrossedRowsTPC;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
 
  
  AliAnalysisTaskCorrPbPbMC(const AliAnalysisTaskCorrPbPbMC&);
  AliAnalysisTaskCorrPbPbMC& operator=(const AliAnalysisTaskCorrPbPbMC&);  
  ClassDef(AliAnalysisTaskCorrPbPbMC, 1);
};

#endif
