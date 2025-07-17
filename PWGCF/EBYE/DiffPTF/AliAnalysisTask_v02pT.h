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
#ifndef AliAnalysisTask_v02pT_H
#define AliAnalysisTask_v02pT_H

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

class TH1D;
class TH2D;
class TH3D;
class TProfile;
class TProfile2D;
class AliPIDResponse;
class AliPIDCombined;
class AliMultSelection;


class AliAnalysisTask_v02pT : public AliAnalysisTaskSE {
 public:
  AliAnalysisTask_v02pT();
  AliAnalysisTask_v02pT(const char *name);
  virtual ~AliAnalysisTask_v02pT();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t GetEvent();
  Bool_t PassedTrackQualityCuts (AliAODTrack *track);
  Bool_t KaonSelector (AliVTrack *track, Double_t nSigmaCut); 
  Bool_t ProtonSelector (AliVTrack *track, Double_t nSigmaCut); 
  Bool_t PionSelector (AliVTrack *track, Double_t nSigmaCut); 
  Bool_t ElectronRejectionCut(AliVTrack *track, Int_t fCut);
  void GetMCEffCorrectionHist();
  void FilltrackQAplots_beforeCut(Double_t fDcaXY, Double_t fDcaZ, Double_t fEta, Double_t fITSchi2perNcls, Double_t fTPCchi2perNcls, Double_t fTPCcrossedrows);
  void FilltrackQAplots_afterCut(Double_t fDcaXY, Double_t fDcaZ, Double_t fEta, Double_t fITSchi2perNcls, Double_t fTPCchi2perNcls, Double_t fTPCcrossedrows);
  void FillPIDQAplots_beforeCut(AliVTrack *track);
  void FillPIDQAplots_afterCut(AliVTrack *track, Bool_t Pionflag, Bool_t Kaon_flag, Bool_t Proton_flag);
  void GlobalTracksAOD(AliAODEvent *aAOD);
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
  Float_t fQcosGapL_21;
  Float_t fQcosGapL_01;
  Float_t fQsinGapL_21;
  Float_t fQsinGapL_01;
  Float_t fQcosGapR_21;
  Float_t fQcosGapR_01;
  Float_t fQsinGapR_21;
  Float_t fQsinGapR_01;
  Float_t fPt_no_hadron[20];
  Float_t fPt_no_pion[20];
  Float_t fPt_no_kaon[20];
  Float_t fPt_no_proton[20];
  Float_t fPt_no_hadron_cos[20];
  Float_t fPt_no_pion_cos[20];
  Float_t fPt_no_kaon_cos[20];
  Float_t fPt_no_proton_cos[20];
  Float_t fPt_no_hadron_sin[20];
  Float_t fPt_no_pion_sin[20];
  Float_t fPt_no_kaon_sin[20];
  Float_t fPt_no_proton_sin[20];
  
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

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //analysis central profiles----------------------->
  
  TProfile *hProfile_X_real;
  TProfile *hProfile_X_imag;
  TProfile *hProfile_Y_real;
  TProfile *hProfile_Y_imag;
  TProfile *hProfile_XY_real;
  TProfile *hProfile_XY_imag;

  TProfile2D *hProfile2D_Z_hadrons;
  TProfile2D *hProfile2D_Z_pions;
  TProfile2D *hProfile2D_Z_kaons;
  TProfile2D *hProfile2D_Z_protons;

  TProfile2D *hProfile2D_YZ_real_hadrons;
  TProfile2D *hProfile2D_YZ_imag_hadrons;
  TProfile2D *hProfile2D_YZ_real_pions;
  TProfile2D *hProfile2D_YZ_imag_pions;
  TProfile2D *hProfile2D_YZ_real_kaons;
  TProfile2D *hProfile2D_YZ_imag_kaons;
  TProfile2D *hProfile2D_YZ_real_protons;
  TProfile2D *hProfile2D_YZ_imag_protons;

  TProfile2D *hProfile2D_XZ_real_hadrons;
  TProfile2D *hProfile2D_XZ_imag_hadrons;
  TProfile2D *hProfile2D_XZ_real_pions;
  TProfile2D *hProfile2D_XZ_imag_pions;
  TProfile2D *hProfile2D_XZ_real_kaons;
  TProfile2D *hProfile2D_XZ_imag_kaons;
  TProfile2D *hProfile2D_XZ_real_protons;
  TProfile2D *hProfile2D_XZ_imag_protons;

  TProfile2D *hProfile2D_XYZ_real_hadrons;
  TProfile2D *hProfile2D_XYZ_imag_hadrons;
  TProfile2D *hProfile2D_XYZ_real_pions;
  TProfile2D *hProfile2D_XYZ_imag_pions;
  TProfile2D *hProfile2D_XYZ_real_kaons;
  TProfile2D *hProfile2D_XYZ_imag_kaons;
  TProfile2D *hProfile2D_XYZ_real_protons;
  TProfile2D *hProfile2D_XYZ_imag_protons;

  //analysis error profiles----------------------->
  
  TProfile *hSubProfile_X_real[10];
  TProfile *hSubProfile_X_imag[10];
  TProfile *hSubProfile_Y_real[10];
  TProfile *hSubProfile_Y_imag[10];
  TProfile *hSubProfile_XY_real[10];
  TProfile *hSubProfile_XY_imag[10];

  TProfile2D *hSubProfile2D_Z_hadrons[10];
  TProfile2D *hSubProfile2D_Z_pions[10];
  TProfile2D *hSubProfile2D_Z_kaons[10];
  TProfile2D *hSubProfile2D_Z_protons[10];

  TProfile2D *hSubProfile2D_YZ_real_hadrons[10];
  TProfile2D *hSubProfile2D_YZ_imag_hadrons[10];
  TProfile2D *hSubProfile2D_YZ_real_pions[10];
  TProfile2D *hSubProfile2D_YZ_imag_pions[10];
  TProfile2D *hSubProfile2D_YZ_real_kaons[10];
  TProfile2D *hSubProfile2D_YZ_imag_kaons[10];
  TProfile2D *hSubProfile2D_YZ_real_protons[10];
  TProfile2D *hSubProfile2D_YZ_imag_protons[10];

  TProfile2D *hSubProfile2D_XZ_real_hadrons[10];
  TProfile2D *hSubProfile2D_XZ_imag_hadrons[10];
  TProfile2D *hSubProfile2D_XZ_real_pions[10];
  TProfile2D *hSubProfile2D_XZ_imag_pions[10];
  TProfile2D *hSubProfile2D_XZ_real_kaons[10];
  TProfile2D *hSubProfile2D_XZ_imag_kaons[10];
  TProfile2D *hSubProfile2D_XZ_real_protons[10];
  TProfile2D *hSubProfile2D_XZ_imag_protons[10];

  TProfile2D *hSubProfile2D_XYZ_real_hadrons[10];
  TProfile2D *hSubProfile2D_XYZ_imag_hadrons[10];
  TProfile2D *hSubProfile2D_XYZ_real_pions[10];
  TProfile2D *hSubProfile2D_XYZ_imag_pions[10];
  TProfile2D *hSubProfile2D_XYZ_real_kaons[10];
  TProfile2D *hSubProfile2D_XYZ_imag_kaons[10];
  TProfile2D *hSubProfile2D_XYZ_real_protons[10];
  TProfile2D *hSubProfile2D_XYZ_imag_protons[10];

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
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
  Double_t fPIDbayesPion;
  Double_t fPIDbayesKaon;
  Double_t fPIDbayesProton;
  
  TExMap *fGlobalTracksAOD; //! global tracks in AOD for FB128 **Ante**

  
  AliAnalysisTask_v02pT(const AliAnalysisTask_v02pT&);
  AliAnalysisTask_v02pT& operator=(const AliAnalysisTask_v02pT&);  
  ClassDef(AliAnalysisTask_v02pT, 1);
};

#endif
