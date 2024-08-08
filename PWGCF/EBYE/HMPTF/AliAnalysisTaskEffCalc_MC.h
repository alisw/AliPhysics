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
#ifndef AliAnalysisTaskEffCalc_MC_H
#define AliAnalysisTaskEffCalc_MC_H

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
class AliPIDCombined;
class AliMultSelection;


class AliAnalysisTaskEffCalc_MC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEffCalc_MC();
  AliAnalysisTaskEffCalc_MC(const char *name);
  virtual ~AliAnalysisTaskEffCalc_MC();
  
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
  void SetEtaCut(Double_t etamax)
  {
    this->fEtaMax = etamax;
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
  AliMCEvent *fMCevent;
  AliStack *fMCstack;
  AliVEvent *fInputEvent;
  AliPIDResponse   *fPIDResponse;
  AliPIDCombined   *fPIDCombined;
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
  TTree *fTreeEvent_gen;
  TTree *fTreeEvent_rec;
  Float_t fTreeVariableCentrality;
  

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
  TH1F *hist_GenHadronPlus;
  TH1F *hist_GenHadronMinus;
  TH1F *hist_RecHadronPlus;
  TH1F *hist_RecHadronMinus;

  //2D histograms for calculating efficiency
  TH2D *hist2D_GenProtonPlus;
  TH2D *hist2D_GenProtonMinus;
  TH2D *hist2D_RecProtonPlus;
  TH2D *hist2D_RecProtonMinus;



  //Event and track parameters to be coordinated 
  Double_t fVertexZMax;
  Int_t fFBNo;
  Double_t fChi2TPC;
  Double_t fChi2ITS;
  Double_t fPIDnSigmaPionCut;
  Double_t fPIDnSigmaKaonCut;
  Double_t fPIDnSigmaProtonCut;
  Double_t fTPCcrossedrows;
  Double_t fEtaMax;
  
  
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

  //Bayesian PID related variables
  Int_t fBayesianPID_flag;
  Double_t fPIDbayesPion;
  Double_t fPIDbayesKaon;
  Double_t fPIDbayesProton;
  
 
  //Efficienc hist per centralities
  TH1D *fEffPionPlus[9];
  TH1D *fEffKaonPlus[9];
  TH1D *fEffProtonPlus[9];
  TH1D *fEffPionMinus[9];
  TH1D *fEffKaonMinus[9];
  TH1D *fEffProtonMinus[9];
  

  
  AliAnalysisTaskEffCalc_MC(const AliAnalysisTaskEffCalc_MC&);
  AliAnalysisTaskEffCalc_MC& operator=(const AliAnalysisTaskEffCalc_MC&);  
  ClassDef(AliAnalysisTaskEffCalc_MC, 1);
};

#endif
