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
#ifndef AliAnalysisTaskHptcorrId_tree_H
#define AliAnalysisTaskHptcorrId_tree_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TTree.h"
#include "AliEventCuts.h"

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
class AliPIDResponse;
class AliMultSelection;


class AliAnalysisTaskHptcorrId_tree : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHptcorrId_tree();
  AliAnalysisTaskHptcorrId_tree(const char *name);
  virtual ~AliAnalysisTaskHptcorrId_tree();
  
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
  void SetEtaCut(Double_t eta_max)
  {
    this->fEtaMax = eta_max;
  
  }
  void SetTreeName(TString TreeName)
  {
    fTreeName = TreeName;
  }
  void SetPtMaxMin(float ptmax, float ptmin)
  {
    fPtMax=ptmax;
    fPtMin=ptmin;
  }
  /*
  void SetDCAXYRangeMax(Double_t dcaxy)          
  {
    fDCAxyMax = dcaxy;
  }
  void SetDCAZRangeMax(Double_t dcaz)
  {
    fDCAzMax = dcaz;
  }
  void SetMinNCrossedRowsTPCRange(Double_t nCrossedRow)
  {
    fNCrossedRowsTPC = nCrossedRow;
  }
  

  void SetEtaCut(Double_t EtaMax)
  {
    fEtaMax=EtaMax;
  }
  */

  
 private:
  
  
  AliESDEvent *fESDevent;
  AliAODEvent *fAODevent;
  AliVEvent *fInputEvent;
  AliPIDResponse   *fPIDResponse;
  AliESDtrackCuts  *fESDtrackCuts;
  AliESDtrackCuts  *fESDtrackCuts_primary;
  AliAnalysisUtils *fUtils;
  AliEventCuts fAODeventCuts;
  TList *fOutputList;
  TList *fQAList;
  Double_t fMultLow;
  Double_t fMultHigh;
  TH1D *hNumberOfEvents;
  
  //ftreeEvent object variables
  TTree *fTreeEvent;
  Float_t fTreeVariableCentrality;
  Float_t fQ1_prot;
  Float_t fQ2_prot;
  Float_t fQ3_prot;
  Float_t fQ4_prot;
  Float_t fNch_prot;
  Float_t fQ1_pi;
  Float_t fQ2_pi;
  Float_t fQ3_pi;
  Float_t fQ4_pi;
  Float_t fNch_pi;
  Float_t fQ1_K;
  Float_t fQ2_K;
  Float_t fQ3_K;
  Float_t fQ4_K;
  Float_t fNch_K;
 
  Float_t fPtMax;
  Float_t fPtMin;
  
  //Efficiency list of histograms
  TList *fListTRKCorr; 
  TH1D *fHistMCEffKaonPlus;
  TH1D *fHistMCEffKaonMinus;
  TH1D *fHistMCEffPionPlus;
  TH1D *fHistMCEffPionMinus;
  TH1D *fHistMCEffProtonPlus;
  TH1D *fHistMCEffProtonMinus;
  TH1D *fEffPionPlus[9];
  TH1D *fEffKaonPlus[9];
  TH1D *fEffProtonPlus[9];
  TH1D *fEffPionMinus[9];
  TH1D *fEffKaonMinus[9];
  TH1D *fEffProtonMinus[9];
 
  Double_t fVertexZMax;
  Int_t fFBNo;
  Double_t fChi2TPC;
  Double_t fChi2ITS;
  Double_t fPIDnSigmaPionCut;
  Double_t fPIDnSigmaKaonCut;
  Double_t fPIDnSigmaProtonCut;
  Double_t fTPCcrossedrows;
  Double_t fEtaMax;
  TString fTreeName;

  //Pileup cut val
  Int_t fPileupCutVal;
  
  //Flag to select which centrality estimator
  Int_t fCentralityEstimator_flag;
  

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
 
  
  AliAnalysisTaskHptcorrId_tree(const AliAnalysisTaskHptcorrId_tree&);
  AliAnalysisTaskHptcorrId_tree& operator=(const AliAnalysisTaskHptcorrId_tree&);  
  ClassDef(AliAnalysisTaskHptcorrId_tree, 1);
};

#endif
