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
#ifndef AliAnalysisTaskDiffPtFluc_PIDhadrons_H
#define AliAnalysisTaskDiffPtFluc_PIDhadrons_H

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


class AliAnalysisTaskDiffPtFluc_PIDhadrons : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDiffPtFluc_PIDhadrons();
  AliAnalysisTaskDiffPtFluc_PIDhadrons(const char *name);
  virtual ~AliAnalysisTaskDiffPtFluc_PIDhadrons();
  
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
  TH1D *hNumberOfKaonEtaLess0;
  TH1D *hNumberOfPionEtaLess0;
  TH1D *hNumberOfProtonEtaLess0;

  //ftreeEvent object variables
  TTree *fTreeEvent;
  Float_t fTreeVariableCentrality;
  Float_t fMeanPt_less0;
  Float_t fMeanPt_greaterEtaMin;
  Float_t fPt_factor[14];
  Float_t fPt_factor_pion[14];
  Float_t fPt_factor_kaon[14];
  Float_t fPt_factor_proton[14];
  
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
  Double_t fPIDnSigmaCut;
  Double_t fTPCcrossedrows;
  Double_t fEtaMin;
  TString fTreeName;
  

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
 
  
  AliAnalysisTaskDiffPtFluc_PIDhadrons(const AliAnalysisTaskDiffPtFluc_PIDhadrons&);
  AliAnalysisTaskDiffPtFluc_PIDhadrons& operator=(const AliAnalysisTaskDiffPtFluc_PIDhadrons&);  
  ClassDef(AliAnalysisTaskDiffPtFluc_PIDhadrons, 1);
};

#endif
