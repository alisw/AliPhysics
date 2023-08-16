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
#ifndef AliAnalysisTaskDiffPtFStrange_H
#define AliAnalysisTaskDiffPtFStrange_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TTree.h"
#include "AliEventCuts.h"

class TList;
class TTree;

class AliESDEvent;
//class AliAODEvent;
class AliESDtrackCuts;
class AliAODv0;
class AliAODcascade;
class AliMCEvent;
class AliVTrack;
class AliAODTrack;

class TH1D;
class TH2D;
class TH3D;
class TProfile;
class AliPIDResponse;
class AliMultSelection;


class AliAnalysisTaskDiffPtFStrange : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDiffPtFStrange();
  AliAnalysisTaskDiffPtFStrange(const char *name);
  virtual ~AliAnalysisTaskDiffPtFStrange();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t GetEvent();
  Bool_t PassedTrackQualityCuts (AliAODTrack *track);
  Bool_t PassedDaughterTrackDCAtoVertexSelectionCutsV0 (AliAODv0 *v0);
  Bool_t PassedDaughterTrackDCAtoVertexSelectionCutsCascade(AliAODcascade *cascade);
  Bool_t PassedV0SelectionTopologicalCuts (AliAODv0 *v0);
  Bool_t PassedCascadeSelectionCuts (AliAODcascade *cascade);
  Bool_t AcceptEventAfterPileUpCut(AliESDEvent* fESD);
  Bool_t IsLambdaCandidate (AliAODv0 *V0, AliAODTrack *pos, AliAODTrack *neg, Double_t PIDcut);
  Bool_t IsAntiLambdaCandidate (AliAODv0 *V0, AliAODTrack *pos, AliAODTrack *neg, Double_t PIDcut);
  Bool_t IsXiCandidate (AliAODcascade *casc, AliAODTrack *pos, AliAODTrack *neg, AliAODTrack *bac, Double_t PIDcut);
  Bool_t IsAntiXiCandidate (AliAODcascade *casc, AliAODTrack *pos, AliAODTrack *neg, AliAODTrack *bac, Double_t PIDcut);
  Bool_t KaonSelector (AliVTrack *track); 
  Bool_t PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type, Double_t PIDcut);
  Bool_t PassedSingleParticlePileUpCuts(AliAODTrack *track);

  /*
  void SetVzRangeMax(Double_t VzMax)
  {
    this->fVertexZMax = VzMax;
  }
  void SetDCAXYRangeMax(Double_t dcaxy)          
  {
    fDCAxyMax = dcaxy;
  }
  void SetDCAZRangeMax(Double_t dcaz)
  {
    fDCAzMax = dcaz;
  }
  void SetMaxChi2PerTPCClusterRange(Double_t chi2tpc)
  {
    fChi2TPC = chi2tpc;
  }
  void SetMaxChi2PerITSClusterRange(Double_t chi2its)
  {
    fChi2ITS = chi2its;
  }
  void SetMinNCrossedRowsTPCRange(Double_t nCrossedRow)
  {
    fNCrossedRowsTPC = nCrossedRow;
  }
  void SetTreeName(TString TreeName)
  {
    fTreeName = TreeName;
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
  TH1D *hNumberOfV0;
  TH2D *histMassLambda_vs_Pt_beforeMasscut;
  TH2D *histMassAntiLambda_vs_Pt_beforeMasscut;
  TH2D *histMassLambda_vs_Pt;
  TH2D *histMassAntiLambda_vs_Pt;
  TH1D *hNumberOfLambda;
  TH1D *hNumberOfAntiLambda;
  TH1D *histLambdaPt;
  TH1D *histAntiLambdaPt;
  TH1D *hNumberOfCascades;
  TH1D *hNumberOfXi;
  TH1D *hNumberOfAntiXi;
  TH1D *histXiPt;
  TH1D *histAntiXiPt;
  TH1D *hNumberOfKaonPlus;
  TH1D *hNumberOfKaonMinus;
  TH2D *histMassXi_vs_Pt_beforeMasscut;
  TH2D *histMassAntiXi_vs_Pt_beforeMasscut;
  TH2D *histMassXi_vs_Pt;
  TH2D *histMassAntiXi_vs_Pt;
  TTree *fTreeEvent;
  Float_t fTreeVariableCentrality;
  Float_t fMeanPt_less0;
  Float_t fMeanPt_greaterEtaMin;
  Float_t fPt_factor_Lambda[14];
  Float_t fPt_factor_Xi[14];
  Float_t fNoXi_ptmax2;
  Float_t fNoAntiXi_ptmax2;
  Float_t fNoXi_ptmax3;
  Float_t fNoAntiXi_ptmax3;
  Float_t fNoKaonPlus_ptmax2;
  Float_t fNoKaonMinus_ptmax2;
  Float_t fNoKaonPlus_ptmax3;
  Float_t fNoKaonMinus_ptmax3;
  

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
 
  
  AliAnalysisTaskDiffPtFStrange(const AliAnalysisTaskDiffPtFStrange&);
  AliAnalysisTaskDiffPtFStrange& operator=(const AliAnalysisTaskDiffPtFStrange&);  
  ClassDef(AliAnalysisTaskDiffPtFStrange, 1);
};

#endif
