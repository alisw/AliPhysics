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
#ifndef AliAnalysisTaskMeanptFluctuationPP_H
#define AliAnalysisTaskMeanptFluctuationPP_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TTree.h"
#include "AliEventCuts.h"

class TList;
class TTree;

class AliESDEvent;
class AliESDtrackCuts;
class AliMCEvent;

class TH1F;
class TH2F;
class TH3F;
class TProfile;
class AliPIDResponse;
class AliMultSelection;


class AliAnalysisTaskMeanptFluctuationPP : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMeanptFluctuationPP();
  AliAnalysisTaskMeanptFluctuationPP(const char *name);
  virtual ~AliAnalysisTaskMeanptFluctuationPP();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void           SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts);
  Double_t       MyRapidity(Double_t rE, Double_t rPz) const;
  //---------------------------------------------------------------------------------------
  Double_t      GetTOFBeta(AliVTrack *esdtrack);
  Bool_t        MatchTOF(AliVTrack *vtrack);

  //Functions for AddTask
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  void SetTreeName(TString TreeName)
  {
    fTreeName = TreeName;
  }
  void SetVzRangeMax(Double_t VzMax)
  {
    fVertexZMax = VzMax;
  }
 void SetDCAXYRangeMax(Double_t dcaxy)          
  {
    fDCAxyMax = dcaxy;
  }
  void SetDCAXYPtDep_choice(Int_t dcaxy_choice)          
  {
    fDCAxy_ptdep = dcaxy_choice;
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

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 
 private:
  enum
  {
    kMaxTrack=500000
  };

  TTree *fTreeEvent;
  TList *fListHist; 
  TList *fQAList;
  AliPIDResponse   *fPIDResponse;
  AliESDtrackCuts  *fESDtrackCuts;
  AliEventCuts fEventCuts;
  Float_t fTreeVariableCentrality;
  Float_t fvertex;
  Float_t fQ1_ptmax2;
  Float_t fQ2_ptmax2;
  Float_t fQ3_ptmax2;
  Float_t fQ4_ptmax2;
  Float_t fNch_ptmax2;
  Float_t fQ1_ptmax3;
  Float_t fQ2_ptmax3;
  Float_t fQ3_ptmax3;
  Float_t fQ4_ptmax3;
  Float_t fNch_ptmax3;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TString fTreeName;
  Double_t fVertexZMax;
  Int_t fDCAxy_ptdep;
  Double_t fDCAxyMax;
  Double_t fDCAzMax;
  Double_t fChi2TPC;
  Double_t fChi2ITS;
  Double_t fNCrossedRowsTPC;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Int_t fTreeEventNTrack;
  TProfile *Profile_mean_term1_ptmax2;
  TProfile *Profile_var_term1_ptmax2;
  TProfile *Profile_var_term2_ptmax2;
  TProfile *Profile_skewness_term1_ptmax2;
  TProfile *Profile_skewness_term2_ptmax2;
  TProfile *Profile_skewness_term3_ptmax2;
  TProfile *Profile_kurtosis_term1_ptmax2;
  TProfile *Profile_kurtosis_term2_ptmax2;
  TProfile *Profile_kurtosis_term3_ptmax2;
  TProfile *Profile_kurtosis_term4_ptmax2;
  TProfile *Profile_mean_term1_ptmax3;
  TProfile *Profile_var_term1_ptmax3;
  TProfile *Profile_var_term2_ptmax3;
  TProfile *Profile_skewness_term1_ptmax3;
  TProfile *Profile_skewness_term2_ptmax3;
  TProfile *Profile_skewness_term3_ptmax3;
  TProfile *Profile_kurtosis_term1_ptmax3;
  TProfile *Profile_kurtosis_term2_ptmax3;
  TProfile *Profile_kurtosis_term3_ptmax3;
  TProfile *Profile_kurtosis_term4_ptmax3;
  TH2D *fHistTPConlyVsV0MBeforeAliEventCut;
  TH2D *fHistTPConlyVsV0MAfterAliEventCut;
  TH2D *fHistCL0VsV0MBeforeAliEventCut;
  TH2D *fHistCL0VsV0MAfterAliEventCut;
  TH2D *fHistTPCtracksVsITStrkltsBeforeAliEventCut;
  TH2D *fHistTPCtracksVsITStrkltsAfterAliEventCut;
  TH2D *fHistTPConlyVsV0MBefore;
  TH2D *fHistTPConlyVsV0MAfter;
  TH2D *fHistCL0VsV0MBefore;
  TH2D *fHistCL0VsV0MAfter;
  TH2D *fHistTPCtracksVsITStrkltsBefore;
  TH2D *fHistTPCtracksVsITStrkltsAfter;
  TH2D *fHistTPConlyVsV0MAfterITSTPCcorelationCut;
  TH2D *fHistCL0VsV0MAfterITSTPCcorelationCut;
  TH2D *fHistTPCtracksVsITStrkltsAfterITSTPCcorelationCut;
  TH2D *fHistTPCrefitVsITSrefitBeforeAliEventCut;
  TH2D *fHistTPCrefitVsITSrefitAfterAliEventCut;
  TH2D *fHistTPCrefitVsITSrefitBefore;
  TH2D *fHistTPCrefitVsITSrefitAfter;
  TH2D *fHistTPCrefitVsITSrefitAfterITSTPCcorelationCut;
  TF1 *fCenCutLowPU;
  TF1 *fCenCutHighPU;
  TF1 *fSPDCutPU;
  TF1 *fV0CutPU;
  TF1 *fMultCutPU;
  TH1D *Hist_fTreeTrackVariableDcaXY;
  TH1D *Hist_fTreeTrackVariableDcaZ;
  TH1D *Hist_fTreeTrackVariableTpcNCls;
  TH1D *Hist_fTreeTrackVariableTpcNCrossedRows;
  TH1D *Hist_fTreeTrackVariableLeastRatioCrossedRowsOverFindable;
  TH1D *Hist_fTreeTrackVariableChiSqrPerTpcCls;
  TH1D *Hist_fTreeTrackVariableChiSqrPerItsCls;
  TH1D *Hist_Eta;
  TH1D *Hist_Pt;

  Float_t fTreeTrackVariableDcaXY[kMaxTrack]; //!
  Float_t fTreeTrackVariableDcaZ[kMaxTrack]; //!
  Float_t fTreeTrackVariableTpcNCls[kMaxTrack]; //!
  Float_t fTreeTrackVariableTpcNCrossedRows[kMaxTrack]; //!
  Float_t fTreeTrackVariableLeastRatioCrossedRowsOverFindable[kMaxTrack]; //!
  Float_t fTreeTrackVariableChiSqrPerTpcCls[kMaxTrack]; //!
  Float_t fTreeTrackVariableChiSqrPerItsCls[kMaxTrack]; //!

  
  

  //Custom Functions:
  Bool_t AcceptEventAfterPileUpCut(AliESDEvent* fESD);
  
  AliAnalysisTaskMeanptFluctuationPP(const AliAnalysisTaskMeanptFluctuationPP&);
  AliAnalysisTaskMeanptFluctuationPP& operator=(const AliAnalysisTaskMeanptFluctuationPP&);  
  ClassDef(AliAnalysisTaskMeanptFluctuationPP, 1);
};

#endif
