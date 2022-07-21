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
#ifndef AliAnalysisTaskResonanceVsMultiplicityROOT6_H
#define AliAnalysisTaskResonanceVsMultiplicityROOT6_H

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


class AliAnalysisTaskResonanceVsMultiplicityROOT6 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskResonanceVsMultiplicityROOT6();
  AliAnalysisTaskResonanceVsMultiplicityROOT6(const char *name);
  virtual ~AliAnalysisTaskResonanceVsMultiplicityROOT6();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void           SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts);
  Double_t       MyRapidity(Double_t rE, Double_t rPz) const;
  //---------------------------------------------------------------------------------------
  Double_t      GetTOFBeta(AliVTrack *esdtrack);
  Bool_t        MatchTOF(AliVTrack *vtrack);


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

  
 private:
  
  TTree *fTreeEvent;
  TList *fListHist;
  TList *fQAList;
  AliPIDResponse   *fPIDResponse;
  AliESDtrackCuts  *fESDtrackCuts;
  AliEventCuts fEventCuts;
  //Float_t fTreeVariableCentrality;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Variables used for this code ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Float_t fTreeVariableCentrality; 
  Float_t fvertex;
  Float_t fNch_eta0pt5;
  Float_t fQ1_ptmax3;
  Float_t fQ2_ptmax3;
  Float_t fQ3_ptmax3;
  Float_t fQ4_ptmax3;
  Float_t fNch_ptmax3;
  Float_t fQ1_ptmax2;
  Float_t fQ2_ptmax2;
  Float_t fQ3_ptmax2;
  Float_t fQ4_ptmax2;
  Float_t fNch_ptmax2;
  Float_t fQ1_gen[4];
  Float_t fQ2_gen[4];
  Float_t fQ3_gen[4];
  Float_t fQ4_gen[4];
  Float_t fNch_gen[4];
  Float_t fQ1_rec[4];
  Float_t fQ2_rec[4];
  Float_t fQ3_rec[4];
  Float_t fQ4_rec[4];
  Float_t fNch_rec[4];
  		
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TH2F *hist2D_pt_gen_centrality;
  TH2F *hist2D_pt_rec_centrality;
  TProfile *Profile_mean_term1;
  TProfile *Profile_var_term1;
  TProfile *Profile_var_term2;
  TProfile *Profile_skewness_term1;
  TProfile *Profile_skewness_term2;
  TProfile *Profile_skewness_term3;
  TProfile *Profile_kurtosis_term1;
  TProfile *Profile_kurtosis_term2;
  TProfile *Profile_kurtosis_term3;
  TProfile *Profile_kurtosis_term4;
  TH1D *hist_centrality_beforecut;
  TH1D *Hist_fTreeTrackVariableDcaXY;
  TH1D *Hist_fTreeTrackVariableDcaZ;
  TH1D *Hist_fTreeTrackVariableTpcNCls;
  TH1D *Hist_fTreeTrackVariableTpcNCrossedRows;
  TH1D *Hist_fTreeTrackVariableLeastRatioCrossedRowsOverFindable;
  TH1D *Hist_fTreeTrackVariableChiSqrPerTpcCls;
  TH1D *Hist_fTreeTrackVariableChiSqrPerItsCls;
  TH1D *Hist_Eta;
  TH1D *Hist_Pt;

  //Pileup Histograms
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
  

  //Custom Functions:
  TF1 *fCenCutLowPU;
  TF1 *fCenCutHighPU;
  TF1 *fSPDCutPU;
  TF1 *fV0CutPU;
  TF1 *fMultCutPU;
  Bool_t AcceptEventAfterPileUpCut(AliESDEvent* fESD);


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

 
  
  AliAnalysisTaskResonanceVsMultiplicityROOT6(const AliAnalysisTaskResonanceVsMultiplicityROOT6&);
  AliAnalysisTaskResonanceVsMultiplicityROOT6& operator=(const AliAnalysisTaskResonanceVsMultiplicityROOT6&);  
  ClassDef(AliAnalysisTaskResonanceVsMultiplicityROOT6, 1);
};

#endif
