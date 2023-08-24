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
#ifndef AliAnalysisTaskResonanceVsMultiplicityMC_H
#define AliAnalysisTaskResonanceVsMultiplicityMC_H

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


class AliAnalysisTaskResonanceVsMultiplicityMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskResonanceVsMultiplicityMC();
  AliAnalysisTaskResonanceVsMultiplicityMC(const char *name);
  virtual ~AliAnalysisTaskResonanceVsMultiplicityMC();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void           SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts);
  Double_t       MyRapidity(Double_t rE, Double_t rPz) const;
  //---------------------------------------------------------------------------------------
  Double_t      GetTOFBeta(AliVTrack *esdtrack);
  Bool_t        MatchTOF(AliVTrack *vtrack);

//Functions required for AddTask
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  void SetVzRangeMax(Double_t VzMax)
  {
    fVertexZ_max = VzMax;
  }
  void SetDCAXYRangeMax(Double_t dcaxy)          
  {
    fDCAxy_max = dcaxy;
  }
  void SetDCAZRangeMax(Double_t dcaz)
  {
    fDCAz_max = dcaz;
  }
  void SetMaxChi2PerTPCClusterRange(Double_t chi2tpc)
  {
    fchi2TPC_max = chi2tpc;
  }
  void SetMaxChi2PerITSClusterRange(Double_t chi2its)
  {
    fchi2ITS_max = chi2its;
  }
  void SetMinNCrossedRowsTPCRange(Double_t nCrossedRow)
  {
    fNCR_min = nCrossedRow;
  }
 void SetTreeName(TString TreeName)
 {
    fTreeName = TreeName;
 }
 void SetEtaCut(Double_t EtaMax)
{
    fEta_max=EtaMax;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
 private:
  
  TTree *fTreeEvent;
  TList *fListHist; 
  AliPIDResponse   *fPIDResponse;
  AliESDtrackCuts  *fESDtrackCuts;
  AliEventCuts fEventCuts;
  //Float_t fTreeVariableCentrality;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
//Variables used for this code 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Float_t fTreeVariableCentrality; 
  Float_t fvertex;
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
  TH2F *hist2D_pt_gen_centrality;
  TH2F *hist2D_pt_rec_centrality;

  Double_t fEta_max;
  Double_t fDCAxy_max;
  Double_t fDCAz_max;
  Double_t fchi2TPC_max;
  Double_t fchi2ITS_max;
  Double_t fNCR_min;
  Double_t fVertexZ_max;
  TString fTreeName;			
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //Float_t fvertex;
  Float_t fQ1;
  Float_t fQ2;
  Float_t fQ3;
  Float_t fQ4;
  Float_t fNch;
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
  
  AliAnalysisTaskResonanceVsMultiplicityMC(const AliAnalysisTaskResonanceVsMultiplicityMC&);
  AliAnalysisTaskResonanceVsMultiplicityMC& operator=(const AliAnalysisTaskResonanceVsMultiplicityMC&);  
  ClassDef(AliAnalysisTaskResonanceVsMultiplicityMC, 1);
};

#endif
