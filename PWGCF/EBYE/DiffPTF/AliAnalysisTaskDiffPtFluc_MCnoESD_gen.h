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
#ifndef AliAnalysisTaskDiffPtFluc_MCnoESD_gen_H
#define AliAnalysisTaskDiffPtFluc_MCnoESD_gen_H

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


class AliAnalysisTaskDiffPtFluc_MCnoESD_gen : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDiffPtFluc_MCnoESD_gen();
  AliAnalysisTaskDiffPtFluc_MCnoESD_gen(const char *name);
  virtual ~AliAnalysisTaskDiffPtFluc_MCnoESD_gen();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void           SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts);
  Double_t       MyRapidity(Double_t rE, Double_t rPz) const;
  //---------------------------------------------------------------------------------------
  Double_t      GetTOFBeta(AliVTrack *esdtrack);
  Bool_t        MatchTOF(AliVTrack *vtrack);
  
  void SetMCGeneratorChoice(int MC_choice)
  {
    fMCchoice = MC_choice;
  }
  
 private:
  
  TTree *fTreeEvent;
  TList *fListHist; 
  AliPIDResponse   *fPIDResponse;
  AliESDtrackCuts  *fESDtrackCuts;
  AliEventCuts fEventCuts;
  
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Variables used for this code ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Float_t fTreeVariableCentrality; 
  Float_t fvertex;
  Float_t fNch_eta0pt5;
  Float_t fNpart_1;
  Float_t fNpart_2;
  Float_t fQ1_gen[2];
  Float_t fQ2_gen[2];
  Float_t fQ3_gen[2];
  Float_t fQ4_gen[2];
  Float_t fNch_gen[2];
  Float_t fQ1_rec[2];
  Float_t fQ2_rec[2];
  Float_t fQ3_rec[2];
  Float_t fQ4_rec[2];
  Float_t fNch_rec[2];
  TH2F *hist2D_pt_gen_centrality;
  TH2F *hist2D_pt_rec_centrality;
  Int_t fMCchoice;
  Float_t fMeanPt_less0;
  Float_t fMeanPt_greaterEtaMin;
  Float_t fPt_factor[14];
  Float_t fPt_factor_pion[14];
  Float_t fPt_factor_kaon[14];
  Float_t fPt_factor_proton[14];		
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //Float_t fvertex;
  //Float_t fQ1;
  //Float_t fQ2;
  //Float_t fQ3;
  //Float_t fQ4;
  //Float_t fNch;
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
  
  AliAnalysisTaskDiffPtFluc_MCnoESD_gen(const AliAnalysisTaskDiffPtFluc_MCnoESD_gen&);
  AliAnalysisTaskDiffPtFluc_MCnoESD_gen& operator=(const AliAnalysisTaskDiffPtFluc_MCnoESD_gen&);  
  ClassDef(AliAnalysisTaskDiffPtFluc_MCnoESD_gen, 1);
};

#endif
