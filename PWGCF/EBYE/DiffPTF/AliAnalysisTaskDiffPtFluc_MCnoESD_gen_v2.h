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
#ifndef AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2_H
#define AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2_H

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


class AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2();
  AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2(const char *name);
  virtual ~AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void           SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts);
  Double_t       MyRapidity(Double_t rE, Double_t rPz) const;
  //---------------------------------------------------------------------------------------
  Double_t      GetTOFBeta(AliVTrack *esdtrack);
  Bool_t        MatchTOF(AliVTrack *vtrack);
  
  void SetMCGeneratorChoice(Int_t MC_choice)
  {
    fMCchoice = MC_choice;
  }
  void SetEtaLeftCut(Double_t eta_leftval)
  {
    fEtaLeftCut = eta_leftval;
  }
  void SetEtaRightCut(Double_t eta_rightval)
  {
    fEtaRightCut = eta_rightval;
  }
  
 private:
  
  TTree *fTreeEvent;
  TList *fListHist; 
  AliPIDResponse   *fPIDResponse;
  AliESDtrackCuts  *fESDtrackCuts;
  AliEventCuts fEventCuts;
  
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Variables used for this code ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Float_t fTreeVariableCentrality; 
  Float_t fNch_eta0pt5;
  Float_t fNpart_1;
  Float_t fNpart_2;
  Double_t fEtaLeftCut;
  Double_t fEtaRightCut;
  TH2F *hist2D_pt_gen_centrality;
  TH2F *hist2D_pt_rec_centrality;
  Int_t fMCchoice;
  Float_t fPtsum_hadrons_less0;
  Float_t fPtsum_hadrons_greaterEtaMin;
  Float_t fNsum_hadrons_less0;
  Float_t fNsum_hadrons_greaterEtaMin;
  Float_t fNsum_pions_less0;
  Float_t fNsum_kaons_less0;
  Float_t fNsum_protons_less0;
  Float_t fPt_no_hadron[14];
  Float_t fPt_no_pion[14];
  Float_t fPt_no_kaon[14];
  Float_t fPt_no_proton[14];
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  
  TH1D *hist_centrality_beforecut;
  
  AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2(const AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2&);
  AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2& operator=(const AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2&);  
  ClassDef(AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2, 1);
};

#endif
