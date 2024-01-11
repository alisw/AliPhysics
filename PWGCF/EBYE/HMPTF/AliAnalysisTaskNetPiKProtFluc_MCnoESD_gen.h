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
#ifndef AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen_H
#define AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen_H

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


class AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen();
  AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen(const char *name);
  virtual ~AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen();
  
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
  
  void SetNoResonanceChoice(int choice)
  {
    fNoResoChoice = choice;
  }
  
 private:
  
  TTree *fTreeEvent;
  TList *fListHist; 
  AliPIDResponse   *fPIDResponse;
  AliESDtrackCuts  *fESDtrackCuts;
  AliEventCuts fEventCuts;
  
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Variables used for this code ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Float_t fTreeVariableCentrality; 
  //generated
  Float_t fNoGenKaonPlus_ptmax2;
  Float_t fNoGenKaonMinus_ptmax2;
  Float_t fNoGenKaonPlus_ptmax3;
  Float_t fNoGenKaonMinus_ptmax3;
  Float_t fNoGenPionPlus_ptmax2;
  Float_t fNoGenPionMinus_ptmax2;
  Float_t fNoGenPionPlus_ptmax3;
  Float_t fNoGenPionMinus_ptmax3;
  Float_t fNoGenProtonPlus_ptmax2;
  Float_t fNoGenProtonMinus_ptmax2;
  Float_t fNoGenProtonPlus_ptmax3;
  Float_t fNoGenProtonMinus_ptmax3;
  Float_t fvertex;
  Float_t fNch_eta0pt5;
  Float_t fNpart_1;
  Float_t fNpart_2;
  
  TH2F *hist2D_pt_gen_centrality;
  TH2F *hist2D_pt_rec_centrality;
  Int_t fMCchoice;
  Int_t fNoResoChoice;
  		
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //Float_t fvertex;
  //Float_t fQ1;
  //Float_t fQ2;
  //Float_t fQ3;
  //Float_t fQ4;
  //Float_t fNch;
  
  TH1D *hist_centrality_beforecut;
  
  AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen(const AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen&);
  AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen& operator=(const AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen&);  
  ClassDef(AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen, 1);
};

#endif
