#ifndef ALIANALYSISTASKCHECKVERTEXAOD
#define ALIANALYSISTASKCHECKVERTEXAOD

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskCheckVertexAOD
// AliAnalysisTaskSE to extract QA plots for vertices in AODs
// 
//
// Author:
//          F. Prino, prino@to.infn.it
//          
//*************************************************************************

class TList;
class TH1F;
class TH2F;
class TH3F;
class TString;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckVertexAOD : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskCheckVertexAOD();
  virtual ~AliAnalysisTaskCheckVertexAOD();

  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);

  void SetReadMC(Bool_t optMC=kTRUE){
    fReadMC=optMC;
  }
  void SelectedGeneratorName(TString name){
    fGenerToKeep=name.Data(); fSelectOnGenerator=kTRUE;}
  void ExcludedGeneratorName(TString name){
    fGenerToExclude=name.Data(); fSelectOnGenerator=kTRUE;}
  void SetUsePhysicsSelection(Bool_t opt=kTRUE){
    fUsePhysSel=opt;
  }
  void SetTriggerMask(Int_t mask){
    fTriggerMask=mask;
  }
  void SetUpperMultiplicity(Double_t maxMult){
    fMaxMult=maxMult;
  }
  void SetCutOnContribToSPDPileupVert(Int_t cutc) {
    fSPDContributorsCut=cutc;
  }
  void SetCutOnSPDZDiff(Double_t cutz) {
    fSPDZDiffCut=cutz;
  }
  void ApplyPbPbOutOfBunchPileupCut(Bool_t opt){
    fApplyPbPbOutOfBunchPileupCut=opt;
  }

 private:

  AliAnalysisTaskCheckVertexAOD(const AliAnalysisTaskCheckVertexAOD &source);
  AliAnalysisTaskCheckVertexAOD& operator=(const AliAnalysisTaskCheckVertexAOD &source);
  
  TList*  fOutput;                    //!<!  list of output histos
  TH1F* fHistNEvents;                 //!<!  histo with N of events  
  TH1F* fHistAllVtxType;              //!<!  histo of vtx type
  TH1F* fHistPrimVtxType;             //!<!  histo of vtx type
  TH2F* fHistXspdVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistYspdVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistZspdVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistZspdOnlyZVsContrib;      //!<!  histo of vtx coord.
  TH2F* fHistXtrkVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistYtrkVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistZtrkVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistXtpcVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistYtpcVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistZtpcVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistXspdVsMult;              //!<!  histo of vtx coord.
  TH2F* fHistYspdVsMult;              //!<!  histo of vtx coord.
  TH2F* fHistZspdVsMult;              //!<!  histo of vtx coord.
  TH2F* fHistZspdOnlyZVsMult;         //!<!  histo of vtx coord.
  TH2F* fHistXtrkVsMult;              //!<!  histo of vtx coord.
  TH2F* fHistYtrkVsMult;              //!<!  histo of vtx coord.
  TH2F* fHistZtrkVsMult;              //!<!  histo of vtx coord.
  TH2F* fHistXtpcVsMult;              //!<!  histo of vtx coord.
  TH2F* fHistYtpcVsMult;              //!<!  histo of vtx coord.
  TH2F* fHistZtpcVsMult;              //!<!  histo of vtx coord.
  TH2F* fHistPrimVtxTypeVsCent;       //!<!  histo of vtx type
  TH2F* fHistXspdVsCent;              //!<!  histo of vtx coord.
  TH2F* fHistYspdVsCent;              //!<!  histo of vtx coord.
  TH2F* fHistZspdVsCent;              //!<!  histo of vtx coord.
  TH2F* fHistZspdOnlyZVsCent;         //!<!  histo of vtx coord.
  TH2F* fHistContrSpdVsCent;          //!<!  histo of vtx contrib
  TH2F* fHistContrSpdOnlyZVsCent;     //!<!  histo of vtx contrib
  TH2F* fHistXtrkVsCent;              //!<!  histo of vtx coord.
  TH2F* fHistYtrkVsCent;              //!<!  histo of vtx coord.
  TH2F* fHistZtrkVsCent;              //!<!  histo of vtx coord.
  TH2F* fHistContrTrkVsCent;          //!<!  histo of vtx contrib
  TH2F* fHistXtpcVsCent;              //!<!  histo of vtx coord.
  TH2F* fHistYtpcVsCent;              //!<!  histo of vtx coord.
  TH2F* fHistZtpcVsCent;              //!<!  histo of vtx coord.
  TH2F* fHistContrTpcVsCent;          //!<!  histo of vtx contrib
  TH2F* fHistXtrkResidVsMult;         //!<!  histos for resolution (MC)
  TH2F* fHistYtrkResidVsMult;         //!<!  histos for resolution (MC)
  TH2F* fHistZtrkResidVsMult;         //!<!  histos for resolution (MC)
  TH2F* fHistXtrkResidVsCent;         //!<!  histos for resolution (MC)
  TH2F* fHistYtrkResidVsCent;         //!<!  histos for resolution (MC)
  TH2F* fHistZtrkResidVsCent;         //!<!  histos for resolution (MC)
  TH2F* fHistNtracklVsZtrue;          //!<!  histo of vtx coord.
  TH1F* fHistoNOfPileupVertSPD;       //!<! histo of SPD pileup
  TH1F* fHistoNOfSelPileupVertSPD;    //!<! histo of SPD pileup
  TH1F* fHistoNOfPileupVertMV;        //!<! histo of SPD pileup
  TH1F* fHistoNOfSelPileupVertMV;     //!<! histo of SPD pileup
  TH2F* fHistoV0MultVsNclsTPC;        //!<! control histo for TPC pileup
  Bool_t fUsePhysSel;                 /// flag use/not use phys sel
  Int_t  fTriggerMask;                /// mask used in physics selection
  Bool_t fSelectOnGenerator;          /// flag to select events with generator name
  TString fGenerToKeep;               /// generator name to analyse
  TString fGenerToExclude;            /// generator name to exclude
  Double_t fMaxMult;                  /// upper limit of multiplicity plots
  Double_t fSPDContributorsCut;       /// cut on SPD vertex for pileup
  Double_t fSPDZDiffCut;              /// cut on SPD vertex for pileup
  Int_t  fMVContributorsCut;          /// cut on MV pileup vertex
  Double_t fMVCChi2Cut;               ///  cut on MV pileup vertex
  Float_t  fMVWeiZDiffCut;            /// cut on MV pileup vertex
  Bool_t  fMVCheckPlpFromDifferentBC; /// cut on MV pileup vertex
  Bool_t  fApplyPbPbOutOfBunchPileupCut; /// switch for cut on TPC clus vs. V0
  Bool_t  fReadMC;                    /// flag read/not-read MC truth info


  ClassDef(AliAnalysisTaskCheckVertexAOD,11);
};


#endif
