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
  TH2F* fHistXtrkVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistYtrkVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistZtrkVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistXtpcVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistYtpcVsContrib;           //!<!  histo of vtx coord.
  TH2F* fHistZtpcVsContrib;           //!<!  histo of vtx coord.
  TH1F* fHistoNOfPileupVertSPD;       //!<! histo of SPD pileup
  TH1F* fHistoNOfSelPileupVertSPD;    //!<! histo of SPD pileup
  TH1F* fHistoNOfPileupVertMV;       //!<! histo of SPD pileup
  TH1F* fHistoNOfSelPileupVertMV;    //!<! histo of SPD pileup
  Bool_t  fUsePhysSel;                // flag use/not use phys sel
  Int_t   fTriggerMask;               // mask used in physics selection
  Double_t fMaxMult;                  // upper limit of multiplicity plots
  Double_t fSPDContributorsCut;       // cut on SPD vertex for pileup
  Double_t fSPDZDiffCut;              // cut on SPD vertex for pileup
  Int_t  fMVContributorsCut;          // cut on MV pileup vertex
  Double_t fMVCChi2Cut;               //  cut on MV pileup vertex
  Float_t  fMVWeiZDiffCut;            // cut on MV pileup vertex
  Bool_t  fMVCheckPlpFromDifferentBC; // cut on MV pileup vertex
  Bool_t  fReadMC;                    // flag read/not-read MC truth info


  ClassDef(AliAnalysisTaskCheckVertexAOD,1);
};


#endif
