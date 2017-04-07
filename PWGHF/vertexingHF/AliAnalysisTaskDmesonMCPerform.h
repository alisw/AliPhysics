#ifndef ALIANALYSISTASKDMESONMCPERFORM_H
#define ALIANALYSISTASKDMESONMCPERFORM_H

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
/// \class Class AliAnalysisTaskDmesonMCPerform
/// \brief AliAnalysisTaskSE for performance studies of 
///        D meson hadronic decays in MC simulations
/// \author F. Prino, prino@to.infn.it
//*************************************************************************

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "TH3F.h"
#include "TH1F.h"

class AliAnalysisTaskDmesonMCPerform : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskDmesonMCPerform();
  virtual ~AliAnalysisTaskDmesonMCPerform();
  
  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
  void SetDplusAnalysisCuts(AliRDHFCutsDplustoKpipi* cts){
    fRDHFCutsDplus=cts;
  }
  void SetUseCentrality(Int_t flag){
    fRDHFCuts->SetUseCentrality(flag);
  }
  virtual void UserCreateOutputObjects();
  virtual void Init(){};
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 private:

  void FillGenLevelHistos(AliAODEvent* aod, TClonesArray *arrayMC, AliAODMCHeader *mcHeader);
  void FillCandLevelHistos(Int_t idCase, AliAODEvent* aod, TClonesArray *arrayDcand, TClonesArray *arrayMC);

  Bool_t CheckAcceptance(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau);
  AliAODRecoDecayHF* GetRecoDecay(AliAODEvent* aod, Int_t nProng, Int_t *labDau);
  void MapTrackLabels(AliAODEvent* aod);

  AliAnalysisTaskDmesonMCPerform(const AliAnalysisTaskDmesonMCPerform &source);
  AliAnalysisTaskDmesonMCPerform& operator=(const AliAnalysisTaskDmesonMCPerform &source);
  enum {kDecays=5, kMaxLabel=1000000};

  TList* fOutput;                   //!<! list send on output slot 0
  TH1F* fHistNEvents;               //!<! hist. for N. of events
  TH1F* fHistNGenD;                 //!<! hist. for N. of D's
  TH2F* fHistNCand;                 //!<! hist. for N. of D candidates
  TH3F* fHistPtYMultGen[2*kDecays];         //!<! hist. for generated D
  TH3F* fHistPtYMultGenDauInAcc[2*kDecays]; //!<! hist. for generated D in acc
  TH3F* fHistPtYMultReco[2*kDecays];        //!<! hist. for D with reco daught.
  TH3F* fHistPtYMultRecoFilt[2*kDecays];    //!<! hist. for D candidates
  TH3F* fHistPtYMultRecoSel[2*kDecays];     //!<! hist. for D candidates (sel)

  TH2F* fHistXvtxResVsPt[2*kDecays];    //!<! hist. for sec vert x resol
  TH2F* fHistYvtxResVsPt[2*kDecays];    //!<! hist. for sec vert x resol
  TH2F* fHistZvtxResVsPt[2*kDecays];    //!<! hist. for sec vert x resol
  TH2F* fHistInvMassVsPt[2*kDecays];    //!<! hist. of inv mass (meas)
  TH2F* fHistDecLenVsPt[2*kDecays];     //!<! hist. of decay length (meas)
  TH2F* fHistNormDLxyVsPt[2*kDecays];     //!<! hist. of decay length (meas)
  TH2F* fHistCosPointVsPt[2*kDecays];   //!<! hist. of cos(theta_p) (meas)

  Int_t fAODProtection;         /// flag to activate protection against AOD-dAOD mismatch.

  AliRDHFCutsD0toKpi *fRDHFCuts;             /// Cuts for event selection
  AliRDHFCutsDplustoKpipi *fRDHFCutsDplus;   /// Cuts for Dplus
  TString fPartName[kDecays];                /// Names of hadron species
  Int_t fMapTrLabel[kMaxLabel];              /// map of track labels

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskDmesonMCPerform,1); 
  /// \endcond
};
#endif
