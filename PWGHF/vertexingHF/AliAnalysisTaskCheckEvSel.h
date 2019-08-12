#ifndef ALIANALYSISTASKCHECKEVSEL_H
#define ALIANALYSISTASKCHECKEVSEL_H

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

//*************************************************************************
// Class AliAnalysisTaskCheckEvSel
// Task to check event selection with D2H code
// Author: Francesco Prino
//*************************************************************************

#include "AliAnalysisTaskSE.h"

class TH1F;
class TH2F;
class TH3F;
class AliNormalizationCounter;
class AliRDHFCutsD0toKpi;

class AliAnalysisTaskCheckEvSel : public AliAnalysisTaskSE
{
public:
  
  AliAnalysisTaskCheckEvSel();
  AliAnalysisTaskCheckEvSel(Bool_t readMC, Int_t system, AliRDHFCutsD0toKpi* cuts);
  virtual ~AliAnalysisTaskCheckEvSel();
  
  virtual void UserCreateOutputObjects();
  virtual void Init(){};
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void FinishTaskOutput(){};

  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
  void SetCutOnzVertexSPD(Int_t opt) {
    if(opt>=0 && opt<=3) fCutOnzVertexSPD=opt;
    else AliError("Wrong option for cut on zVertexSPD");
  }
  void SetUseAliEventCuts(Bool_t opt=kTRUE){fUseAliEventCuts=opt;}
  void SetEnableEvPropNtuple(Bool_t dontuple) {fEnableEvPropNtuple=dontuple;}
  void SetEnableVertexNtuple(Bool_t dontuple) {fEnableVertexNtuple=dontuple;}

  AliRDHFCuts* GetCutObject() { return (AliRDHFCuts*)fAnalysisCuts;}
  
private:

  void ConfigureEvSelAxis(TAxis* ax);

  AliAnalysisTaskCheckEvSel(const AliAnalysisTaskCheckEvSel &source);
  AliAnalysisTaskCheckEvSel& operator=(const AliAnalysisTaskCheckEvSel& source);
  
  TList  *fOutput;        //!<! list send on output slot 0
  Bool_t fReadMC;         /// flag for MC
  Int_t  fSystem;         /// 0=pp, 1=Pb-Pb, 2=p-Pb
  Int_t fCutOnzVertexSPD; /// cut on zSPD vertex to remove outliers in centrality vs. tracklets (0=no cut, 1= cut at 12 cm, 2= cut on difference to z of vtx tracks, 3=cut on nsigma distance between SPD and track vertices

  Int_t fAODProtection;   /// flag to activate protection against AOD-dAOD mismatch.
                          /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
  TH1F *fHistNEvents;                    //!<! hist. for No. of events
  TH2F *fHistNEventsVsCent;              //!<! hist. for No. of events
  TH2F *fHistNEventsVsCL1;               //!<! hist. for No. of events
  TH1F *fHistWhyRej;                     //!<! hist. for No. of events
  TH2F *fHistNEventsVsWhyRej;            //!<! hist. for No. of events
  TH1F *fHistNEventsVsTime;              //!<! hist. for No. of events
  TH1F *fHistNTrackletsBeforePileup;     //!<! hist. for No. of tracklets
  TH1F *fHistNTrackletsAfterPileup;      //!<! hist. for No. of tracklets
  TH1F *fHistNCL1BeforePileup;           //!<! hist. for No. of tracklets
  TH1F *fHistNCL1AfterPileup;            //!<! hist. for No. of tracklets
  TH1F *fHistCentrality;                 //!<! hist. of centrality distribution
  TH2F *fHistCL0vsV0MCentrality;         //!<! hist. of centrality (CL0 vs V0)
  TH2F *fHistNTracksTPCoutVsV0Cent;      //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksFB4VsV0Cent;         //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksFB4EtaPosVsV0Cent;   //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksFB4EtaNegVsV0Cent;   //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksBC0VsV0Cent;         //!<! Centrality-multiplicity correl
  TH2F *fHistNTrackletsVsV0Cent;         //!<! Centrality-multiplicity correl
  TH2F *fHistNTrackletsGoldenVsV0Cent;   //!<! Centrality-multiplicity correl
  TH3F *fHistNTrackletsGoldenVsV0CentVsZvert; //!<! Centrality-multiplicity correl
  TH1F *fHistPhiTrackelts;               //!<! Control plot
  TH2F *fHistNCL0VsV0Cent;                //!<! pileup control plot
  TH2F *fHistNCL1VsV0Cent;                //!<! pileup control plot
  TH2F *fHistT0AmplVsV0Ampl;             //!<! pileup control plot
  TH2F *fHistT0AmplVsV0Cent;             //!<! pileup control plot
  TH2F *fHistT0AmplVsNCL0;               //!<! pileup control plot
  TH2F *fHistT0AmplVsCL0Cent;               //!<! pileup control plot
  TH2F *fHistNTracksTPCoutVsNTracklets;  //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksFB4VsNTracklets;     //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksBC0VsNTracksFB4;     //!<! Centrality-multiplicity correl
  TH2F* fHistZVertexSPDBeforeCuts;       //!<! z-vertex distr.
  TH2F* fHistZVertexSPDBeforeSPDCut;     //!<! z-vertex distr.
  TH2F* fHistZVertexSPDAfterCuts;        //!<! z-vertex distr.
  TH2F* fHistZVertexSPDBadTrackVert;     //!<! z-vertex distr.
  THnSparseF* fEventProp;                //!<! event properties
  TNtuple* fNtupleEvProp;                //!<! ntuple of event props
  Bool_t fEnableEvPropNtuple;            /// flag to enable ntuple
  TNtuple* fNtupleZvtxDistVsWhyRej;      //!<! ntuple of ZvtxTRK vs. ZvtxSPD vs. Ncontributors vs. whyrej flag
  Bool_t fEnableVertexNtuple;            /// flag to enable ntuple for primary vertex studies
  Bool_t fUseAliEventCuts;               /// flag to use AliEventCuts for selection
  
  AliNormalizationCounter* fCounter;  //!<!Counter for normalization

  AliRDHFCutsD0toKpi *fAnalysisCuts;  /// Cuts for candidates

  ClassDef(AliAnalysisTaskCheckEvSel,13);
};

#endif
