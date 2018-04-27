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
  TH2F *fHistNEventsVsWhyRej;            //!<! hist. for No. of events
  TH1F *fHistNTrackletsBeforePileup;     //!<! hist. for No. of tracklets
  TH1F *fHistNTrackletsAfterPileup;      //!<! hist. for No. of tracklets
  TH1F *fHistNCL1BeforePileup;           //!<! hist. for No. of tracklets
  TH1F *fHistNCL1AfterPileup;            //!<! hist. for No. of tracklets
  TH1F *fHistCentrality;                 //!<! hist. of centrality distribution
  TH2F *fHistNTracksTPCoutVsV0Cent;      //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksFB4VsV0Cent;         //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksBC0VsV0Cent;         //!<! Centrality-multiplicity correl
  TH2F *fHistNTrackletsVsV0Cent;         //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksTPCoutVsNTracklets;  //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksFB4VsNTracklets;     //!<! Centrality-multiplicity correl
  TH2F *fHistNTracksBC0VsNTracksFB4;     //!<! Centrality-multiplicity correl
  TH2F* fHistZVertexSPDBeforeCuts;       //!<! z-vertex distr.
  TH2F* fHistZVertexSPDBeforeSPDCut;     //!<! z-vertex distr.
  TH2F* fHistZVertexSPDAfterCuts;        //!<! z-vertex distr.
  TH2F* fHistZVertexSPDBadTrackVert;     //!<! z-vertex distr.

  AliNormalizationCounter* fCounter;  //!<!Counter for normalization

  AliRDHFCutsD0toKpi *fAnalysisCuts;  /// Cuts for candidates

  ClassDef(AliAnalysisTaskCheckEvSel,1); 
};

#endif
