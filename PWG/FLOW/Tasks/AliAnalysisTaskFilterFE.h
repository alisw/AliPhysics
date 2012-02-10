/////////////////////////////////////////////////////
// AliAnalysisTaskFilterFE:
// analysis task to (re)tag RFP and POI of flow event 
//////////////////////////////////////////////////////

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef ALIANALYSISTASKFILTERFE_H
#define ALIANALYSISTASKFILTERFE_H

#include "AliFlowTrackSimpleCuts.h"
#include "AliFlowEventSimple.h"

class AliAnalysisTaskSE;

class AliAnalysisTaskFilterFE : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFilterFE();
  AliAnalysisTaskFilterFE(const char *name, AliFlowTrackSimpleCuts *cutsRFP, AliFlowTrackSimpleCuts *cutsPOI);
  virtual ~AliAnalysisTaskFilterFE();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  void SetSubeventEtaRange(Double_t minA, Double_t maxA, Double_t minB, Double_t maxB)
    {this->fMinA = minA; this->fMaxA = maxA; this->fMinB = minB; this->fMaxB = maxB; }
  
 private:

  AliAnalysisTaskFilterFE(const AliAnalysisTaskFilterFE& aAnalysisTask);
  AliAnalysisTaskFilterFE& operator=(const AliAnalysisTaskFilterFE& aAnalysisTask); 

  AliFlowTrackSimpleCuts* fCutsRFP; //cuts for RFPs
  AliFlowTrackSimpleCuts* fCutsPOI; //cuts for POIs
  Double_t fMinA; //minimum of eta range for subevent A
  Double_t fMaxA; //maximum of eta range for subevent A
  Double_t fMinB; //minimum of eta range for subevent B
  Double_t fMaxB; //maximum of eta range for subevent B
  AliFlowEventSimple* fFlowEvent; //flowevent
  
  ClassDef(AliAnalysisTaskFilterFE, 1); // example of analysis
};

#endif

