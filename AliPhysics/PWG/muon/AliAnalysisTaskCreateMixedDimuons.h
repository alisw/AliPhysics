#ifndef AliAnalysisTaskCreateMixedDimuons_H
#define AliAnalysisTaskCreateMixedDimuons_H

/* $Id$ */ 

// Example of an analysis task creating aod events filled with mixed muon pairs
//
// Authors Alessandro De Falco and Antonio Uras, INFN Cagliari
// alessandro.de.falco@ca.infn.it  antonio.uras@ca.infn.it

#include "AliAnalysisTaskME.h"
#include "AliEventPoolMuon.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "TTree.h"
#include "AliEventPoolMuon.h"

//===========================================================================================

class AliAnalysisTaskCreateMixedDimuons : public AliAnalysisTaskME {

 public:
  AliAnalysisTaskCreateMixedDimuons(const char *name = "AliAnalysisTaskCreateMixedDimuons");
  virtual ~AliAnalysisTaskCreateMixedDimuons() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   SetDebug(Bool_t debug) { fDebug = debug; }
  
 private:
  AliAnalysisTaskCreateMixedDimuons(const AliAnalysisTaskCreateMixedDimuons&);
  AliAnalysisTaskCreateMixedDimuons& operator=(const AliAnalysisTaskCreateMixedDimuons&);
  Int_t              fBufferSize;
  AliAODEvent       *fInputAOD[100];        // AOD input events
  AliAODHandler     *fOutputUserHandler;    // AOD handler for the user-defined output events
  AliAODEvent       *fOutputUserAOD;
  TTree             *fOutputUserAODTree;
  AliEventPoolMuon  *fPoolMuon;
  Bool_t             fDebug;
   
  ClassDef(AliAnalysisTaskCreateMixedDimuons, 1);

};

//===========================================================================================

#endif
