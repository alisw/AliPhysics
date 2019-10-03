#ifndef ALIANALYSISTASKMUONAODCREATION_H
#define ALIANALYSISTASKMUONAODCREATION_H

/* $Id$ */ 

#include <TChain.h>
#include <TTree.h>
#include <TList.h>
#include <TH1.h>
#include <TClonesArray.h>

#include "TMath.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliLog.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliAODHeader.h"

class AliAnalysisTaskMuonAODCreation : public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskMuonAODCreation();
  AliAnalysisTaskMuonAODCreation(const Char_t* name);
  AliAnalysisTaskMuonAODCreation& operator= (const AliAnalysisTaskMuonAODCreation& c);
  AliAnalysisTaskMuonAODCreation(const AliAnalysisTaskMuonAODCreation& c);
  virtual ~AliAnalysisTaskMuonAODCreation();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  void     UserCreateOutputObjects();
  
 protected:
  
  TList  *fOutput;
  TTree *fTree;           //  AOD output Tree
  AliAODEvent *fOutputAOD;       //! AOD out 
  
  ClassDef(AliAnalysisTaskMuonAODCreation,1);
};

#endif
