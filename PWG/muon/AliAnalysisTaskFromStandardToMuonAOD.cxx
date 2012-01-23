
/* $Id$ */

#define AliAnalysisTaskFromStandardToMuonAOD_cxx

#include "TTree.h"
#include "TROOT.h"

#include "AliAODEvent.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskFromStandardToMuonAOD.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODEventInfo.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskFromStandardToMuonAOD)

//________________________________________________________________________
AliAnalysisTaskFromStandardToMuonAOD::AliAnalysisTaskFromStandardToMuonAOD(const char *name):
	AliAnalysisTaskSE(name){}

//________________________________________________________________________
void AliAnalysisTaskFromStandardToMuonAOD::UserCreateOutputObjects() {}

//________________________________________________________________________
void AliAnalysisTaskFromStandardToMuonAOD::UserExec(Option_t *) {
  if(fDebug > 1) AliInfo(Form(" Processing event # %lld",  Entry())) ; 
  return;
}      

//________________________________________________________________________
void AliAnalysisTaskFromStandardToMuonAOD::Terminate(Option_t *) {
  if(fDebug > 1) printf("AliAnalysisTaskFromStandardToMuonAOD: Terminate() \n");
}
