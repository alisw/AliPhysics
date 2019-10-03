#ifndef AliAnalysisTaskFromStandardToMuonAOD_H
#define  AliAnalysisTaskFromStandardToMuonAOD_H

/* $Id$ */ 

#include "AliAODEvent.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODEventInfo.h"

class AliAnalysisTaskFromStandardToMuonAOD : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFromStandardToMuonAOD() : AliAnalysisTaskSE() {}
  AliAnalysisTaskFromStandardToMuonAOD(const char *name);
  virtual ~AliAnalysisTaskFromStandardToMuonAOD() {}
  
  virtual void   UserCreateOutputObjects();				
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);			
  
 private:
  AliAnalysisTaskFromStandardToMuonAOD(const AliAnalysisTaskFromStandardToMuonAOD&); // Not implemented
  AliAnalysisTaskFromStandardToMuonAOD& operator=(const AliAnalysisTaskFromStandardToMuonAOD&); // Not implemented

 ClassDef(AliAnalysisTaskFromStandardToMuonAOD, 1); 
};
#endif

