#if !defined(__CINT__) && !defined(__CLING__)
#include "AliAnalysisTaskEmcalJetEnergyFlow.h"
#endif

  AliAnalysisTaskEmcalJetEnergyFlow* AddTaskEmcalJetEnergyFlow(
       const char *ntracks            = "usedefault",
       const char *nclusters          = "usedefault",
       const char* ncells             = "usedefault",
       Double_t Rstep_EF              = 0.1,
       AnalysisType fAnType           = AliAnalysisTaskEmcalJetEnergyFlow::kppData,
       const char *suffix             = ""
     )
     {
       return AliAnalysisTaskEmcalJetEnergyFlow::AddTaskEmcalJetEnergyFlow(ntracks,
          nclusters,
          ncells,
          Rstep_EF,                                                                 
          fAnType,
          suffix);
    }
 
