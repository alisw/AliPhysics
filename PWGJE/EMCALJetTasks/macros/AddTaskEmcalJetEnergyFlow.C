#if !defined(__CINT__) && !defined(__CLING__)
#include "AliAnalysisTaskEmcalJetEnergyFlow.h"
#endif

  AliAnalysisTaskEmcalJetEnergyFlow* AddTaskEmcalJetEnergyFlow(
       const char *ntracks            = "usedefault",
       const char *nclusters          = "usedefault",
       const char* ncells             = "usedefault",
       Double_t Rstep_EF              = 0.1,
       Double_t Max_match_dr          = 0.2,
       Double_t Lead_pt_cut           = 0.0,
       AliAnalysisTaskEmcalJetEnergyFlow::AnalysisType fAnType = AliAnalysisTaskEmcalJetEnergyFlow::kppData,
       const char *suffix             = ""
     )
     {
       return AliAnalysisTaskEmcalJetEnergyFlow::AddTaskEmcalJetEnergyFlow(ntracks,
          nclusters,
          ncells,
          Rstep_EF,                                                                 
          Max_match_dr,
          Lead_pt_cut,
          fAnType,
          suffix);
    }
 
