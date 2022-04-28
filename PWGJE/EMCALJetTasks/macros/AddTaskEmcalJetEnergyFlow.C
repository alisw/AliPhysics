  AliAnalysisTaskEmcalJetEnergyFlow* AddTaskEmcalJetEnergyFlow(
       const char *ntracks            = "usedefault",
       const char *nclusters          = "usedefault",
       const char* ncells             = "usedefault",
       Double_t Rstep_EF               = 0.1,
       Bool_t      IsMCprod           = kTRUE,
       const char *suffix             = ""
     )
     {
       return AliAnalysisTaskEmcalJetEnergyFlow::AddTaskEmcalJetEnergyFlow(ntracks,
          nclusters,
          ncells,
          Rstep_EF,                                                                 
          IsMCprod,
          suffix);
    }
 
