  AliAnalysisTaskEmcalJetEnergyFlow* AddTaskEmcalJetEnergyFlow(
       const char *ntracks            = "usedefault",
       const char *nclusters          = "usedefault",
       const char* ncells             = "usedefault",
       Bool_t      IsMCprod           = kTRUE,
        const char *suffix             = ""
     )
     {
       return AliAnalysisTaskEmcalJetEnergyFlow::AddTaskEmcalJetEnergyFlow(ntracks,
          nclusters,
          ncells,
          IsMCprod,
          suffix);
    }
