  AliAnalysisTaskEmcalJetEnergyFlow* AddTaskEmcalJetEnergyFlow(
       const char *ntracks            = "usedefault",
       const char *nclusters          = "usedefault",
       const char* ncells             = "usedefault",
       const char *suffix             = ""
        Bool_t      IsMCprod           = kTRUE
     )
     {
       return AliAnalysisTaskEmcalJetEnergyFlow::AddTaskEmcalJetEnergyFlow(ntracks,
          nclusters,
          ncells,
          suffix,
          IsMCprod);
    }
