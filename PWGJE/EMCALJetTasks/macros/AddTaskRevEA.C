PWGJE::EMCALJetTasks::AliAnalysisTaskRevEA* AddTaskRevEA(
       Int_t       mode               = 0, // analysis mode   normal=0, mc=1 or embedded=2 
       const char* jetarrayname       = "Jet_AKTChargedR040_tracks_pT0150_pt_scheme", //name of jet TClones array for detector level jets (or combined event jets when emb)
       const char* jetarraynamePartMC = "Jet_AKTChargedR040_mcparticles_pT0150_pt_scheme", //name of jet TClones array for MC particle level jets
       const char* trackarrayname     = "tracks",            //name of track TClonesArray for detector level tracks  (or tracks in the combined event when embedding)  
       const char* mcpariclearraynamePartMC = "mcparticles", //name of mcparticle TClonesArray array for MC particle level jets
       const char* ktjetarrayname     = "", //name of kT jets real data 
       const char* ktjetarraynamePartMC = "", //name of kT jets task mc data
       Double_t    jetRadius          = 0.4,  //radius of analyzed jets
       UInt_t      trigger            = AliVEvent::kAny,  //trigger
       Double_t    trackEtaWindow     = 0.9,   //pseudorapidity range for tracks
       Bool_t      useVertexCut       = kTRUE,  // vertex cut
       Bool_t      usePileUpCut       = kTRUE, // discard pile up event
       Double_t    acut               = 0.,   //cut on relative jet area
       const char* suffix = ""                              //SUBWAGON has to be the l
){

   return PWGJE::EMCALJetTasks::AliAnalysisTaskRevEA::AddTaskRevEA(
       mode, 
       jetarrayname, 
       jetarraynamePartMC,
       trackarrayname,
       mcpariclearraynamePartMC,
       ktjetarrayname,
       ktjetarraynamePartMC,
       jetRadius, 
       trigger,
       trackEtaWindow,
       useVertexCut,
       usePileUpCut,
       acut,
       suffix 
          ); 
}
