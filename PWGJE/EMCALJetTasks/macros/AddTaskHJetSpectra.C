AliAnalysisTaskHJetSpectra* AddTaskHJetSpectra(
      Int_t collisionSystem        = 0,  // 0=pp,  1=pPb
      Int_t typeOfAnal             = 0,  // 0= Realdata, 1 = Eff with PYTHIA MC
      const char* jetarrayname   = "Jet_AKTChargedR040_tracks_pT0150_pt_scheme", //name of jet TClones array for detector level jets
      const char* jetarraynamePartMC = "Jet_AKTChargedR040_mcparticles_pT0150_pt_scheme", //name of jet TClones array for MC particle level jets
      const char* trackarrayname     = "tracks",            //name of track TClonesArray for detector level tracks  (or tracks in the combined event when embedding)
      const char* mcpariclearraynamePartMC = "mcparticles", //name of mcparticle TClonesArray array for MC particle level jets
      const char* ktjetarrayname     = "", //name of rho task real data
      const char* ktjetarraynamePartMC = "", //name of rho task mc data
      Double_t    jetRadius          = 0.4,  //radius of analyzed jets
      UInt_t      trigger            = AliVEvent::kAny,  //trigger
      Double_t    trackEtaWindow     = 0.9,   //pseudorapidity range for tracks
      Bool_t      useVertexCut       = kTRUE,  // vertex cut
      Bool_t      usePileUpCut       = kTRUE, // discard pile up event
      Double_t    acut               = 0.   //cut on relative jet area
   ){

   return AliAnalysisTaskHJetSpectra::AddTaskHJetSpectra(
      collisionSystem, 
      typeOfAnal,
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
      acut
    );

}

