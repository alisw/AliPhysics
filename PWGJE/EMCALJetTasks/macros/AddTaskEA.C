AliAnalysisTaskEA* AddTaskEA(
  Int_t       system                          = 1,  //1=ppb, 0=pp
  const char*         jetarrayname            = "Jet_AKTChargedR040_tracks_pT0150_pt_scheme", //name of jet TClones array for detector level jets
  const char*         jetarraynameMC          = "Jet_AKTChargedR040_mcparticles_pT0150_pt_scheme", //name of jet TClones array for MC particle level jets
  const char*         trackarrayname          = "tracks", //name of track TClonesArray for detector level jets
  const char*         mcpariclearrayname      = "mcparticles", //name of track TClonesArray array for MC particle level jets
  const char*         clusterarrayname        = "caloClusters", //name of EMCAL cluster TClonesArray array  for detector level
  const char*         rhoname                 = "", //name of track TClonesArray for detector level jets
  const char*         mcrhoname               = "", //name of track TClonesArray array for MC particle level jets
  Double_t            jetRadius               = 0.4,  //radius of analyzed jets
  UInt_t              trigger                 = AliVEvent::kAny,  //trigger
  Int_t               isMC                    = 0,     // 0=real data    , 1= particle+detector level simulation 
  Double_t            trackEtaWindow          = 0.9,   //pseudorapidity range for tracks
  Bool_t              useVertexCut            = kTRUE,  // vertex cut
  Bool_t              usePileUpCut            = kTRUE, // discard pile up event
  Double_t            acut                    = 0.6,   //cut on relative jet area
  Double_t            emcaltofcut             = 30e-9,   //cut on relative jet area
  const char* suffix = ""                              //SUBWAGON has to be the last parameter
){

   return AliAnalysisTaskEA::AddTaskEA(
             system,
             jetarrayname,     
             jetarraynameMC,  
             trackarrayname,  
             mcpariclearrayname, 
             clusterarrayname,  
             rhoname,          
             mcrhoname,       
             jetRadius,      
             trigger,       
             isMC,          
             trackEtaWindow, 
             useVertexCut,   
             usePileUpCut,   
             acut,
             emcaltofcut,          
             suffix
          ); 
}
