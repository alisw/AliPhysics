AliAnalysisTaskEA* AddTaskEA(
  const char*         jetarrayname            = "Jet_AKTChargedR040_tracks_pT0150_pt_scheme", //name of jet TClones array for detector level jets
  const char*         jetarraynameMC          = "Jet_AKTChargedR040_mcparticles_pT0150_pt_scheme", //name of jet TClones array for MC particle level jets
  const char*         trackarrayname          = "tracks", //name of track TClonesArray for detector level jets
  const char*         mcpariclearrayname      = "mcparticles", //name of track TClonesArray array for MC particle level jets
  const char*         rhoname                 = "", //name of track TClonesArray for detector level jets
  const char*         mcrhoname               = "", //name of track TClonesArray array for MC particle level jets
  Double_t            jetRadius               = 0.4,  //radius of analyzed jets
  UInt_t              trigger                 = AliVEvent::kINT7,  //trigger
  Int_t               isMC                    = 0,     // 0=real data    , 1= particle+detector level simulation 
  Double_t            trackEtaWindow          = 0.9,   //pseudorapidity range for tracks
  Bool_t              useVertexCut            = kTRUE,  // vertex cut
  Bool_t              usePileUpCut            = kTRUE, // discard pile up event
  Double_t            acut                    = 0.6,   //cut on relative jet area
  const char* suffix = ""                              //SUBWAGON has to be the last parameter
){

   //typeOfData   
   TString kClusName     = "CaloClusters";
   TString kCorrClusName = "CaloClustersCorr";

   Double_t jetEtaRange   = TMath::Abs(trackEtaWindow - jetRadius);


   // #### DEFINE MANAGER AND DATA CONTAINER NAMES
   AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
   if(!manager){
      ::Error("AliAnalysisTaskEA.cxx", "No analysis manager to connect to.");
      return NULL;
   }
 
 
   //__________________________________________________________________________________
   // #### DEFINE MY ANALYSIS TASK

   TString myContName("");
   myContName = Form("JetAnalysisR%02d", TMath::Nint(jetRadius*10));
   myContName.Append(suffix);



   AliAnalysisTaskEA *task = new AliAnalysisTaskEA(myContName.Data());

   if(isMC){  //for PYTHIA
      task->SetIsPythia(kTRUE);  //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
      task->SetMakeGeneralHistograms(kTRUE); //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
   }   

   //inspired by AliAnalysisTaskEmcalQGTagging
   //_____________________________________________
   //TRACK/PARTICLE CONTAINTERS
   AliTrackContainer    *trackCont      = 0x0; //detector level track container 
   AliParticleContainer *trackContTrue  = 0x0; //mc particle container

   trackCont = task->AddTrackContainer(trackarrayname);  //detector level tracks 
   trackCont->SetMinPt(0.15); 
   trackCont->SetEtaLimits(-trackEtaWindow, trackEtaWindow);
 
   if(isMC){ 
      trackContTrue = task->AddMCParticleContainer(mcpariclearrayname); //particle level MC particles   
      trackContTrue->SetClassName("AliAODMCParticle"); 
      trackContTrue->SetMinPt(0.15); 
      trackContTrue->SetEtaLimits(-trackEtaWindow,trackEtaWindow);
   }
   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRec    = 0x0; //jet container with detector level tracks
   AliJetContainer *jetContTrue   = 0x0; //jet container with mc particles

   jetContRec   = task->AddJetContainer(jetarrayname,"TPC",jetRadius);

   if(jetContRec) { //DETECTOR LEVEL JET
      jetContRec->ConnectParticleContainer(trackCont);
      jetContRec->SetPercAreaCut(acut);
      jetContRec->SetMinPt(0.150);
      jetContRec->SetMaxTrackPt(1000);
      jetContRec->SetJetAcceptanceType(AliEmcalJet::kUser);
      jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
    }

    if(isMC){
      //AKT JETS PARTICLE LEVEL
      jetContTrue = task->AddJetContainer(jetarraynameMC,"TPC",jetRadius);
      if(jetContTrue){
         jetContTrue->ConnectParticleContainer(trackContTrue);
         jetContTrue->SetPercAreaCut(acut);
         jetContTrue->SetMinPt(0.15);
         jetContTrue->SetMaxTrackPt(1000);
         jetContTrue->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }
   }
 
   // #### Task configuration 
   task->SetMC(isMC);
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SetAcceptanceWindows(trackEtaWindow);
   task->SelectCollisionCandidates(trigger);
   task->SetExternalRhoTaskName(rhoname);
   task->SetExternalRhoTaskNameMC(mcrhoname);
   task->SetTrackContainerName(trackarrayname);
   task->SetMCParticleContainerName(mcpariclearrayname);
   task->SetJetContainerName(jetarrayname);
   task->SetMCJetContainerName(jetarraynameMC);
   task->SetUseNewCentralityEstimation(kTRUE);  //CENTRALITY

   task->SetDebugLevel(0); //No debug messages 0

   // output container
   contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChJetSpectra%s", AliAnalysisManager::GetCommonFileName(), myContName.Data()));
 
 
   // #### ADD ANALYSIS TASK
   manager->AddTask(task);
   manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
   manager->ConnectOutput(task, 1, contHistos);

   return task;
}
