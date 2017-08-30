AliAnalysisTaskJetPP* AddTaskJetPP(
   Bool_t	      ismc                    = kFALSE,
   const char*         hybridtracks            = "",
   const char*         mcparticless            = "",
   const char*         jetfinderAkt            = "",
   const char*         jetfinderKt             = "",
   const char*         jetfinderAktmc          = "",
   const char*         jetfinderKtmc           = "",
   //UInt_t              trigger                 = AliVEvent::kINT7,  //trigger
   UInt_t              trigger                 = AliVEvent::kINT7,  //trigger
   const char*         containerSuffix         = "",   //tag to the name of container
   const char*         centralityType          = "V0A",   //centrality
   Double_t            trackEtaWindow          = 0.9,   //pseudorapidity range for tracks
   Bool_t              useVertexCut            = kTRUE,  // vertex cut
   Bool_t              usePileUpCut            = kTRUE // discard pile up event
){

   Double_t jetAKTRadius;
   Double_t jetKTRadius;

   TString jetFinderNameAKT = jetfinderAkt;
   TString jetFinderNameKT = jetfinderKt;
   for(int i=2;i<6;i++){
      if(jetFinderNameAKT.Contains(Form("R0%d0",i))) jetAKTRadius = i*0.1;
      if(jetFinderNameKT.Contains(Form("R0%d0",i)))  jetKTRadius = i*0.1;
   }
   enum MyContainer {
      kContainerOne = 0,
      kContainerTwo   = 1,
      kContainerThree = 2,
      kContainerFour  = 3
   };

   //typeOfData   
   Double_t jetEtaRange   = TMath::Abs(trackEtaWindow - jetAKTRadius);  //fiducial cut on jet axis
   Double_t jetEtaRangeKT   = TMath::Abs(trackEtaWindow - jetKTRadius);  //fiducial cut on jet axis

   // #### Detect the demanded trigger with its readable name
   TString triggerName(Form("Trigger_%i", trigger));
   if(trigger == AliVEvent::kAnyINT)
      triggerName = "kAnyINT";
   else if(trigger == AliVEvent::kAny)
      triggerName = "kAny";
   else if(trigger == AliVEvent::kINT7)
      triggerName = "kINT7";
   else if(trigger == AliVEvent::kMB)
      triggerName = "kMB";
   else if(trigger == AliVEvent::kEMC7)
      triggerName = "kEMC7";
   else if(trigger == AliVEvent::kEMCEJE)
      triggerName = "kEMCEJE";
   else if(trigger == AliVEvent::kEMCEGA)
      triggerName = "kEMCEGA";
   else 
      triggerName = "kUNKW";

 
   // #### DEFINE MANAGER AND DATA CONTAINER NAMES
   AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
   if(!manager){
      ::Error("AddTaskJetPP.C", "No analysis manager to connect to.");
      return NULL;
   }
 
   TString containerNameSuffix("");
   if(strcmp(containerSuffix,""))  containerNameSuffix = Form("_%s", containerSuffix);
 
   TString myContName(""); //name of my container
   myContName = Form("AnalysisAKTR%02dKTR%02d_%s_%s", 
      TMath::Nint(jetAKTRadius*10),TMath::Nint(jetKTRadius*10), triggerName.Data(), containerNameSuffix.Data());

   //_________________________________________________________
   //  JET FINDER TASKS
   enum AlgoType {kkt, kantikt};
   enum JetType  {kfulljets, kchargedjets, kneutraljets};

   if(jetAKTRadius < 0.1) return NULL;
 
//__________________________________________________________________________________
   //MC Tagger match jets
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetTagger.C");
   AliAnalysisTaskEmcalJetTagger* tagr = NULL;
   //cout<<"jetfinderAkt="  <<jetfinderAkt<<endl;
   if(ismc){
      //cout <<"jetfinderAktmc =" <<jetfinderAktmc<<endl;
      tagr = AddTaskEmcalJetTagger(jetfinderAkt,jetfinderAktmc,jetAKTRadius,
                                                               "","", //nrhoBase, nrhoTa
                                                                hybridtracks, //ntracks
                                                                "", //nclusters
                                                                "TPC",
                                                                "V0A",//   centralityType,
                                                                trigger,
                                                                "",""  //trigClass, kEmcalTriggers
                                                                );


      tagr->SetJetTaggingType(AliAnalysisTaskEmcalJetTagger::kClosest); //GEN-REC JET MATCHING DONE HERE
      tagr->SetJetTaggingMethod(AliAnalysisTaskEmcalJetTagger::kGeo);
      tagr->SelectCollisionCandidates(trigger);
      tagr->SetDebugLevel(0);
      AliJetContainer *cont  = tagr->GetJetContainer(kContainerOne); //0
      AliJetContainer *cont2 = tagr->GetJetContainer(kContainerTwo);//1
      if(!cont) cout <<"AddTaskPP: Missing cont1"<<endl;
      if(!cont2) cout <<"AddTaskPP: Missing cont2"<<endl;
      cont->SetMinPt(0.150);
      cont2->SetMinPt(0.150);
      cont->SetMaxTrackPt(1000);
      cont2->SetMaxTrackPt(1000);
      cont->SetJetPhiLimits(-10.,10.);
      cont2->SetJetPhiLimits(-10.,10.);
      cont->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      cont2->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
   }

   //__________________________________________________________________________________
   // #### DEFINE MY ANALYSIS TASK

   AliAnalysisTaskJetPP *task = new AliAnalysisTaskJetPP(
                                             Form("MyJetTask_%s",myContName.Data()));
   if (ismc){
      task->SetIsPythia(kTRUE);
      task->SetMakeGeneralHistograms(kTRUE); //fill XSection and nTrials
   }

   //_____________________________________________
   //TRACK CONTAINTERS
   AliParticleContainer *trackCont      = 0x0; //reconstructed tracks

   trackCont   =  task->AddParticleContainer(hybridtracks);  //reconstructed tracks 
   trackCont->SetClassName("AliAODTrack");
   //trackCont->SetFilterHybridTracks(kTRUE);

   //_____________________________________________
   //MC PARTICLE  CONTAINTERS

   AliParticleContainer *particleCont  = 0x0; //mc particles

   if(ismc){
      particleCont = task->AddMCParticleContainer(mcparticless); //gen particles   
      particleCont->SetClassName("AliAODMCParticle");  //???????????????///
   }

   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRecAKT = 0x0; //AKT jets with reconstructed tracks
   AliJetContainer *jetContRecKT  = 0x0;     //KT jets with reconstructed tracks

   //AKT jet containers settings
   jetContRecAKT   = task->AddJetContainer(jetfinderAkt,"TPC",jetAKTRadius);
   if(jetContRecAKT) {
      jetContRecAKT->ConnectParticleContainer(trackCont);
      //jetContRecAKT->SetPercAreaCut(0.6); // cut out all jets that have area less than 0.6*pi*R^2 (typically low momentum jets)
      jetContRecAKT->SetMaxTrackPt(1000); //it is good to cut away high pT tracks (mostly tracking fakes)
      jetContRecAKT->SetMinTrackPt(0.15); //it is good to cut away high pT tracks (mostly tracking fakes)
      jetContRecAKT->SetJetPtCut(0.149);
      jetContRecAKT->SetJetAcceptanceType(AliJetContainer::kUser);
      jetContRecAKT->SetJetEtaLimits(-jetEtaRange,jetEtaRange); //fiducial cut on jet axis acceptance
   }

   //KT jet containers settings
   jetContRecKT = task->AddJetContainer(jetfinderKt,"TPC",jetKTRadius);
   if(jetContRecKT){
      jetContRecKT->ConnectParticleContainer(trackCont);
      //jetContRecKT->SetPercAreaCut(0.6);
      jetContRecKT->SetMaxTrackPt(1000);
      jetContRecKT->SetMinTrackPt(0.);
      jetContRecKT->SetJetPtCut(0.);
      jetContRecKT->SetJetAcceptanceType(AliJetContainer::kUser);
      jetContRecKT->SetJetEtaLimits(-jetEtaRangeKT,jetEtaRangeKT);  // RANGE   
   }

   //_____________________________________________________________
   //MC
   AliJetContainer *jetContTrue   = 0x0; //jets from mc particles
   AliJetContainer *jetContTrueKT = 0x0; //KT jets from mc particles

   //MC AKT jet containers settings
   if(ismc){
      jetContTrue = task->AddJetContainer(jetfinderAktmc,"TPC",jetAKTRadius);
      if(jetContTrue) {
         jetContTrue->ConnectParticleContainer(particleCont);
         //jetContTrue->SetPercAreaCut(0.6);//0.6
         jetContTrue->SetMinPt(0.15);
         jetContTrue->SetMaxTrackPt(1000);
         jetContTrue->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }
      //MC KT jet containers settings
      jetContTrueKT = task->AddJetContainer(jetfinderKtmc,"TPC",jetKTRadius);
      if(jetContTrueKT) {
         jetContTrueKT->ConnectParticleContainer(particleCont);
         jetContTrueKT->SetMinPt(0.);
         jetContTrueKT->SetMaxTrackPt(1000);
         jetContTrueKT->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContTrueKT->SetJetEtaLimits(-jetEtaRangeKT,jetEtaRangeKT);
      }
   }

 
   // #### PP Task preferences
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SetAcceptanceWindows(trackEtaWindow, jetAKTRadius);
   task->SelectCollisionCandidates(trigger);
   task->SetCentralityType(centralityType); 
   task->SetMC(ismc);
   
   task->SetDebugLevel(0); //No debug messages 0

   // output container
   contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChJetSpectra%s", AliAnalysisManager::GetCommonFileName(), myContName.Data()));
 
 
   // #### ADD ANALYSIS TASK
   manager->AddTask(task);
   manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
   manager->ConnectOutput(task, 1, contHistos);

   return task;
}
