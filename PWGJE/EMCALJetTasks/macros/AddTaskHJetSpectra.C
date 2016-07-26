enum AlgoType {kKTxx, kANTIKTxx};
enum JetType  {kFULLJETSxx, kCHARGEDJETSxx, kNEUTRALJETSxx};

enum MyContainer {
   kContainerOne = 0,  
   kContainerTwo   = 1,  
   kContainerThree = 2, 
   kContainerFour  = 3 
};

enum MyDataType {
  kReal   = 0,  // reconstructed real data 
  kPythia = 1,  // pythia simulation 
  kHijing = 2   // hijing simulation
};
enum MyAnalType {
  kRec    = 0,  // reconstructed real data 
  kEff    = 1,  // MC true+recontructed 
  kEmb    = 2,  // embedding pythia
  kEmbSingl = 3,  // embedding single track
  kKine   = 4   // kine 
};

enum MySystem {  //collision system
  kpp    = 0,
  kpPb   = 1,
  kPbPb  = 2
};
 
//__________________________________________________________________________________
//                REAL
//__________________________________________________________________________________

AliAnalysisTaskHJetSpectra* AddTaskHJetSpectra(
  const char*         aktJetFinderName        ="",   //AKT jet finder
  const char*         ktJetFinderName         ="",    //KT jet finder
  Double_t            jetRadius               = 0.4,  //radius of analyzed jets
  UInt_t              trigger                 = AliVEvent::kINT7,  //trigger
  Int_t               collisionSystem         = 0,  // 0=pp,  1=pPb, 2= PbPb  
  Int_t               typeOfData              = 0,  // 0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
  Int_t               typeOfAnal              = 0,  // 0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine
  const char*         containerSuffix         = "",   //tag to the name of container
  const char*         centralityType          = "V0A",   //centrality
  Double_t            trackEtaWindow          = 0.9,   //pseudorapidity range for tracks
  Bool_t              useVertexCut            = kTRUE,  // vertex cut
  Bool_t              usePileUpCut            = kTRUE, // discard pile up event
  Double_t            ttLow                   = 8.0,      // trigger hardron low pT
  Double_t            ttHigh                  = 50.0,     // trigger hadron high pT
  Int_t               ttType                  = 0,        // 0= single inclusive hadron trigger, else inclusive hadron trigger
  Double_t            dphi                    = 0.6, // |Delta phi_jet, trigger|< pi-0.6
  Bool_t              binning                 = 0,    //binning of jet histograms 0=2GeV width  1=1GeV width
  Double_t            acut                    = 0.6,   //cut on relative jet area
  Int_t               recombscheme            = 5,    //recombination scheme  5=BIpt_scheme  0=E_scheme 
  Int_t               nRandCones              = 2 
){

   //ANALYSIS OF REAL DATA

   TString recoTracks  = "tracks"; //"PicoTracks";

   //_________________________________________________________
   //  JET FINDER TASKS
   if(jetRadius < 0.1) return NULL;
   if(typeOfData != kReal) return NULL;
   if(typeOfAnal != kRec) return NULL;

   Double_t jetRadiusBg = 0.4;
   Double_t jetEtaRange   = TMath::Abs(trackEtaWindow - jetRadius);
   Double_t jetEtaRangeKT = TMath::Abs(trackEtaWindow - jetRadiusBg);

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
      ::Error("AddTaskHJetSpectra.C", "No analysis manager to connect to.");
      return NULL;
   }
 
   TString containerNameSuffix("");
   if(strcmp(containerSuffix,""))  containerNameSuffix = Form("_%s", containerSuffix);

   TString cntype = centralityType;
 
   TString myContName("");
   myContName = Form("AnalysisR%02d_%s_An%d%d_%s_Dphi%02d_T%d_Ptt%d_%d_%s", 
      TMath::Nint(jetRadius*10), triggerName.Data(), typeOfData, typeOfAnal, containerNameSuffix.Data(),
      TMath::Nint(10*dphi), ttType, TMath::Nint(ttLow), TMath::Nint(ttHigh), cntype.Data());

   //__________________________________________________________________________________
   // #### DEFINE MY ANALYSIS TASK


   AliAnalysisTaskHJetSpectra *task = new AliAnalysisTaskHJetSpectra(
                                  Form("HJetSpectra_%s_%s_An%d%d_TT%d%d", 
                                  aktJetFinderName, triggerName.Data(), typeOfData, typeOfAnal,
                                  TMath::Nint(ttLow),TMath::Nint(ttHigh)));

   //inspired by AliAnalysisTaskEmcalQGTagging
   //_____________________________________________
   //TRACK/PARTICLE CONTAINTERS
   AliTrackContainer *trackCont      = 0x0; //reconstructed tracks

   trackCont   =  task->AddTrackContainer(recoTracks.Data());  //reconstructed tracks 
   //trackCont->SetClassName("AliAODTrack");
   //trackCont->SetFilterHybridTracks(kTRUE);
   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRec    = 0x0; //jets with reconstructed tracks
   AliJetContainer *jetContRecKT  = 0x0; //KT jets with reconstructed tracks

   //typeOfData   // 0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
   // typeOfAnal  // 0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine

   jetContRec   = task->AddJetContainer(aktJetFinderName,"TPC",jetRadius);

   if(jetContRec) {
      jetContRec->ConnectParticleContainer(trackCont);
      jetContRec->SetPercAreaCut(acut);//0.6
      jetContRec->SetMinPt(0.15);
      jetContRec->SetMaxTrackPt(1000);
      jetContRec->SetJetAcceptanceType(AliJetContainer::kUser);
      jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
   }

   jetContRecKT = task->AddJetContainer(ktJetFinderName,"TPC",jetRadiusBg);
   if(jetContRecKT){
      jetContRecKT->ConnectParticleContainer(trackCont);
      //jetContRecKT->SetPercAreaCut(acut);//0.6         ?????????   APPLY CUT FOR BG KT JETS
      jetContRecKT->SetMinPt(0.);
      jetContRecKT->SetMaxTrackPt(1000);
      jetContRecKT->SetJetAcceptanceType(AliJetContainer::kUser);
      jetContRecKT->SetJetEtaLimits(-jetEtaRangeKT,jetEtaRangeKT);  // RANGE   
   }

 
   // #### Task preferences
   task->SetAnalysisType(collisionSystem,typeOfData,typeOfAnal);
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SetNofRandomCones(nRandCones);
   task->SetAcceptanceWindows(trackEtaWindow, jetRadius);
   task->SelectCollisionCandidates(trigger);
   task->SetCentralityType(centralityType); 
   task->SetDoubleBinPrecision(binning); 
 
   task->SetTT(ttLow, ttHigh);
   task->SetTTType(ttType); 
   task->SetDphi(dphi); // |Delta phi_jet, trigger|< pi-0.6

   task->SetDebugLevel(0); //No debug messages 0

   // output container
   contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChJetSpectra%s", AliAnalysisManager::GetCommonFileName(), myContName.Data()));
 
 
   // #### ADD ANALYSIS TASK
   manager->AddTask(task);
   manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
   manager->ConnectOutput(task, 1, contHistos);

   return task;
}
//___________________________________________________________________________________________________________
//               MC
//___________________________________________________________________________________________________________
AliAnalysisTaskHJetSpectra* AddTaskHJetSpectra(
  const char*         aktJetFinderName        ="",   //AKT jet finder
  const char*         ktJetFinderName         ="",    //KT jet finder
  const char*         aktJetFinderNameMC      ="",   //AKT jet finder
  const char*         ktJetFinderNameMC       ="",    //KT jet finder
  Double_t            jetRadius               = 0.4,  //radius of analyzed jets
  UInt_t              trigger                 = AliVEvent::kINT7,  //trigger
  Int_t               collisionSystem         = 0,  // 0=pp,  1=pPb, 2= PbPb  
  Int_t               typeOfData              = 0,  // 0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
  Int_t               typeOfAnal              = 0,  // 0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine
  const char*         containerSuffix         = "",   //tag to the name of container
  const char*         centralityType          = "V0A",   //centrality
  Double_t            trackEtaWindow          = 0.9,   //pseudorapidity range for tracks
  Bool_t              useVertexCut            = kTRUE,  // vertex cut
  Bool_t              usePileUpCut            = kTRUE, // discard pile up event
  Double_t            ttLow                   = 8.0,      // trigger hardron low pT
  Double_t            ttHigh                  = 50.0,     // trigger hadron high pT
  Int_t               ttType                  = 0,        // 0= single inclusive hadron trigger, else inclusive hadron trigger
  Double_t            dphi                    = 0.6, // |Delta phi_jet, trigger|< pi-0.6
  Bool_t              binning                 = 0,    //binning of jet histograms 0=2GeV width  1=1GeV width
  Double_t            acut                    = 0.6,   //cut on relative jet area
  Int_t               recombscheme            = 5    //recombination scheme  5=BIpt_scheme  0=E_scheme 
){

   //MC DATA EFF

   if(typeOfAnal!=kEff) return NULL;

   TString recoTracks  = "tracks"; //"PicoTracks";
   TString mcParticles = "mcparticles"; //"MCParticlesSelected"; //MC//


   Double_t jetRadiusBg = 0.4;
   Double_t jetEtaRange   = TMath::Abs(trackEtaWindow - jetRadius);
   Double_t jetEtaRangeKT = TMath::Abs(trackEtaWindow - jetRadiusBg);

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
      ::Error("AddTaskHJetSpectra.C", "No analysis manager to connect to.");
      return NULL;
   }
 
   TString containerNameSuffix("");
   if(strcmp(containerSuffix,""))  containerNameSuffix = Form("_%s", containerSuffix);
 
   TString cntype = centralityType;

   TString myContName("");
   myContName = Form("AnalysisR%02d_%s_An%d%d_%s_Dphi%02d_T%d_Ptt%d_%d_%s", 
      TMath::Nint(jetRadius*10), triggerName.Data(), typeOfData, typeOfAnal, containerNameSuffix.Data(),
      TMath::Nint(10*dphi), ttType, TMath::Nint(ttLow), TMath::Nint(ttHigh), cntype.Data());

   //_________________________________________________________________

   //_________________________________________________________________
   //EXPLANATION OF SETTINGS  
   //typeOfData          0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
   //typeOfAnal           0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine
   //_________________________________________________________________
   //  JET FINDER TASKS
   if(jetRadius < 0.1) return NULL;
 
   //____________________________________________________________________
 
   //Tagger - find closest generator level and reconstructed level jets    
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetTagger.C");
   AliAnalysisTaskEmcalJetTagger* tagr = AddTaskEmcalJetTagger(aktJetFinderName,
                                                               aktJetFinderNameMC,
                                                               jetRadius,
                                                               "","", //nrhoBase, nrhoTa
                                                               recoTracks.Data(), //ntracks
                                                               "", //nclusters
                                                               "TPC",
                                                               centralityType,  
                                                               trigger,
                                                               "",""  //trigClass, kEmcalTriggers
                                                               ); 
         
         
    tagr->SetJetTaggingType(AliAnalysisTaskEmcalJetTagger::kClosest); //GEN-REC JET MATCHING DONE HERE
    tagr->SetJetTaggingMethod(AliAnalysisTaskEmcalJetTagger::kGeo);
    tagr->SelectCollisionCandidates(trigger);
    tagr->SetDebugLevel(0);
    AliJetContainer *cont  = tagr->GetJetContainer(kContainerOne); //0
    AliJetContainer *cont2 = tagr->GetJetContainer(kContainerTwo);//1
    cont->SetMinPt(0.15);
    cont2->SetMinPt(0.15);
    cont->SetMaxTrackPt(1000);
    cont2->SetMaxTrackPt(1000);
    cont->SetJetPhiLimits(-10.,10.); 
    cont2->SetJetPhiLimits(-10.,10.);
    cont->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
    cont2->SetJetEtaLimits(-jetEtaRange,jetEtaRange);



   //__________________________________________________________________________________
   // #### DEFINE MY ANALYSIS TASK

   AliAnalysisTaskHJetSpectra *task = new AliAnalysisTaskHJetSpectra(
                                  Form("HJetSpectra_%s_%s_An%d%d_TT%d%d", 
                                  aktJetFinderName, triggerName.Data(), typeOfData, typeOfAnal,
                                  TMath::Nint(ttLow),TMath::Nint(ttHigh)));

   if(typeOfData == kPythia){  //EFF with PYTHIA
      task->SetIsPythia(kTRUE);  //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
      task->SetMakeGeneralHistograms(kTRUE); //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
   }   

   //inspired by AliAnalysisTaskEmcalQGTagging
   //_____________________________________________
   //TRACK/PARTICLE CONTAINTERS
   AliTrackContainer      *trackCont      = 0x0; //reconstructed tracks
   AliMCParticleContainer *trackContTrue  = 0x0; //mc particles

   trackCont   =  task->AddTrackContainer(recoTracks.Data());  //reconstructed tracks 
   //trackCont->SetClassName("AliAODTrack");
   //trackCont->SetFilterHybridTracks(kTRUE);

   trackContTrue = task->AddMCParticleContainer(mcParticles.Data()); //gen particles   
   //trackContTrue->SetClassName("AliAODMCParticle");  //???????????????///
   //trackContTrue->SetFilterHybridTracks(kTRUE); //???????????????///
   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRec    = 0x0; //jets with reconstructed tracks
   AliJetContainer *jetContTrue   = 0x0; //jets from mc particles
   AliJetContainer *jetContRecKT  = 0x0; //KT jets with reconstructed tracks
   AliJetContainer *jetContTrueKT = 0x0; //KT jets from mc particles

   //typeOfData   // 0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
   // typeOfAnal  // 0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine

   // MC TRUE+REC  OR  EMBEDDED
   //AKT JET REC  ContainerOne 
   jetContRec = task->AddJetContainer(aktJetFinderName,"TPC",jetRadius);
   if(jetContRec) {
      jetContRec->ConnectParticleContainer(trackCont);
      jetContRec->SetPercAreaCut(acut);//0.6
      //jetContRec->SetPythiaInfoName("PythiaInfo");
      jetContRec->SetMinPt(0.15);
      jetContRec->SetMaxTrackPt(1000);
      jetContRec->SetJetAcceptanceType(AliJetContainer::kUser);
      jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
   }
   //AKT JETS GEN ContainerTwo
   jetContTrue = task->AddJetContainer(aktJetFinderNameMC,"TPC",jetRadius);
   if(jetContTrue) {
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut);//0.6
      //jetContTrue->SetPythiaInfoName("PythiaInfo");
      jetContTrue->SetMinPt(0.15);
      jetContTrue->SetMaxTrackPt(1000);
      jetContTrue->SetJetAcceptanceType(AliJetContainer::kUser);
      jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
   }

   //KT JET CONTAINERS FOR BG REC  ContainerThree
   jetContRecKT = task->AddJetContainer(ktJetFinderName,"TPC",jetRadius);
   if(jetContRecKT) {
      jetContRecKT->ConnectParticleContainer(trackCont);
      //jetContRecKT->SetPercAreaCut(acut);//0.6
      jetContRecKT->SetMinPt(0.);
      //jetContRecKT->SetPythiaInfoName("PythiaInfo");
      jetContRecKT->SetMaxTrackPt(1000);
      jetContRecKT->SetJetAcceptanceType(AliJetContainer::kUser);
      jetContRecKT->SetJetEtaLimits(-jetEtaRangeKT,jetEtaRangeKT);
   }
   //KT JET CONTAINERS FOR BG GEN ContainerFour
   jetContTrueKT = task->AddJetContainer(ktJetFinderNameMC,"TPC",jetRadius);
   if(jetContTrueKT) {
      jetContTrueKT->ConnectParticleContainer(trackContTrue);
      //jetContTrueKT->SetPercAreaCut(acut);//0.6
      //jetContTrueKT->SetPythiaInfoName("PythiaInfo");
      jetContTrueKT->SetMinPt(0.);
      jetContTrueKT->SetMaxTrackPt(1000);
      jetContTrueKT->SetJetAcceptanceType(AliJetContainer::kUser);
      jetContTrueKT->SetJetEtaLimits(-jetEtaRangeKT,jetEtaRangeKT);
   }


   // #### Task preferences
   task->SetAnalysisType(collisionSystem,typeOfData,typeOfAnal);
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SetNofRandomCones(1);
   task->SetAcceptanceWindows(trackEtaWindow, jetRadius);
   task->SelectCollisionCandidates(trigger);
   task->SetCentralityType(centralityType); 
   task->SetDoubleBinPrecision(binning); 
 
   task->SetTT(ttLow, ttHigh);
   task->SetTTType(ttType); 
   task->SetDphi(dphi); // |Delta phi_jet, trigger|< pi-0.6

   task->SetDebugLevel(0); //No debug messages 0

   // output container
   contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChJetSpectra%s", AliAnalysisManager::GetCommonFileName(), myContName.Data()));
 
 
   // #### ADD ANALYSIS TASK
   manager->AddTask(task);
   manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
   manager->ConnectOutput(task, 1, contHistos);

   return task;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
AliAnalysisTaskHJetSpectra* AddTaskHJetSpectra(
  Double_t            jetRadius               = 0.4,  //radius of analyzed jets
  UInt_t              trigger                 = AliVEvent::kINT7,  //trigger
  Int_t               collisionSystem         = 0,  // 0=pp,  1=pPb, 2= PbPb  
  Int_t               typeOfData              = 0,  // 0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
  Int_t               typeOfAnal              = 0,  // 0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine
  const char*         containerSuffix         = "",   //tag to the name of container
  const char*         centralityType          = "V0A",   //centrality
  Double_t            trackEtaWindow          = 0.9,   //pseudorapidity range for tracks
  Bool_t              useVertexCut            = kTRUE,  // vertex cut
  Bool_t              usePileUpCut            = kTRUE, // discard pile up event
  Double_t            ttLow                   = 8.0,      // trigger hardron low pT
  Double_t            ttHigh                  = 50.0,     // trigger hadron high pT
  Int_t               ttType                  = 0,        // 0= single inclusive hadron trigger, else inclusive hadron trigger
  Double_t            dphi                    = 0.6, // |Delta phi_jet, trigger|< pi-0.6
  Bool_t              binning                 = 0,    //binning of jet histograms 0=2GeV width  1=1GeV width
  Double_t            acut                    = 0.6,   //cut on relative jet area
  Int_t               recombscheme            = 5,    //recombination scheme  5=BIpt_scheme  0=E_scheme 
  Int_t               nRandCones              = 2,   //number of leading jets to be excluded from bg
  Int_t               pytuneEmb               = 2,    // pythia tune used for embedding
  Double_t            ptHardMinEmb            = 50., //    pt hard min in embedding
  Double_t            ptHardMaxEmb            = 1000.,//  pt hard max in embedding
  Double_t            ecmsGeVEmb              = 5020.,//  E cms  in embedding    
  Float_t             ptWeightEmb             = 0.,    //weighting power of the embedded spectrum 
  Double_t            trackeff                = 1.1   //artificial reduction of tracking efficiency
){

   //typeOfData   

   TString kTracksName   =  (typeOfAnal==kRec || typeOfAnal==kEff) ? "tracks" : "PicoTracks";
   TString kTracksNameMC =  (typeOfAnal==kEff) ?  "mcparticles" : "MCParticlesSelected"; //MC//
   TString kGenPartices  = "GenParticles"; //particles generated for embedding
   TString kClusName     = "CaloClusters";
   TString kCorrClusName = "CaloClustersCorr";

   Double_t jetRadiusBg = 0.4;
   Double_t jetEtaRange   = TMath::Abs(trackEtaWindow - jetRadius);
   Double_t jetEtaRangeKT = TMath::Abs(trackEtaWindow - jetRadiusBg);

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
      ::Error("AddTaskHJetSpectra.C", "No analysis manager to connect to.");
      return NULL;
   }
 
   TString containerNameSuffix("");
   if(strcmp(containerSuffix,""))  containerNameSuffix = Form("_%s", containerSuffix);
 
   TString cntype = centralityType;

   TString myContName("");
   myContName = Form("AnalysisR%02d_%s_An%d%d_%s_Dphi%02d_T%d_Ptt%d_%d_%s", 
      TMath::Nint(jetRadius*10), triggerName.Data(), typeOfData, typeOfAnal, containerNameSuffix.Data(),
      TMath::Nint(10*dphi), ttType, TMath::Nint(ttLow), TMath::Nint(ttHigh), cntype.Data());

   //_________________________________________________________________
   TString recoTracks  = ""; //DETECTOR LEVEL TRACKS NAME
   TString mcParticles = ""; //GENERATOR LEVEL PARTICLE NAME

   if(typeOfAnal < kKine)                    recoTracks = kTracksName.Data();
   if(typeOfAnal==kEff || typeOfAnal==kKine) mcParticles = kTracksNameMC.Data();


   //_________________________________________________________________
   //EXPLANATION OF SETTINGS  
   //typeOfData          0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
   //typeOfAnal           0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine
   //_________________________________________________________________

   if(typeOfAnal == kEmb){ //EMBEDDING   to real data
      kTracksName = "PicoTracks";
      recoTracks  = kTracksName.Data();

      gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbeddingFromGen.C");
      AliJetEmbeddingFromGenTask* embTask = AddTaskJetEmbeddingFromGen(
          0, //use Pythia as default
         ptHardMinEmb, //ptHardMin          //  ???  WHAT RANGE OF PTHARD to EMBEDD 
         ptHardMaxEmb,//ptHardMax
         ecmsGeVEmb,//Double_t        ecms 
         kGenPartices.Data(),//tracksName
         Form("JetEmbeddingFromGenTask_TT%d%d_AN%d%d_R%02d",TMath::Nint(ttLow), TMath::Nint(ttHigh),typeOfData, typeOfAnal,  TMath::Nint(10*jetRadius)),//taskName
         0.15, //const Double_t  minPt 
         1000.,//const Double_t  maxPt          
        -0.9, //const Double_t  minEta 
         0.9, //const Double_t  maxEta    
         0, //const Double_t  minPhi       
         TMath::Pi() * 2, //const Double_t  maxPhi         
         0, //const Bool_t    copyArray    before modelling ??? kTRUE 
         1, //   const Bool_t    drawQA         = kTRUE,
        "PythiaInfo", //  const char     *pythiaInfoName = ,
         ptWeightEmb, //Float_t            ptWeight        = //????
         pytuneEmb, //   Int_t               kTune            =2,
         1 //    Int_t               kColorReco     =1,
         //   Float_t            kQuench        =4.4e6,
         //   Int_t               kAnglePyquen = 2
      );
      embTask->SetChargedOnly(kTRUE);
      embTask->SetMarkMC(99999);  //Set embedded MC track label
      //?//embTask->SetCopyArray(kFALSE);
      if(trigger>0) embTask->SelectCollisionCandidates(trigger);
   }

   if(typeOfAnal == kEmbSingl){ //EMBEDDING SINGLE TRACK  to real data 
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbedding.C");
      AliJetEmbeddingTask *embSingle = AddTaskJetEmbedding(kGenPartices.Data(), "", 
                           Form("SigleTrackEmb_TT%d%d_AN%d%d_R%02d",TMath::Nint(ttLow), TMath::Nint(ttHigh),typeOfData, typeOfAnal, TMath::Nint(10*jetRadius)), 
                                       ptHardMinEmb, ptHardMaxEmb, //min pT max pT
                                       -jetEtaRange, jetEtaRange, //min Eta. max Eta  ???????????? What range
                                       0.,TMath::TwoPi(),//min phi max phi
                                       1, //ntracks
                                       0,kFALSE);
      embSingle->SetMarkMC(99999);  //SET EMBEDDED MC TRACK LABEL
      if(trigger>0) embSingle->SelectCollisionCandidates(trigger);
      //FK//?// embSingle->SetMasslessParticles(kTRUE);

      embSingle->SetCopyArray(kTRUE);
      embSingle->SetSuffix(Form("TT%d%d_AN%d%d",TMath::Nint(ttLow), TMath::Nint(ttHigh),typeOfData, typeOfAnal));
      kGenPartices = embSingle->GetOutTrackName(); 
   }
 
   if(typeOfAnal == kEmb || typeOfAnal == kEmbSingl){ 
      //__________________________________________________ 
      // Branch merger  merge reconstructed tracks with embedded pythia generated 
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskMergeBranches.C");
      AliJetModelMergeBranches* brmergTask = AddTaskMergeBranches(kTracksName.Data(),kGenPartices.Data(),
      Form("Emb_TT%d%d_AN%d%d",TMath::Nint(ttLow), TMath::Nint(ttHigh),typeOfData, typeOfAnal),"");
      brmergTask->SetCopyArray(kTRUE);
      brmergTask->SelectCollisionCandidates(trigger);

      recoTracks.Append(Form("Emb_TT%d%d_AN%d%d",TMath::Nint(ttLow), TMath::Nint(ttHigh),typeOfData, typeOfAnal)); //reconstructed level are  tracks +ebeded tracks
      mcParticles = kGenPartices.Data(); //generator level are embedded particles
   }

   //_________________________________________________________
   //  JET FINDER TASKS

   if(jetRadius < 0.1) return NULL;
 
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");

   TString apx="";
   if(trackeff<1.0) apx = Form("_EFF%03d",TMath::Nint(trackeff*100)); 
  
   TString note = Form("_TT%d%d_AN%d%d_%s%s", TMath::Nint(ttLow), TMath::Nint(ttHigh), typeOfData, typeOfAnal, cntype.Data(), apx.Data());
 
   //REAL TRACKS - JET CLUSTERIZER 
   AliEmcalJetTask* jetFinderTask = 0x0;
   AliEmcalJetTask* jetFinderRhoKT = 0x0; 

   if(typeOfAnal != kKine){
       //ANTIKT  DETECTOR LEVEL 
      jetFinderTask = AddTaskEmcalJet(recoTracks.Data(),"",kANTIKTxx,jetRadius,  kCHARGEDJETSxx,0.150,0.300,0.005,recombscheme,
         Form("JetAKT_%s", note.Data()),0,0,0);
      jetFinderTask->SetMinJetPt(0.150);

      //jetFinderTask->GetParticleContainer(0)->SetFilterHybridTracks(kTRUE);

      //KT DETECTOR LEVEL
      jetFinderRhoKT = AddTaskEmcalJet(recoTracks.Data(),"", kKTxx, jetRadiusBg, kCHARGEDJETSxx,0.150,0.300,0.005,recombscheme,  
         Form("JetKT_%s_SigR%02d", note.Data(), TMath::Nint(10*jetRadius)),0.,0,0);
      jetFinderRhoKT->SetMinJetPt(0);

      //jetFinderRhoKT->GetParticleContainer(0)->SetFilterHybridTracks(kTRUE);

      if(trackeff<1.0){  //artificial reduction of tracking efficiency for real data jet finders
         jetFinderTask->SetTrackEfficiency(trackeff);
         jetFinderRhoKT->SetTrackEfficiency(trackeff);
      }
   }

   //____________________________________________________________________

   AliEmcalJetTask* jetFinderTaskMC = NULL;
   AliEmcalJetTask* jetFinderRhoKTMC = NULL; 
 
   if( typeOfAnal == kEff || typeOfAnal == kEmb || typeOfAnal == kEmbSingl || typeOfAnal == kKine ){

      //ANTIKT GENERATOR LEVEL
      jetFinderTaskMC = AddTaskEmcalJet(mcParticles.Data(),"", kANTIKTxx, jetRadius,  kCHARGEDJETSxx,0.150,0.300,0.005,recombscheme, 
         Form("JetAKTMC_%s", note.Data()),0.,0,0); 
      jetFinderTaskMC->SetMinJetPt(0.150);

      //if(typeOfAnal == kEff) jetFinderTaskMC->GetParticleContainer(0)->SelectPhysicalPrimaries(kTRUE);

      //KT GENERATOR LEVEL
      jetFinderRhoKTMC = AddTaskEmcalJet(mcParticles.Data(),"", kKTxx,   jetRadiusBg, kCHARGEDJETSxx,0.150,0.300,0.005,recombscheme, 
         Form("JetKTMC_%s_SigR%02d", note.Data(), TMath::Nint(10*jetRadius)),0.,0,0); 
      jetFinderRhoKTMC->SetMinJetPt(0);

      //if(typeOfAnal == kEff) jetFinderRhoKTMC->GetParticleContainer(0)->SelectPhysicalPrimaries(kTRUE);
 
      if(typeOfAnal == kEff || typeOfAnal == kEmb || typeOfAnal == kEmbSingl){ //EFFICIENCY OR EMBEDDING 
         //Tagger - find closest generator level and reconstructed level jets    
         gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetTagger.C");
         AliAnalysisTaskEmcalJetTagger* tagr = AddTaskEmcalJetTagger(jetFinderTask->GetName(),
                                                                     jetFinderTaskMC->GetName(),
                                                                     jetRadius,
                                                                     "","", //nrhoBase, nrhoTa
                                                                     recoTracks.Data(), //ntracks
                                                                     "", //nclusters
                                                                     "TPC",
                                                                     centralityType,  
                                                                     trigger,
                                                                     "",""  //trigClass, kEmcalTriggers
                                                                     ); 
         
         
         tagr->SetJetTaggingType(AliAnalysisTaskEmcalJetTagger::kClosest); //GEN-REC JET MATCHING DONE HERE
         tagr->SetJetTaggingMethod(AliAnalysisTaskEmcalJetTagger::kGeo);
         //tagr->SetTypeAcceptance(3);
   //      tagr->SetIsPythia(kTRUE); //?????// 
         tagr->SelectCollisionCandidates(trigger);
         tagr->SetDebugLevel(0);
         AliJetContainer *cont  = tagr->GetJetContainer(kContainerOne); //0
         AliJetContainer *cont2 = tagr->GetJetContainer(kContainerTwo);//1
         cont->SetMinPt(0.150);
         cont2->SetMinPt(0.150);
         cont->SetMaxTrackPt(1000);
         cont2->SetMaxTrackPt(1000);
         cont->SetJetPhiLimits(-10.,10.); 
         cont2->SetJetPhiLimits(-10.,10.);
         cont->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
         cont2->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }
   }

   //__________________________________________________________________________________
   // #### DEFINE MY ANALYSIS TASK

   TString tname = "";
   if(jetFinderTaskMC) tname = jetFinderTaskMC->GetName();
   if(jetFinderTask)   tname = jetFinderTask->GetName();

   AliAnalysisTaskHJetSpectra *task = new AliAnalysisTaskHJetSpectra(
                                  Form("HJetSpectra_%s_%s_%s", 
                                  tname.Data(), triggerName.Data(), note.Data()));

   if(typeOfAnal == kKine) task->SetNeedEmcalGeom(kFALSE); //KINE
   if(typeOfAnal == kEff && typeOfData == kPythia){  //EFF with PYTHIA
      task->SetIsPythia(kTRUE);  //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
      task->SetMakeGeneralHistograms(kTRUE); //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
   }   

   //inspired by AliAnalysisTaskEmcalQGTagging
   //_____________________________________________
   //TRACK/PARTICLE CONTAINTERS
   AliTrackContainer    *trackCont      = 0x0; //reconstructed tracks
   AliParticleContainer *trackContTrue  = 0x0; //mc particles

   if(typeOfAnal != kKine){ //not filled when kine analyzed only
      trackCont   =  task->AddTrackContainer(recoTracks.Data());  //reconstructed tracks 
      //if(typeOfAnal==kRec || typeOfAnal==kEff){ 
      //   trackCont->SetClassName("AliAODTrack");
      //}else{
      //   trackCont->SetClassName("AliVTrack");  ///AliPicoTrack
      //} 
      //trackCont->SetFilterHybridTracks(kTRUE);

   }
   if(typeOfAnal > kRec){ //embedding + eff + kine 
      if(typeOfAnal == kEff){ //not filled when kine analyzed only
         trackContTrue = task->AddMCParticleContainer(mcParticles.Data()); //gen particles   
         trackContTrue->SetClassName("AliAODMCParticle");  //???????????????///
      }else{
         trackContTrue = task->AddParticleContainer(mcParticles.Data()); //gen particles   
         trackContTrue->SetClassName("AliVParticle");  //???????????????///
      }
   } 
   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRec    = 0x0; //jets with reconstructed tracks
   AliJetContainer *jetContTrue   = 0x0; //jets from mc particles
   AliJetContainer *jetContRecKT  = 0x0; //KT jets with reconstructed tracks
   AliJetContainer *jetContTrueKT = 0x0; //KT jets from mc particles

   //typeOfData   // 0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
   // typeOfAnal  // 0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine

   if(typeOfAnal == kRec){
      //REAL DATA
      jetContRec   = task->AddJetContainer(jetFinderTask->GetName(),"TPC",jetRadius);

      if(jetContRec) {
         jetContRec->ConnectParticleContainer(trackCont);
         jetContRec->SetPercAreaCut(acut);//0.6
         jetContRec->SetMinPt(0.150);
         jetContRec->SetMaxTrackPt(1000);
         jetContRec->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }

      jetContRecKT = task->AddJetContainer(jetFinderRhoKT->GetName(),"TPC",jetRadiusBg);
      if(jetContRecKT){
         jetContRecKT->ConnectParticleContainer(trackCont);
         //jetContRecKT->SetPercAreaCut(acut);//0.6         ?????????   APPLY CUT FOR BG KT JETS
         jetContRecKT->SetMinPt(0.);
         jetContRecKT->SetMaxTrackPt(1000);
         jetContRecKT->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContRecKT->SetJetEtaLimits(-jetEtaRangeKT,jetEtaRangeKT);  // RANGE   
      }

   }else if(typeOfAnal==kEff || typeOfAnal==kEmb || typeOfAnal == kEmbSingl){
      // MC TRUE+REC  OR  EMBEDDED
      //AKT JET REC  ContainerOne 
      jetContRec = task->AddJetContainer(jetFinderTask->GetName(),"TPC",jetRadius);
      if(jetContRec) {
         jetContRec->ConnectParticleContainer(trackCont);
         jetContRec->SetPercAreaCut(acut);//0.6
         //jetContRec->SetPythiaInfoName("PythiaInfo");
         jetContRec->SetMinPt(0.15);
         jetContRec->SetMaxTrackPt(1000);
         jetContRec->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }
      //AKT JETS GEN ContainerTwo
      jetContTrue = task->AddJetContainer(jetFinderTaskMC->GetName(),"TPC",jetRadius);
      if(jetContTrue) {
         jetContTrue->ConnectParticleContainer(trackContTrue);
         jetContTrue->SetPercAreaCut(acut);//0.6
         //jetContTrue->SetPythiaInfoName("PythiaInfo");
         jetContTrue->SetMinPt(0.15);
         jetContTrue->SetMaxTrackPt(1000);
         jetContTrue->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }

      //KT JET CONTAINERS FOR BG REC  ContainerThree
      jetContRecKT = task->AddJetContainer(jetFinderRhoKT->GetName(),"TPC",jetRadius);
      if(jetContRecKT) {
         jetContRecKT->ConnectParticleContainer(trackCont);
         //jetContRecKT->SetPercAreaCut(acut);//0.6
         //jetContRecKT->SetPythiaInfoName("PythiaInfo");
         jetContRecKT->SetMinPt(0.);
         jetContRecKT->SetMaxTrackPt(1000);
         jetContRecKT->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContRecKT->SetJetEtaLimits(-jetEtaRangeKT,jetEtaRangeKT);
      }
      //KT JET CONTAINERS FOR BG GEN ContainerFour
      jetContTrueKT = task->AddJetContainer(jetFinderRhoKTMC->GetName(),"TPC",jetRadius);
      if(jetContTrueKT) {
         jetContTrueKT->ConnectParticleContainer(trackContTrue);
         //jetContTrueKT->SetPercAreaCut(acut);//0.6
         //jetContTrueKT->SetPythiaInfoName("PythiaInfo");
         jetContTrueKT->SetMinPt(0.);
         jetContTrueKT->SetMaxTrackPt(1000);
         jetContTrueKT->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContTrueKT->SetJetEtaLimits(-jetEtaRangeKT,jetEtaRangeKT);
      }


   }else if(typeOfAnal == kKine){
      //ANALYZE KINE
      //AKT JET CONTAINERS  GEN  ContainterOne
      jetContTrue = task->AddJetContainer(jetFinderTaskMC->GetName(),"TPC",jetRadius);
      if(jetContTrue){
         jetContTrue->ConnectParticleContainer(trackContTrue);
         jetContTrue->SetPercAreaCut(acut);//0.6
         //jetContTrue->SetPythiaInfoName("PythiaInfo");
         jetContTrue->SetMinPt(0.15);
         jetContTrue->SetMaxTrackPt(1000);
         jetContTrue->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }
      //KT JET CONTAINERS FOR BG GEN  ContainterTwo
      jetContTrueKT = task->AddJetContainer(jetFinderRhoKTMC->GetName(),"TPC",jetRadius);
      if(jetContTrueKT){
         jetContTrueKT->ConnectParticleContainer(trackContTrue);
         //jetContTrueKT->SetPercAreaCut(acut);//0.6
         //jetContTrueKT->SetPythiaInfoName("PythiaInfo");
         jetContTrueKT->SetMinPt(0.);
         jetContTrueKT->SetMaxTrackPt(1000);
         jetContTrueKT->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContTrueKT->SetJetEtaLimits(-jetEtaRangeKT,jetEtaRangeKT);
      }
   }
 
   // #### Task preferences
   task->SetAnalysisType(collisionSystem,typeOfData,typeOfAnal);
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SetNofRandomCones(nRandCones);
   task->SetAcceptanceWindows(trackEtaWindow, jetRadius);
   task->SelectCollisionCandidates(trigger);
   task->SetCentralityType(centralityType); 
   if(typeOfAnal==kEmb || typeOfAnal == kEmbSingl){
      task->SetMinFractionShared(0.5);   //min pT fraction shared by  MC truth and rec jet
   }
   task->SetDoubleBinPrecision(binning); 
 
   task->SetTT(ttLow, ttHigh);
   task->SetTTType(ttType); 
   task->SetDphi(dphi); // |Delta phi_jet, trigger|< pi-0.6

   task->SetDebugLevel(0); //No debug messages 0

   // output container
   contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChJetSpectra%s", AliAnalysisManager::GetCommonFileName(), myContName.Data()));
 
 
   // #### ADD ANALYSIS TASK
   manager->AddTask(task);
   manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
   manager->ConnectOutput(task, 1, contHistos);

   return task;
}
