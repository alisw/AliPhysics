AliAnalysisTaskHJetSpectra* AddTaskHJetSpectra(
  Double_t            jetRadius               = 0.4,  //radius of analyzed jets
  UInt_t              trigger                 = AliVEvent::kINT7,  //trigger
  Int_t               collisionSystem         = 0,  // 0=pp,  1=pPb, 2= PbPb  
  Int_t               typeOfData              = 0,  // 0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
  Int_t               typeOfAnal              = 0,  // 0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine
  Double_t            perpConeR             = 0.4,    //perp cone for deltaPt + perp cone bg
  const char*         containerSuffix         = "",   //tag to the name of container
  const char*         centralityType          = "V0A",   //centrality
  Float_t             centMin                 = 0,      //lower centrality percentil 
  Float_t             centMax                 = 100,    //upper centrality percentil 
  Double_t            trackEtaWindow          = 0.9,   //pseudorapidity range for tracks
  Bool_t              useVertexCut            = kTRUE,  // vertex cut
  Bool_t              usePileUpCut            = kTRUE, // discard pile up event
  Double_t            ttLow                   = 8.0,      // trigger hardron low pT
  Double_t            ttHigh                  = 50.0,     // trigger hadron high pT
  Int_t               ttType                  = 0,        // 0= single inclusive hadron trigger, else inclusive hadron trigger
  Double_t            dphi                    = 0.6, // |Delta phi_jet, trigger|< pi-0.6
  Bool_t              binning                 = 0,    //binning of jet histograms 0=2GeV width  1=1GeV width
  Double_t            acut                    = 0.6   //cut on relative jet area 
){

   //typeOfData   

   TString kTracksName   = "PicoTracks";
   TString kTracksNameMC = "MCParticlesSelected"; //MC//
   TString kGenPartices  = "GenParticles"; //particles generated for embedding
   TString kClusName     = "CaloClusters";
   TString kCorrClusName = "CaloClustersCorr";


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


  enum MyContainer {
     kContainerOne = 0,  
     kContainerTwo   = 1  
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
 
   // #### DEFINE MANAGER AND DATA CONTAINER NAMES
   AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
   if(!manager){
      ::Error("AddTaskHJetSpectra.C", "No analysis manager to connect to.");
      return NULL;
   }
 
   TString containerNameSuffix("");
   if(strcmp(containerSuffix,""))  containerNameSuffix = Form("_%s", containerSuffix);
 
   TString myContName("");
   myContName = Form("AnalysisR%02d_%s_An%d%d_%s_cent%.0f_%.0f_Dphi%02d_T%d_Ptt%d_%d", 
      TMath::Nint(jetRadius*10), triggerName.Data(), typeOfData, typeOfAnal, containerNameSuffix.Data(),
      centMin, centMax,
      TMath::Nint(10*dphi), ttType, TMath::Nint(ttLow), TMath::Nint(ttHigh));

   //_________________________________________________________________
   TString recoTracks  = "";// kTracksName.Data();
   TString mcParticles = "";

   if(typeOfAnal < kKine)                    recoTracks = kTracksName.Data();
   if(typeOfAnal==kEff || typeOfAnal==kKine) mcParticles = kTracksNameMC.Data();


   //_________________________________________________________________
   //EXPLANATION OF SETTINGS  
   //typeOfData          0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
   //typeOfAnal           0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine
   //_________________________________________________________________

   if(typeOfAnal == kEmb){ //EMBEDDING   to real data (or also to sim?)

      gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbeddingFromGen.C");
      AliJetEmbeddingFromGenTask* embTask = AddTaskJetEmbeddingFromGen(
          0, //use Pythia as default
         50., //ptHardMin          //  ???  WHAT RANGE OF PTHARD to EMBEDD 
         1000.,//ptHardMax
         5020.,//Double_t        ecms 
         kGenPartices.Data(),//tracksName
         "JetEmbeddingFromGenTask",//taskName
         0.15, //const Double_t  minPt 
         1000.,//const Double_t  maxPt          
        -0.9, //const Double_t  minEta 
         0.9, //const Double_t  maxEta    
         0, //const Double_t  minPhi       
         TMath::Pi() * 2, //const Double_t  maxPhi         
         0, //const Bool_t    copyArray    before modelling ??? kTRUE 
         1, //   const Bool_t    drawQA         = kTRUE,
        "PythiaInfo", //  const char     *pythiaInfoName = ,
         1, //Float_t            ptWeight        = //????
         2, //   Int_t               kTune            =2,
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
      AliJetEmbeddingTask *embSingle = AddTaskJetEmbedding(kGenPartices.Data(), "", "SigleTrackEmb", 
                                       5.,120., //min pT max pT
                                       -0.5,0.5, //min Eta. max Eta`
                                       0.,TMath::TwoPi(),//min phi max phi
                                       1, //ntracks
                                       0,kFALSE);
      embSingle->SetMarkMC(99999);  //Set embedded MC track label
      if(trigger>0) embSingle->SelectCollisionCandidates(trigger);
      //FK// embSingle->SetMasslessParticles(kTRUE);
   }
 
   if(typeOfAnal == kEmb || typeOfAnal == kEmbSingl){ 
      //__________________________________________________ 
      // Branch merger  merge reconstructed tracks with embedded pythia generated 
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskMergeBranches.C");
      AliJetModelMergeBranches* brmergTask = AddTaskMergeBranches(kTracksName.Data(),kGenPartices.Data(),"Emb","");
      brmergTask->SetCopyArray(kTRUE);
      brmergTask->SelectCollisionCandidates(trigger);

      recoTracks.Append("Emb");          //reconstructed level are  tracks +ebeded tracks
      mcParticles = kGenPartices.Data(); //generator level are embedded particles
   }

   //_________________________________________________________
   //  JET FINDER TASKS
   enum AlgoType {kKT, kANTIKT};
   enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};

   if(jetRadius < 0.1) return NULL;
 
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
   //REAL TRACKS - JET CLUSTERIZER 
   AliEmcalJetTask* jetFinderTask = 0x0;
   AliEmcalJetTask* jetFinderRho   = 0x0; 
   AliEmcalJetTask* jetFinderRhoKT = 0x0; 
   TString myRhoName="";
   //Float_t jetPtThresholdRho = 100.; //GEV cut on jet pT to be excluded from bg median calculus 
                               //??? HOW LARGE THIS CUT SHOULD BE ???
                       //MULTUIPLE HARD SCATTERINGS MAY HAPPEN
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");

   if(typeOfAnal != kKine){ 
      jetFinderTask = AddTaskEmcalJet(recoTracks.Data(),"",kANTIKT,jetRadius,  kCHARGEDJETS,0.150,0.300); //FK

      //EXTERN CMS RHO TASK for real tracks
      myRhoName = "ExternalRhoTask";
      jetFinderRho   = AddTaskEmcalJet(recoTracks.Data(),"", kANTIKT, 0.4, kCHARGEDJETS,0.150,0.300,0.005,1,"JetAKT"); // anti-kt
      jetFinderRhoKT = AddTaskEmcalJet(recoTracks.Data(),"", kKT,     0.4, kCHARGEDJETS,0.150,0.300,0.005,1,"JetKT",0.,0,0); // kt
      jetFinderRhoKT->SetMinJetPt(0);

      AliAnalysisTaskRhoSparse* rhotask = AddTaskRhoSparse(jetFinderRhoKT->GetName(), 
                                                        jetFinderRho->GetName(), 
                                                        recoTracks.Data(),   //pico trakcs
                                                                "",   //calo clusters
                                                        myRhoName.Data(), 
                                                               0.4,  //jet radius
                                                             "TPC",  //cut type
                                                              0.01,  //jet area cut
                                                               0.0,  //jet pt cut ????????
                                                                 0,  //enareacut 
                                                                 0,  //sfunc
                                                                 2,  //excl Jets  //FK// How many jets to exclude ????????
                                                            kFALSE,   //no histo
                                                   myRhoName.Data(),  //task name
                                                             kTRUE); //claculate rho CMS

   }
   // #### DEFINE EXTERN CMS RHO TASK for true MC particles 
   AliEmcalJetTask* jetFinderTaskMC = NULL;
   TString myRhoNameMC = "";
   AliEmcalJetTask* jetFinderRhoMC   = NULL; 
   AliEmcalJetTask* jetFinderRhoKTMC = NULL; 
   AliAnalysisTaskRhoSparse* rhotaskMC =NULL; 
 
   if( typeOfAnal == kEff || typeOfAnal == kEmb || typeOfAnal == kEmbSingl || typeOfAnal == kKine ){ //EFFICIENCY OR EMBEDDING ?????????? KINE 
      //clusterizer for MC true particles
      jetFinderTaskMC = AddTaskEmcalJet(mcParticles.Data(),"", kANTIKT, jetRadius,  kCHARGEDJETS,0.150,0.300,0.005,1,"JetMC"); //FK

   // #### DEFINE EXTERN CMS RHO TASK for real tracks
      myRhoNameMC = "ExternalRhoTaskMC";
      jetFinderRhoMC   = AddTaskEmcalJet(mcParticles.Data(),"", kANTIKT, 0.4, kCHARGEDJETS,0.150,0.300,0.005,1,"JetAKTMC"); // anti-kt
      jetFinderRhoKTMC = AddTaskEmcalJet(mcParticles.Data(),"", kKT,     0.4, kCHARGEDJETS,0.150,0.300,0.005,1,"JetKTMC",0.,0,0); // kt
      jetFinderRhoKTMC->SetMinJetPt(0);

      //gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
      rhotaskMC = AddTaskRhoSparse(jetFinderRhoKTMC->GetName(), 
                                   jetFinderRhoMC->GetName(), 
                                   mcParticles.Data(),   //pico trakcs
                                   "",   //calo clusters
                                   myRhoNameMC.Data(), 
                                   0.4,  //jet radius
                                  "TPC",  //cut type
                                    0.,  //jet area cut
                                    0.0,  //jet pt cut ????????
                                     0,  //enareacut 
                                     0,  //sfunc
                                     2,  //excl Jets  //FK// ????????
                                kFALSE,   //no histo
                                   myRhoNameMC.Data(),  //task name
                                kTRUE); //claculate rho CMS

  
      if(typeOfAnal == kKine) rhotaskMC->SetNeedEmcalGeom(kFALSE); //KINE

 
      if( typeOfAnal == kEff || typeOfAnal == kEmb || typeOfAnal == kEmbSingl){ //EFFICIENCY OR EMBEDDING 
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
                                                                     ); //FK
         
         
         tagr->SetJetTaggingType(AliAnalysisTaskEmcalJetTagger::kClosest); //GEN-REC JET MATCHING DONE HERE
         tagr->SetJetTaggingMethod(AliAnalysisTaskEmcalJetTagger::kGeo);
         //tagr->SetTypeAcceptance(3);
   //      tagr->SetIsPythia(kTRUE); //?????// 
         tagr->SelectCollisionCandidates(trigger);
         tagr->SetDebugLevel(0);
         AliJetContainer *cont  = tagr->GetJetContainer(kContainerOne); //0
         AliJetContainer *cont2 = tagr->GetJetContainer(kContainerTwo);//1
         cont->SetMaxTrackPt(1000);
         cont2->SetMaxTrackPt(1000);
         cont->SetJetPhiLimits(0.,10.); //FK//
         cont2->SetJetPhiLimits(0.,10.);//FK//
      }
   }

   //__________________________________________________________________________________
   // #### DEFINE ANALYSIS TASK

   TString tname = "";
   if(jetFinderTaskMC) tname = jetFinderTaskMC->GetName();
   if(jetFinderTask)   tname = jetFinderTask->GetName();

   AliAnalysisTaskHJetSpectra *task = new AliAnalysisTaskHJetSpectra(
                                  Form("HJetSpectra_%s_%s_TT%d%d", tname.Data(), triggerName.Data(),TMath::Nint(ttLow),TMath::Nint(ttHigh))); //KINE

   if(typeOfAnal == kKine) task->SetNeedEmcalGeom(kFALSE); //KINE
   //inspired by AliAnalysisTaskEmcalQGTagging
   //_____________________________________________
   //TRACK/PARTICLE CONTAINTERS
   AliParticleContainer *trackCont      = 0x0; //reconstructed tracks
   AliParticleContainer *trackContTrue  = 0x0; //mc particles
   if(typeOfAnal != kKine){ //not filled when kine analyzed only
      trackCont   =  task->AddParticleContainer(recoTracks.Data());  //reconstructed tracks 
   }
   if(typeOfAnal > kRec){ //embedding + eff + kine 
      trackContTrue = task->AddParticleContainer(mcParticles.Data()); //gen particles    
   } 
   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRec =0x0; //jets with reconstructed tracks
   AliJetContainer *jetContTrue=0x0; //jets from mc particles
   Double_t jetEtaRange = TMath::Abs(trackEtaWindow - jetRadius);
   //typeOfData   // 0 real data, 1=pythia 2=hijing trigger,if MC handler should be used
   // typeOfAnal  // 0= Realdata, 1 = Eff with MC, 2 = delta pT Embedding, 3=pure kine
   if(typeOfAnal == kRec){
      //REAL DATA
      jetContRec = task->AddJetContainer(jetFinderTask->GetName(),"TPC",jetRadius);
      if(jetContRec) {
         jetContRec->SetRhoName(myRhoName.Data());
         jetContRec->ConnectParticleContainer(trackCont);
         jetContRec->SetPercAreaCut(acut);//0.6
         jetContRec->SetMaxTrackPt(1000);
         jetContRec->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }
   }else if(typeOfAnal==kEff || typeOfAnal==kEmb || typeOfAnal == kEmbSingl){
      // MC TRUE+REC  OR  EMBEDDED
      jetContRec = task->AddJetContainer(jetFinderTask->GetName(),"TPC",jetRadius);
      if(jetContRec) {
         jetContRec->SetRhoName(myRhoName.Data());
         jetContRec->ConnectParticleContainer(trackCont);
         jetContRec->SetPercAreaCut(acut);//0.6
         jetContRec->SetPythiaInfoName("PythiaInfo");
         jetContRec->SetMaxTrackPt(1000);
         jetContRec->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }

      jetContTrue = task->AddJetContainer(jetFinderTaskMC->GetName(),"TPC",jetRadius);
      if(jetContTrue) {
         jetContTrue->SetRhoName(myRhoNameMC.Data());
         jetContTrue->ConnectParticleContainer(trackContTrue);
         jetContTrue->SetPercAreaCut(acut);//0.6
         jetContTrue->SetPythiaInfoName("PythiaInfo");
         jetContTrue->SetMaxTrackPt(1000);
         jetContTrue->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }
   }else if(typeOfAnal == kKine){
      //ANALYZE KINE
      jetContTrue = task->AddJetContainer(jetFinderTaskMC->GetName(),"TPC",jetRadius);
      if(jetContTrue){
         jetContTrue->SetRhoName(myRhoNameMC.Data());
         jetContTrue->ConnectParticleContainer(trackContTrue);
         jetContTrue->SetPercAreaCut(acut);//0.6
         jetContTrue->SetPythiaInfoName("PythiaInfo");
         jetContTrue->SetMaxTrackPt(1000);
         jetContTrue->SetJetAcceptanceType(AliJetContainer::kUser);
         jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }
   }
 
   // #### Task preferences
   task->SetAnalysisType(collisionSystem,typeOfData,typeOfAnal);
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SetPerpConeRadius(perpConeR);
   task->SetNofRandomCones(1);
   task->SetAcceptanceWindows(trackEtaWindow, jetRadius);
   //task->SetSignalJetMinArea(acut*jetRadius*jetRadius*TMath::Pi()); //The same used by Mata
   task->SelectCollisionCandidates(trigger);
   task->SetCentralityType(centralityType,centMin,centMax); 
   task->SetExternalRhoTaskName(myRhoName.Data());
   task->SetExternalRhoTaskNameMC(myRhoNameMC.Data());
   if(typeOfAnal==kEff || typeOfAnal==kEmb || typeOfAnal == kEmbSingl){
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
