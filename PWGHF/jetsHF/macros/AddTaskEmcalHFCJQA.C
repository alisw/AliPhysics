AliAnalysisTaskEmcalHFCJQA* AddTaskEmcalHFCJQA(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Int_t       nCentBins          = 1,
  Double_t    jetradius          = 0.4,
  Double_t    jetptcut           = 0,
  Double_t    jetareacut         = 0.6,
  const char *type               = "TPC",
  Int_t       leadhadtype        = 0,
  const char *taskname           = "AliAnalysisTaskEmcalHFCJQA",
  TString cutfile				 ="HFCJCuts.root",
  UInt_t triggerMask			 =-1,//AliVEvent::kEMC1 | AliVEvent::kEMC7 | AliVEvent::kEMC8,/*kMB kEMC7 (kEMC8) kEMCEJE kEMCEGA*/
  Bool_t isMC					 = kFALSE,
  Float_t minC					 = 0.,
  Float_t maxC					 = 7.5
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalHFCJQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalHFCJQA", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(taskname);
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  if (strcmp(type,"")) {
    name += "_";
    name += type;
  }

  Printf("name: %s",name.Data());

  AliAnalysisTaskEmcalHFCJQA* jetTask = new AliAnalysisTaskEmcalHFCJQA(name);
  jetTask->SetReadMC(isMC);
  jetTask->SetDebug(-1);
  jetTask->SetFilterBit(AliAODTrack::kTrkGlobalNoDCA);
  //Defaut parameters
  
  AliParticleContainer *trackCont  = jetTask->AddParticleContainer(ntracks);
  trackCont->SetClassName("AliVTrack");
  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);

  TString strType(type);
  AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);
  if(jetCont) {
    jetCont->SetRhoName(nrho);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
    jetCont->SetZLeadingCut(0.98,0.98);
    jetCont->SetPercAreaCut(0.6);
    jetCont->SetJetPtCut(jetptcut);    
    jetCont->SetLeadingHadronType(leadhadtype);
  }

//=========================CUTS=========================
AliRDHFJetsCuts *cuts;

bool kFileExists=kFALSE;
//if(!gSystem->AccessPathName(cutfile.Data(),kFileExists)){
  
if(gSystem->AccessPathName(cutfile.Data(),kFileExists))
{
	Printf("\n==CutObject not Defined in .root File. Using Standard Cuts==\n");

// possible (!not standard!) selection for pp2012 data triggered with EMCAL
    cuts=new AliRDHFJetsCuts();
    hfTask = new AliAnalysisTaskEmcalHFCJQA("AliAnalysisTaskEmcalHFCJQA");
   
    // AliRDHFJetsCutsVertex *cuts2=new AliRDHFJetsCutsVertex("jetCuts");

    //cuts for jets
    //    cuts->SetJetRadius(0.4);
    //    cuts->SetMaxEtaJet(0.5);//0.9-R
    //    cuts->SetMinPtJet(0);
    //    cuts->SetMaxPtJet(200);
    //    cuts->ResetMaskAndEnableMBTrigger();
    //cuts->SetUseAnyTrigger();
    //cuts->SetTriggerMask(0);
    cuts->SetTriggerMask(AliVEvent::kEMC1 | AliVEvent::kEMC7 | AliVEvent::kEMC8);
    cuts->SetTriggerClass("");

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetMinNClustersITS(2);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetEtaRange(-0.8,0.8);
    esdTrackCuts->SetPtRange(1,1.e10);

    if(minC>0&&minC<maxC)
    {
      // Pb-Pb
      cuts->SetTriggerClass("");
      cuts->ResetMaskAndEnableMBTrigger();
      cuts->EnableCentralTrigger();
      cuts->EnableSemiCentralTrigger();
      cuts->SetUseCentrality(AliRDHFCuts::kCentV0M);
      cuts->SetMinCentrality(minC);
      cuts->SetMaxCentrality(maxC);  
    }

    cuts->AddTrackCuts(esdTrackCuts);
    //cuts for vertexing
    //    ::Error("AddTaskSEHFJets","No Cut Object");
}
else
{
	TFile *f=TFile::Open(cutfile.Data());
	//cuts= (AliRDHFCutsD0toKpi*)f->Get("EventTrackCuts");
	cuts= (AliRDHFJetsCuts*)f->Get("HFCJCuts");

	cout<<"\n==========================================\n Cutfile used:\n"<<cutfile.Data()<<endl;
	//cuts->PrintAll();
	//jetTask->SetJetCuts(cuts);
}
if(triggerMask>0) cuts->SetTriggerMask(triggerMask);

jetTask->SetJetCuts(cuts);
delete cuts;
//========================================================== 
  
//-------------------------------------------------------
// Final settings, pass to manager and set the containers
//-------------------------------------------------------
  
  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );
  
  return jetTask;
}

AliAnalysisTaskEmcalHFCJQA* AddTaskEmcalHFCJQA( AliEmcalJetTask* jetFinderTask,
  Int_t       nCentBins          = 1,
  Double_t    jetareacut         = 0.6,
  const char *type               = "EMCAL",
  Int_t       leadhadtype        = 0,
  const char *taskname           = "AliAnalysisTaskEmcalHFCJQA"
)
    {
    const char* ntracks            = jetFinderTask->GetTracksName();
    const char* nclusters          = jetFinderTask->GetClusName();
    const char* njets              = jetFinderTask->GetJetsName();
    const char* nrho               = jetFinderTask->GetRhoName();
    Double_t    jetradius          = jetFinderTask->GetRadius();
    Double_t    jetptcut           = jetFinderTask->GetMinJetPt();

    AliAnalysisTaskEmcalHFCJQA* jetTask = AddTaskEmcalHFCJQA(ntracks , nclusters, njets, nrho, nCentBins, jetradius, jetptcut, jetareacut, type, leadhadtype, taskname);

    return jetTask;
    }
