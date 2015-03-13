AliAnalysisTaskSEHFJets* AddTaskSEHFJets(const char *trackcontname = "Tracks", const char *jetcontname = "Jets", const char *mctrackcontname = "MCParticles", const char *mcjetcontname = "MCJets", Double_t jetradius = 0.2, const char *type = "TPC",TString fileout="standard",Bool_t corrections_mode=kFALSE, Int_t jetflavormeth=-1, Bool_t pthardsetting=kFALSE, const char* pthardmin="", const char * evtflavor="", Double_t taggingradius=0.4, TString cutfile="HFJetVertexCuts.root",Float_t minC=0., Float_t maxC=7.5,const char *taskname = "HFjetsContainer")
{

  // Mailto: andrea.rossi@ts.infn.it, svallero@to.infn.it

  Int_t last=0;

  // Get the AnalysisManager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSEHFJets", "No analysis manager to connect to.");
    return NULL;
  }

  TString name(taskname);
  if (strcmp(pthardmin,"")) {
    name += "_";
    name += pthardmin;
  }
  TString str,containername;
  if(fileout=="standard"){
    fileout=AliAnalysisManager::GetCommonFileName();
    fileout+=":PWGHF_HFCJ_";
    fileout+="HFJetsVertex";
    // if(containerprefix!="c")fileout+=containerprefix;
    str="HFJetsVertex";
  }
  else {
    str=fileout;
    str.ReplaceAll(".root","");
  }
  str.Prepend("_");

  // Configure analysis task
    AliAnalysisTaskSEHFJets *hfTask;
  hfTask = new AliAnalysisTaskSEHFJets(name);

  AliParticleContainer *trackCont  = hfTask->AddParticleContainer(trackcontname);
  TString strType(type);
  AliJetContainer *jetCont = hfTask->AddJetContainer(jetcontname,strType,jetradius);

  if (corrections_mode){
  AliParticleContainer *mctrackCont  = hfTask->AddParticleContainer(mctrackcontname);
  TString strType(type);
  AliJetContainer *mcjetCont = hfTask->AddJetContainer(mcjetcontname,strType,jetradius);
  }
  
  // Set data or corrections mode
  hfTask->SetCorrectionsMode(corrections_mode); // kFALSE for real data

  //Set the jets branch
  // radius is hardcoded here, should be changed TODO

  // FOR EMCALJETS SLL
  hfTask->SetJetContName(jetcontname);
  hfTask->SetTrackContName(trackcontname);

  // choose the method to select on the flavor of the jet and, if needed, select the pt-hard bin via the minimum pt-hard
  // also set mc jet and particle container
  if(corrections_mode){
  hfTask->SetMcJetContName(mcjetcontname);
  hfTask->SetMcTrackContName(mctrackcontname);  
  hfTask->SetJetIdMethod(jetflavormeth); // 2 = jet flavour via mesons
  hfTask->SetPtHardSelection(pthardsetting,pthardmin,evtflavor); // kTRUE for selecting events with pt hard set
  hfTask->SetJetTaggingRadius(taggingradius);
  }
  // Define the tagger
  AliHFJetsTaggingVertex *tagger=new AliHFJetsTaggingVertex();

  // Set analysis cuts
  if(!gSystem->AccessPathName(cutfile.Data(),kFileExists)){
    // read cuts from file
    ::Info(Form("Reading cuts from file: %s", cutfile.Data()));
    TFile *f=TFile::Open(cutfile.Data());
    AliRDHFJetsCuts *cuts= (AliRDHFCutsD0toKpi*)f->Get("HFJetsCutsVertex");
    // Set centrality 
    cuts->SetMinCentrality(minC);
    cuts->SetMaxCentrality(maxC);
    cuts->SetTriggerMask(AliVEvent::kAny);
    // cuts->SetTriggerClass(trigClass);
    cuts->PrintAll();
    hfTask->SetCuts(cuts);
  } else {
    // define your cuts here
    DefineCutsTask(hfTask, minC,maxC);
    DefineCutsTagger(tagger);
  }
 
  // Add task to manager
  hfTask->SetTagger(tagger);
  mgr->AddTask(hfTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput =   mgr->GetCommonInputContainer();
  mgr->ConnectInput(hfTask,0,cinput);
 
  // All containers in a list
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(name.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(hfTask,1,coutput);
  
  delete tagger;
  return hfTask;
}

//------------------------------------------------------
void DefineCutsTask(AliAnalysisTaskSEHFJets *task, Float_t minC, Float_t maxC){

    // define cuts for task
    AliRDHFJetsCuts *cuts=new AliRDHFJetsCuts();
    // jets
    cuts->SetJetRadius(0.4); // this cut does nothing
    cuts->SetMaxEtaJet(0.5);//0.9-R
    cuts->SetMinPtJet(20.);
    cuts->SetMaxPtJet(60.);
    // Set centrality 
    cuts->SetMinCentrality(minC);
    cuts->SetMaxCentrality(maxC);
    cuts->SetTriggerMask(AliVEvent::kAny);
    // cuts->SetTriggerClass("CEMC7");
    // cuts->ApplySPDMisalignedCutPP2012(); // important for LHC12d,e,f
    task->SetCuts(cuts);
    delete cuts; 
}

//------------------------------------------------------
void DefineCutsTagger(AliHFJetsTaggingVertex *tg){
   
    // define cuts for tagger:
    // 1) AliRDHFJetsCuts: event selections, basic jet cuts (pT jet, accpetance, emcal, pT leading-track)
    // 2) AliRDHFJetsCutsVertex: cuts to reconstruct vertices
    // (nprong, pT jet, eta, R, pTmin tracks, electron ID, displacement and other cut variables)
    AliRDHFJetsCutsVertex *cuts2=new AliRDHFJetsCutsVertex("jetCuts");

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(90);
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetEtaRange(-0.8,0.8);
    esdTrackCuts->SetPtRange(1.0,1.e10);
  
    cuts2->AddTrackCuts(esdTrackCuts);
    
    // vertexing
    cuts2->SetNprongs(3);
    cuts2->SetIsElec(kFALSE); // kTRUE to select e in jet vertex  
    cuts2->SetMinPtHardestTrack(1.0);//default 0.3
    cuts2->SetSecVtxWithKF(kFALSE);//default with StrLinMinDist
    cuts2->SetImpParCut(0.);//default 0
    cuts2->SetDistPrimSec(0.);//default 0
    cuts2->SetCospCut(-1);//default -1
    cuts2->SetTriggerMask(AliVEvent::kAny);
    // cuts2->SetTriggerClass("CEMC7");
    tg->SetCuts(cuts2);

    delete esdTrackCuts;
    delete cuts2;
}


