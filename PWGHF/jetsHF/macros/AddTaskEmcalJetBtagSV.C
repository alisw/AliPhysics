AliAnalysisTaskEmcalJetBtagSV *AddTaskEmcalJetBtagSV(const char *trackcontname = "Tracks", const char *jetcontname = "Jets", const char *mctrackcontname = "MCParticles", const char *mcjetcontname = "MCJets", Double_t jetradius = 0.2, const char *type = "TPC", TString fileout="standard", Bool_t corrections_mode=kFALSE, Bool_t pthardsetting=kFALSE, Bool_t doBkgRej = kTRUE, const char* pthardmin="", const char * eventflavor="", Double_t taggingradius=0.4, TString cutfile="HFJetVertexCuts.root",Float_t minC=0., Float_t maxC=100.,const char *taskname = "HFjetsContainer", Bool_t QAvtxStep=kFALSE)
{

  // Mailto: andrea.rossi@ts.infn.it, svallero@to.infn.it, s.lapointe@cern.ch

  Int_t last=0;

  // Get the AnalysisManager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalJetBtagSV", "No analysis manager to connect to.");
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
    str="HFJetsVertex";
  }
  else {
    str=fileout;
    str.ReplaceAll(".root","");
  }
  str.Prepend("_");

  // Configure analysis task
  AliAnalysisTaskEmcalJetBtagSV *hfTask;
  hfTask = new AliAnalysisTaskEmcalJetBtagSV(name);

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
  hfTask->SetDoBkgRejection(doBkgRej);
  
  hfTask->SetJetContName(jetcontname);
  hfTask->SetTrackContName(trackcontname);

  // choose the method to select on the flavor of the jet and, if needed, select the pt-hard bin via the minimum pt-hard
  // also set mc jet and particle container
  if(corrections_mode){
  hfTask->SetMcJetContName(mcjetcontname);
  hfTask->SetMcTrackContName(mctrackcontname);  
  hfTask->SetPtHardSelection(pthardsetting,pthardmin,eventflavor); // kTRUE for selecting events with pt-hard, event flavor can be, e.g. bbbar, ccbar, jetjet
  hfTask->SetJetTaggingRadius(taggingradius);
  }
  // output of secondary vertex QA 
  hfTask->DoSecondaryVertexQA(QAvtxStep);
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
    cuts->SetTriggerMask(AliTrigger::kAny);
    cuts->SetTriggerClass(trigClass);
    cuts->PrintAll();
    hfTask->SetCuts(cuts);
  } else {
    // define your cuts here
    DefineCutsTask(hfTask, minC, maxC, corrections_mode);
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
Bool_t DefineCutsTask(AliAnalysisTaskEmcalJetBtagSV *task, Float_t minC, Float_t maxC, Bool_t corrections_mode){

    // define cuts for task
    AliRDHFJetsCuts *cuts=new AliRDHFJetsCuts();
    // jets
    cuts->SetJetRadius(0.4); // this cut does nothing
    cuts->SetMaxEtaJet(0.5);//0.9-R
    cuts->SetMinPtJet(5.);
    cuts->SetMaxPtJet(100.);
    // Set centrality 
    cuts->SetMinCentrality(minC);
    cuts->SetMaxCentrality(maxC);
    cuts->SetUsePhysicsSelection(kFALSE);
    if (corrections_mode) cuts->SetTriggerMask(AliTrigger::kAny);
    else {cuts->SetOptPileup(1);
    cuts->ConfigurePileupCuts(5,0.8);
    // cuts->SetUsePhysicsSelection(kFALSE);
    cuts->SetTriggerClass("CINT7");
    cuts->SetTriggerMask(AliTrigger::kINT7);
    } // pPb minbias only
    task->SetCuts(cuts);
    delete cuts;
    return kTRUE;
}

//------------------------------------------------------
Bool_t DefineCutsTagger(AliHFJetsTaggingVertex *tg){
   
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
    tg->SetCuts(cuts2);

    delete esdTrackCuts;
    delete cuts2;

    return kTRUE;

}



