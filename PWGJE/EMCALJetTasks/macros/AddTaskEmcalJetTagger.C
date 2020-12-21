enum AlgoType {kKT, kANTIKT};
enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};

AliAnalysisTaskEmcalJetTagger* AddTaskEmcalJetTagger(const char * njetsBase, 
    const char * njetsTag,
    const Double_t R,
    const char * nrhoBase,
    const char * nrhoTag,
    const char * ntracks,
    const char * nclusters,
    const char *type,
    const char *CentEst,
    Int_t       pSel,
    TString     trigClass      = "",
    TString     kEmcalTriggers = "");

AliAnalysisTaskRhoBase *AttachRhoTaskTagger(TString     kPeriod             = "LHC13b",
    TString     kTracksName         = "PicoTracks",
    TString     kClusName           = "caloClustersCorr",
    Double_t    R                   = 0.4,
    Double_t    ptminTrack          = 0.15,
    Double_t    etminClus           = 0.3,
    Int_t       recombScheme        = 0,
    TString     tag                 = "Jet"
);

AliAnalysisTaskEmcalJetTagger* AddTaskEmcalJetTagger(TString     kTracksName         = "PicoTracks", 
    TString     kClusName           = "caloClustersCorr",
    Double_t    R                   = 0.4,
    Double_t    ptminTrack          = 0.15,
    Double_t    etminClus           = 0.3,
    Double_t    ptminTag            = 4.,
    Int_t       rhoType             = 1,
    const char *type                = "EMCAL",
    const char *CentEst             = "V0M",
    Int_t       pSel                = AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB,
    TString     trigClass           = "",
    TString     kEmcalTriggers      = "",
    TString     kPeriod             = "LHC11h",
    Int_t       recombScheme        = 0,
    TString     tag1                = "Jet",
    TString     tag2                = ""
) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskEmcalJetTagger","No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetTagger", "This task requires an input event handler");
    return NULL;
  }

  // #### Add necessary jet finder tasks
  TMacro jetfinderadd(Form("%s/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C", gSystem->Getenv("ALICE_PHYSICS")));

  if(tag2.IsNull()) tag2=tag1;

  AliEmcalJetTask* jetFinderTaskBase = 0x0;
  if (strcmp(type,"TPC")==0)
    jetFinderTaskBase = (AliEmcalJetTask *)jetfinderadd.Exec(Form("\"%s\", \"\", %d, %f, %d, %f, %f, 0.005, %d, %s, %f, kFALSE, kFALSE", kTracksName.Data(), kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus,recombScheme,tag1.Data(),ptminTrack));
  else if (strcmp(type,"EMCAL")==0)
    jetFinderTaskBase = (AliEmcalJetTask *)jetfinderadd.Exec(Form("\"%s\", \"%s\", %d, %f, %d, %f, %f, 0.005, %d, %s, %f, kFALSE, kFALSE", kTracksName.Data(), kClusName.Data(), kANTIKT, R, kFULLJETS, ptminTrack, etminClus,recombScheme,tag1.Data(),ptminTrack));
  jetFinderTaskBase->SelectCollisionCandidates(AliVEvent::kAny);

  AliEmcalJetTask* jetFinderTaskTag  = (AliEmcalJetTask *)jetfinderadd.Exec(Form("\"%s\", \"\", %d, %f, %d, %f, %f, 0.005, %d, %s, %f, kFALSE, kFALSE", kTracksName.Data(), kANTIKT, R, kCHARGEDJETS, ptminTag, etminClus,recombScheme,tag2.Data(),ptminTag));
  jetFinderTaskTag->SelectCollisionCandidates(AliVEvent::kAny);

  /*
  if(tag1.EqualTo("JetPythia"))
    jetFinderTaskBase->SetTrackLabelRange(99999,9999999);
  //    jetFinderTaskBase->SelectConstituents(TObject::kBitMask, 0);
  if(tag2.EqualTo("JetPythia"))
    jetFinderTaskTag->SetTrackLabelRange(99999,9999999);
    */
  //    jetFinderTaskTag->SelectConstituents(TObject::kBitMask, 0);

  TString strJetsBase = jetFinderTaskBase->GetName();
  TString strJetsTag  = jetFinderTaskTag->GetName();

  AliAnalysisTaskRhoBase *rhoTaskBase;
  AliAnalysisTaskRhoBase *rhoTaskTag;
  TString rhoNameBase = "";
  TString rhoNameTag  = "";
  if(rhoType==1) {
    rhoTaskBase = AttachRhoTaskTagger(kPeriod,kTracksName,kClusName,R,ptminTrack,etminClus,recombScheme,tag1);
    if(rhoTaskBase) {
      rhoTaskBase->SetCentralityEstimator(CentEst);  
      rhoTaskBase->SelectCollisionCandidates(AliVEvent::kAny);
      if (strcmp(type,"TPC")==0)
        rhoNameBase = rhoTaskBase->GetOutRhoName();
      if (strcmp(type,"EMCAL")==0)
        rhoNameBase = rhoTaskBase->GetOutRhoScaledName();
    }
    if(rhoTaskTag) {
      rhoTaskTag = AttachRhoTaskTagger(kPeriod,kTracksName,kClusName,R,ptminTag,0.,recombScheme,tag2);
      rhoTaskTag->SetCentralityEstimator(CentEst); 
      rhoTaskTag->SelectCollisionCandidates(AliVEvent::kAny);
      rhoNameTag  = rhoTaskTag->GetOutRhoName();
    }
  }

  //Configure jet tagger task
  AliAnalysisTaskEmcalJetTagger *task = AddTaskEmcalJetTagger(jetFinderTaskBase->GetName(),
      jetFinderTaskTag->GetName(),
      R,
      rhoNameBase,
      rhoNameTag,
      kTracksName.Data(),
      kClusName.Data(),
      type,
      CentEst,
      pSel,
      trigClass,
      kEmcalTriggers
  );

  return task;  

}

AliAnalysisTaskEmcalJetTagger* AddTaskEmcalJetTagger(const char * njetsBase, 
    const char * njetsTag,
    const Double_t R,
    const char * nrhoBase,
    const char * nrhoTag,
    const char * ntracks,
    const char * nclusters,
    const char * type,
    const char * CentEst,
    Int_t        pSel,
    TString      trigClass,
    TString      kEmcalTriggers
) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskEmcalJetTagger","No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetTagger", "This task requires an input event handler");
    return NULL;
  }

  TString wagonName = Form("JetTagger_%s_%s_TC%s",njetsBase,njetsTag,trigClass.Data());

  //Configure jet tagger task
  AliAnalysisTaskEmcalJetTagger *task = new AliAnalysisTaskEmcalJetTagger(wagonName);

  task->SetNCentBins(4);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);

  task->SetJetContainerBase(0);
  task->SetJetContainerTag(1);

  TString strType(type);
  AliJetContainer *jetContBase = task->AddJetContainer(njetsBase,strType,R);
  if(jetContBase) {
    jetContBase->SetRhoName(nrhoBase);
    jetContBase->ConnectParticleContainer(trackCont);
    jetContBase->ConnectClusterContainer(clusterCont);
    jetContBase->SetMaxTrackPt(10000.);
  }

  AliJetContainer *jetContTag = task->AddJetContainer(njetsTag,"TPC",R);
  if(jetContTag) {
    jetContTag->SetRhoName(nrhoTag);
    jetContTag->ConnectParticleContainer(trackCont);
    jetContTag->ConnectClusterContainer(clusterCont);
    jetContTag->SetMaxTrackPt(10000.);
  }
  for(Int_t i=0; i<2; i++) {
    task->SetPercAreaCut(0.6, i); //keep?
  }
  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName(wagonName);
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);

  return task;  
}

//this AddTask only creates the task, need to add manually all the rest (for subwagons)
AliAnalysisTaskEmcalJetTagger* AddTaskEmcalJetTagger(const char *wagonName){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskEmcalJetTagger","No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetTagger", "This task requires an input event handler");
    return NULL;
  }


  //Define jet tagger task
  AliAnalysisTaskEmcalJetTagger *task = new AliAnalysisTaskEmcalJetTagger(wagonName);
  mgr->AddTask(task);
  // all the configuration is missing, do it in the wagon/subwagon
  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName(wagonName);
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);

  return task;
}

//Attach rho task
AliAnalysisTaskRhoBase *AttachRhoTaskTagger(TString kPeriod,
    TString     kTracksName,
    TString     kClusName,
    Double_t    R,
    Double_t    ptminTrack,
    Double_t    etminClus,
    Int_t       recombScheme,
    TString     tag
) {

  AliAnalysisTaskRhoBase *rhoTaskBase;

  kPeriod.ToLower();

  // Add kt jet finder and rho task in case we want background subtraction
  TMacro jetfinderadd(Form("%s/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C", gSystem->Getenv("ALICE_PHYSICS")));
  TMacro rhosparseadd(Form("%s/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C", gSystem->Getenv("ALICE_PHYSICS")));
  AliEmcalJetTask *jetFinderKt;
  AliEmcalJetTask *jetFinderAKt;
  jetFinderAKt  = (AliEmcalJetTask *)jetfinderadd.Exec(Form("\"%s\", \"\", %d, %f, %d, %f, %f, 0.005, %d, %s, %f", kTracksName.Data(), kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus,recombScheme,tag.Data(),1.));
  jetFinderAKt->SelectCollisionCandidates(AliVEvent::kAny);
  if(kPeriod.EqualTo("lhc13b") || kPeriod.EqualTo("lhc13c") || kPeriod.EqualTo("lhc13d") || kPeriod.EqualTo("lhc13e") || kPeriod.EqualTo("lhc13f")) {
    jetFinderKt   = (AliEmcalJetTask *)jetfinderadd.Exec(Form("\"%s\", \"\", %d, %f, %d, %f, %f, 0.005, %d, %s, %f", kTracksName.Data(), kKT, R, kCHARGEDJETS, ptminTrack, etminClus,recombScheme,tag.Data(),0.));

    TF1 *fScale = new TF1("fScale","1.42",0.,100.); //scale factor for pPb
    AliAnalysisTaskRhoSparse *rhoTaskSparse = (AliAnalysisTaskRhoSparse *)rhosparseadd.Exec(Form(
        "\"%s\", \"%s\", \"%s\", \"%s\", \"%s\", %f, \"TPC\", 0.01, 0.15, 0, (TF1 *)%p, 0, kTRUE, \"%s\", kTRUE",
        jetFinderKt->GetName(),
        jetFinderAKt->GetName(),
        kTracksName.Data(),
        kClusName.Data(),
        Form("RhoSparseR%03d",(int)(100*R)),
        R,
        fScale,
        Form("RhoSparseR%03d",(int)(100*R))
    ));
    rhoTaskSparse->SetUseAliAnaUtils(kTRUE);
    rhoTaskBase = dynamic_cast<AliAnalysisTaskRhoBase*>(rhoTaskSparse);
  }
  else if(kPeriod.EqualTo("lhc10h") || kPeriod.EqualTo("lhc11h") ) {
    jetFinderKt   = (AliEmcalJetTask *)jetfinderadd.Exec(Form("\"%s\", \"\", %d, %f, %d, %f, %f, 0.005, %d, %s, 0.1", kTracksName.Data(), kKT, R, kCHARGEDJETS, ptminTrack, etminClus,recombScheme,tag.Data()));

    TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
    sfunc->SetParameter(2,1.76458);
    sfunc->SetParameter(1,-0.0111656);
    sfunc->SetParameter(0,0.000107296);
    TString rhoname = Form("%sRhoR%03dptmin%3.0f%s",tag.Data(),(int)(100*R),ptminTrack*1000.0,kTracksName.Data());
    Printf("rhoname: %s",rhoname.Data());
    AliAnalysisTaskRho *rhoTask = (AliAnalysisTaskRho *)rhosparseadd.Exec(Form(
        "\"%s\", \"%s\", \"%s\", \"%s\", %f, \"TPC\", 0.01, 0, (TF1 *)%p, 2, kTRUE",
        jetFinderKt->GetName(),
        kTracksName.Data(),
        kClusName.Data(),
        rhoname.Data(),
        R,
        sfunc
        ));
    rhoTask->SetHistoBins(100,0,250);

    rhoTaskBase = dynamic_cast<AliAnalysisTaskRhoBase*>(rhoTask);
  }
  jetFinderKt->SelectCollisionCandidates(AliVEvent::kAny);
  /*
  if(tag.EqualTo("JetPythia")) {
    jetFinderKt->SetTrackLabelRange(99999,9999999);
    jetFinderAKt->SetTrackLabelRange(99999,9999999);
  }
  */

  return rhoTaskBase;
}
