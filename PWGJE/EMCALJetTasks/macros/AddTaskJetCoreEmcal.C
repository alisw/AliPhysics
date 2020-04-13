AliAnalysisTaskJetCoreEmcal* AddTaskJetCoreEmcal(
								 const char *njetsBase,
								 const char *njetsTrue,
								 const char *njetsPart,
								 const char *nRho,
						     const Double_t R,
								 const char *type,
								 Int_t jetShapeType = AliAnalysisTaskJetCoreEmcal::kData,
								 Float_t kTTminr=8,
								 Float_t kTTmaxr=9,
								 Float_t kTTmins=20,
								 Float_t kTTmaxs=50,
								 const char *listName = ""
                 ) {
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      ::Error("AddTaskJetCoreEmcal","No analysis manager found.");
      return NULL;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskJetCoreEmcal", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName1 = Form("JetCore_%s_TC",njetsBase);
  TString wagonName2 = Form("JetCore_%s_TCTree",njetsBase);
  //Configure jet tagger task
  AliAnalysisTaskJetCoreEmcal *task = new AliAnalysisTaskJetCoreEmcal(wagonName1.Data());
	task->SetJetContName(njetsBase);
	if(jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart) {
		task->SetJetContPartName(njetsPart);
		task->SetJetContTrueName(njetsTrue);
	}
	if(jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetPart || jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) task->SetJetContPartName(njetsPart);
	task->SetJetShapeType(jetShapeType);
	task->SetTTLowRef(kTTminr);
	task->SetTTUpRef(kTTmaxr);
	task->SetTTLowSig(kTTmins);
	task->SetTTUpSig(kTTmaxs);

	TString name = "JetCore";
	TString clusName = "caloClusters";

	AliParticleContainer *trackCont = 0x0;
	AliParticleContainer *trackContPartLevel = 0x0;
	AliParticleContainer *trackContTrueLevel = 0x0;

	if(jetShapeType == AliAnalysisTaskJetCoreEmcal::kMCTrue) {
		trackCont = task->AddMCParticleContainer("mcparticles");
	}
	else if (jetShapeType == AliAnalysisTaskJetCoreEmcal::kData) {
		trackCont = task->AddTrackContainer("tracks");
	}
	else if (jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart) {
		trackCont = task->AddTrackContainer("tracks");
    trackContTrueLevel = task->AddTrackContainer("tracks");
    trackContPartLevel = task->AddMCParticleContainer("mcparticles");
	}
	else if (jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPartCorr) {
		trackCont = task->AddTrackContainer("tracks");
    trackContTrueLevel = task->AddTrackContainer("tracks");
	}
	else if (jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetPart) {
		trackCont = task->AddTrackContainer("tracks");
    trackContPartLevel = task->AddMCParticleContainer("mcparticles");
	}
	else if (jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) {
		trackCont = task->AddTrackContainer("tracks");
    trackContPartLevel = task->AddTrackContainer("tracks");
	}


  task->AddClusterContainer(clusName);

  // connect jet container 
  TString typeStr = TString(type);
  AliJetContainer* jetContBase = 0x0;
  AliJetContainer* jetContTrue = 0x0;
  AliJetContainer* jetContPart = 0x0;
	if(jetShapeType == AliAnalysisTaskJetCoreEmcal::kData ||
			jetShapeType == AliAnalysisTaskJetCoreEmcal::kMCTrue ||
			jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPartCorr) { 

		jetContBase = task->AddJetContainer(njetsBase,typeStr,R);
		jetContBase->SetRhoName(nRho);
		jetContBase->ConnectParticleContainer(trackCont);
		jetContBase->SetPercAreaCut(0.0);
	}
	if(jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart) {

		jetContBase = task->AddJetContainer(njetsBase,typeStr,R);
		jetContBase->SetRhoName(nRho);
		jetContBase->ConnectParticleContainer(trackCont);
		jetContBase->SetPercAreaCut(0.0);

		jetContTrue = task->AddJetContainer(njetsTrue,typeStr,R);
		jetContTrue->SetRhoName(nRho);
		jetContTrue->ConnectParticleContainer(trackContTrueLevel);
		jetContTrue->SetPercAreaCut(0.0); 

		jetContPart = task->AddJetContainer(njetsPart,typeStr,R);
		jetContPart->SetRhoName(nRho);
		jetContPart->ConnectParticleContainer(trackContPartLevel);
		jetContPart->SetIsParticleLevel(kTRUE);
		jetContPart->SetPercAreaCut(0.0);
	}

	if(jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetPart) {

		jetContBase = task->AddJetContainer(njetsBase,typeStr,R);
		jetContBase->SetRhoName(nRho);
		jetContBase->ConnectParticleContainer(trackCont);
		jetContBase->SetPercAreaCut(0.0);

		jetContPart = task->AddJetContainer(njetsPart,typeStr,R);
		jetContPart->SetRhoName(nRho);
		jetContPart->ConnectParticleContainer(trackContPartLevel);
		jetContPart->SetIsParticleLevel(kTRUE);
		jetContPart->SetPercAreaCut(0.0);
	}

	if(jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) {

		jetContBase = task->AddJetContainer(njetsBase,typeStr,R);
		jetContBase->SetRhoName(nRho);
		jetContBase->ConnectParticleContainer(trackCont);
		jetContBase->SetPercAreaCut(0.0);

		jetContPart = task->AddJetContainer(njetsPart,typeStr,R);
		jetContPart->SetRhoName(nRho);
		jetContPart->ConnectParticleContainer(trackContPartLevel);
		jetContPart->SetPercAreaCut(0.0);
	}

  cout<<"jet shape type = "<<jetShapeType<<endl;


  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  TString contname1(name);
  TString contname2(name);
  TString contname3(name);
	contname += Form("_%s",TString(listName).Data());
	contname += Form("_%02d",Int_t(R*10+0.001));
  contname1 = contname;
  contname2 = contname;
  contname3 = contname;
  contname1 += "_histos";
  contname2 += "_embTreeInclusive";
  contname3 += "_embTreeRecoil";
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname1.Data(),
			TList::Class(),AliAnalysisManager::kOutputContainer,
			Form("%s", AliAnalysisManager::GetCommonFileName()));
	mgr->ConnectInput  (task, 0,  cinput1 );
	mgr->ConnectOutput (task, 1, coutput1 );
  if(jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart || jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetPart || jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) {
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contname2.Data(),
        TTree::Class(),AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput (task, 2, coutput2 );
  }
  if(jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart || jetShapeType == AliAnalysisTaskJetCoreEmcal::kDetPart) {
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(contname3.Data(),
        TTree::Class(),AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput (task, 3, coutput3 );
  }

  return task;

}

