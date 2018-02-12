AliAnalysisTaskJetCoreEmcal* AddTaskJetCoreEmcal(
								 const char *njetsBase,
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
	task->SetTTLowRef(kTTminr);
	task->SetTTUpRef(kTTmaxr);
	task->SetTTLowSig(kTTmins);
	task->SetTTUpSig(kTTmaxs);

	TString name = "JetCore";
	TString trackName = "tracks";
	TString clusName = "caloClusters";
  if (trackName == "mcparticles") {
    task->AddMCParticleContainer(trackName);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    task->AddTrackContainer(trackName);
  }
  else if (!trackName.IsNull()) {
    task->AddParticleContainer(trackName);
  }
  task->AddClusterContainer(clusName);

  // connect jet container 
  TString sRhoChName = "Rho";
  TString typeStr = TString(type);
  AliJetContainer* jetContBase = task->AddJetContainer(njetsBase,typeStr,R);
  // for Pb-Pb
  jetContBase->SetRhoName(sRhoChName);
  jetContBase->SetPercAreaCut(0.0);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
	contname += Form("_%s",TString(listName).Data());
	contname += Form("_%02d",Int_t(R*10+0.001));
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0,  cinput1 );
  mgr->ConnectOutput (task, 1, coutput1 );

  return task;

}

