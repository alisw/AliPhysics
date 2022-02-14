AliAnalysisTaskJetChargeFlavourTemplates* AddTaskJetChargeFlavourTemplates(
										const char * njetsData, //data jets
										const char * njetsTruth, //Truth Jets
								    const Double_t R,
								    const char * nrhoBase,
								    const char * ntracksData,
							 	    const char *type,
								    const char *CentEst,
								    Int_t       pSel,
								    TString     trigClass      = "",
								    TString     kEmcalTriggers = "",
								    TString     tag            = "",
										const char* suffix = ""
								    ) {



  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskJetChargeFlavourTemplates","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AliAnalysisTaskJetChargeFlavourTemplates", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName1 = Form("AliAnalysisTaskJetChargeFlavourTemplates_%s_TC%s%s_%s",njetsData,trigClass.Data(),tag.Data(), suffix);
  TString wagonName2 = Form("AliAnalysisTaskJetChargeFlavourTemplates_%s_TC%s%sTree_%s",njetsData,trigClass.Data(),tag.Data(), suffix);

  //Configure jet tagger task
  AliAnalysisTaskJetChargeFlavourTemplates *task = new AliAnalysisTaskJetChargeFlavourTemplates(wagonName1.Data());

  task->SetJetRadius(R);

	// This creates the Particle continer which i later collect in AliAnalysis Task Jet
  AliParticleContainer *trackContData=0x0;
	AliParticleContainer *trackContTruth=0x0;
  trackContData = task->AddParticleContainer(ntracksData);
	trackContTruth = task->AddParticleContainer("mcparticles");

	//Initialising the Jet Containers
  AliJetContainer *JetContData=0x0;
	AliJetContainer *JetContTruth=0x0;


  TString strType(type);



	// Adds the Jet container to the Task.
    JetContData = task->AddJetContainer(njetsData,strType,R); //Data
    if(JetContData) {
      JetContData->SetRhoName(nrhoBase);
      JetContData->ConnectParticleContainer(trackContData);
      JetContData->SetPercAreaCut(0.6);
      JetContData->SetJetRadius(R);
      JetContData->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
    }

		// Adds the Truth jet container to the Task.
	    JetContTruth = task->AddJetContainer(njetsTruth,strType,R); //Data
	    if(JetContTruth) {
	      JetContTruth->SetRhoName(nrhoBase);
	      JetContTruth->ConnectParticleContainer(trackContTruth);
	      JetContTruth->SetPercAreaCut(0.6);
	      JetContTruth->SetJetRadius(R);
	      JetContTruth->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
	    }


	//Setting up jet cuts stuff
  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName1(wagonName1);
  TString contName2(wagonName2);




  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);

  return task;

}
