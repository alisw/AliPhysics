AliAnalysisTaskJetChargeFlavourTemplates* AddTaskJetChargeFlavourTemplates(
										const char * njetsData, //data jets
								    const Double_t R,
								    const char * nrhoBase,
								    const char * ntracksData,
							 	    const char *type,
								    const char *CentEst,
								    Int_t       pSel,
										AliAnalysisTaskJetChargeFlavourTemplates::JetShapeSub jetShapeSub,
								    TString     trigClass      = "",
								    TString     kEmcalTriggers = "",
								    TString     tag            = ""
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

  TString wagonName1 = Form("AliAnalysisTaskJetChargeFlavourTemplates_%s_TC%s%s",njetsData,trigClass.Data(),tag.Data());
  TString wagonName2 = Form("AliAnalysisTaskJetChargeFlavourTemplates_%s_TC%s%sTree",njetsData,trigClass.Data(),tag.Data());

  //Configure jet tagger task
  AliAnalysisTaskJetChargeFlavourTemplates *task = new AliAnalysisTaskJetChargeFlavourTemplates(wagonName1.Data());

  task->SetJetShapeSub(jetShapeSub);
  task->SetJetRadius(R);


	// This creates the Particle continer which i later collect in AliAnalysis Task Jet
  AliParticleContainer *trackContData=0x0;
  trackContData = task->AddParticleContainer(ntracksData);


	// This is me attempted to create a particle container to contain the MC particles
	//AliParticleContainer *MCContData=0x0;
	//MCcontData->SetIsPythia(kTrue);
	//MCContData = task->AddParticleContainer("mcparticles");

	//Initialising the Jet Container
  AliJetContainer *JetContData=0x0;

  TString strType(type);


	// Adds the Jet container to the Task, njetsData seems to be refering to the name of the Jet Branch.
    JetContData = task->AddJetContainer(njetsData,strType,R); //Data
    if(JetContData) {
      JetContData->SetRhoName(nrhoBase);
			//Connects the prevous Particle container to the Jet COntainer
      JetContData->ConnectParticleContainer(trackContData);
      JetContData->SetPercAreaCut(0.6);
      JetContData->SetJetRadius(R);
      JetContData->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskJetChargeFlavourTemplates::kConstSub) JetContData->SetAreaEmcCut(-2);
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
