AliJetReader *CreateJetReader(Char_t *jr); // Common config
AliJetFinder *CreateJetFinder(Char_t *jf);

AliAnalysisTaskJets *AddTaskJets(Char_t *jr, Char_t *jf); // for the new AF
AliAnalysisTaskJets *AddTaskJets(Char_t *jr, Char_t *jf,AliAnalysisManager* mgr,AliAnalysisDataContainer *cinput); // for the old AF

AliAnalysisTaskJets *AddTaskJets(Char_t *jr, Char_t *jf)
{
  // Creates a jet finder task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJets", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskJets", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   AliAnalysisTaskJets *jetana;
   AliJetReader *er = CreateJetReader(jr);
    // Define jet header and jet finder
   AliJetFinder *jetFinder = CreateJetFinder(jf);

   if (jetFinder){
       if (er) jetFinder->SetJetReader(er);
   }

   jetana = new AliAnalysisTaskJets(Form("JetAnalysis%s%s",jr,jf));
   AliAnalysisDataContainer *cout_jet = mgr->CreateContainer(Form("jethist%s%s",jr,jf), TList::Class(),
							      AliAnalysisManager::kOutputContainer, Form("jethist%s_%s.root",jr,jf));
   // Connect jet finder to task.
   jetana->SetJetFinder(jetFinder);
   jetana->SetConfigFile("");
   jetana->SetDebugLevel(10);
   mgr->AddTask(jetana);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (jetana, 0, mgr->GetCommonInputContainer());
// AOD output slot will be used in a different way in future
   mgr->ConnectOutput (jetana, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (jetana, 1, cout_jet);
   
   return jetana;
}

AliAnalysisTaskJets *AddTaskJets(Char_t *jr, Char_t *jf,AliAnalysisManager* mgr,AliAnalysisDataContainer *cinput)
{
  // This is only for running on PROOF with the old root version 5-22-00 
  // and the older version of the AF

  // Creates a jet fider task, configures it and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================

  if (!mgr) {
    ::Error("AddTaskJets", "No analysis manager to connect to.");
    return NULL;
  }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskJets", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   AliAnalysisTaskJets *jetana;
   AliJetReader *er = CreateJetReader(jr);

    // Define jet header and jet finder
   AliJetFinder *jetFinder = CreateJetFinder(jf);

   if (jetFinder){
       if (er) jetFinder->SetJetReader(er);
   }

   jetana = new AliAnalysisTaskJets(Form("JetAnalysis%s%s",jr,jf));
   AliAnalysisDataContainer *cout_jet = mgr->CreateContainer(Form("jethist%s%s",jr,jf), TList::Class(),AliAnalysisManager::kOutputContainer, Form("jethist%s_%s.root",jr,jf));
   // Connect jet finder to task.
   jetana->SetJetFinder(jetFinder);
   jetana->SetConfigFile("");
   jetana->SetDebugLevel(10);
   mgr->AddTask(jetana);

   // Create the output containers for the data produced by the task.
   //==============================================================================
   // AOD output slot will be used in a different way in future
   // we need to have different names for each usage
   c_aodjet = mgr->CreateContainer(Form("cAODjet%s%s",jr,jf), TTree::Class(),
				   AliAnalysisManager::kExchangeContainer);
   // Connect to data containers
   mgr->ConnectInput (jetana,0, cinput );
   mgr->ConnectOutput (jetana,0, c_aodjet );
   mgr->ConnectOutput (jetana, 1, cout_jet);
   
   return jetana;
}

AliJetFinder *CreateJetFinder(Char_t *jf){

  switch (jf) {
  case "CDF":
    AliCdfJetHeader *jh = new AliCdfJetHeader();
    jh->SetRadius(0.7);
    
    jetFinder = new AliCdfJetFinder();
    jetFinder->SetOutputFile("jets.root");
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "DA":
    AliDAJetHeader *jh=new AliDAJetHeader();
    jh->SetComment("DA jet code with default parameters");
    jh->SelectJets(kTRUE);
    jh->SetNclust(10);
    
    jetFinder = new AliDAJetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;


  case "Fastjet":
    AliFastJetHeader *jh = new AliFastJetHeader();
    jh->SetRparam(0.7); // setup parameters                                  

    jetFinder = new AliFastJetFinder();
    jetFinder->SetOutputFile("jets.root");
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "UA1":
    AliUA1JetHeaderV1 *jh=new AliUA1JetHeaderV1();
    jh->SetComment("UA1 jet code with default parameters");
    jh->BackgMode(2);
    jh->SetRadius(1.0);
    jh->SetEtSeed(4.);
    jh->SetLegoNbinPhi(432);
    jh->SetLegoNbinEta(274);
    jh->SetLegoEtaMin(-2);
    jh->SetLegoEtaMax(+2);
    jh->SetJetEtaMax(0.5);
    jh->SetJetEtaMin(-0.5);
    jh->SetMinJetEt(10.);

    jetFinder = new AliUA1JetFinderV1();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case default:
    ::Error("AddTaskJets", "Wrong jet finder selected\n");
    return 0;
  }

  return jetFinder;

}

AliJetReader *CreateJetReader(Char_t *jr){
  AliJetReader *er = 0;

  switch (jr) {
  case "MC":
    AliJetKineReaderHeader *jrh = new AliJetKineReaderHeader();
    jrh->SetComment("MC full Kinematics");
    jrh->SetFastSimTPC(kFALSE);
    jrh->SetFastSimEMCAL(kFALSE);
    jrh->SetPtCut(0.);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0 .9                                                                             
    // Define reader and set its header                                     
    er = new AliJetKineReader();
    er->SetReaderHeader(jrh);
    break;

  case "ESD":
    AliJetESDReaderHeader *jrh = new AliJetESDReaderHeader();
    jrh->SetComment("Testing");
    jrh->SetFirstEvent(0);
    jrh->SetLastEvent(1000);
    jrh->SetPtCut(0.);
    jrh->SetReadSignalOnly(kFALSE);
    // Define reader and set its header                                     
    er = new AliJetESDReader();
    er->SetReaderHeader(jrh);
    break;

  case "AOD":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD Reader");
    jrh->SetPtCut(0.);
    jrh->SetTestFilterMask(1<<0);
    // Define reader and set its header
    er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
    break;

  default:
    ::Error("AddTaskJets", "Wrong jet reader selected\n");
    return 0;
  }

  return er;

}
