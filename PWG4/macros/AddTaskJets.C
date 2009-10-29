AliJetReader *CreateJetReader(Char_t *jr); // Common config
AliJetFinder *CreateJetFinder(Char_t *jf,Float_t radius = -1);

AliAnalysisTaskJets *AddTaskJets(Char_t *jr, Char_t *jf,Float_t radius = -1); // for the new AF

AliAnalysisTaskJets *AddTaskJets(){
  // fills the standard "jets" branch in the AOD
  // need the ESDFilter to run before, to access the AODtracks
  // Tracks selected by the first Filter (1<<0)
  // needs to be adapted for different cuts
  
  // UA1 as standard chosen, since it is the most robust and simple JF
  // R = 0.4 suffficient to provide accurate jet axis for correlation studies
  // energy resolution suffers a little
  // Acceptance of jets not limited by the Jet Finder but should be done
  // by user to abs(eta) < 0.5 

  return AddTaskJets("AOD","UA1",0.4);

}



AliAnalysisTaskJets *AddTaskJets(Char_t *jr, Char_t *jf, Float_t radius)
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

  AliAODHandler *aodh = (AliAODHandler*)mgr->GetOutputEventHandler();
  if (!aodh) {
    ::Error("AddTaskJets", "This task needs an output event handler");
    return NULL;
  }   


   // Create the task and configure it.
   //===========================================================================
   AliAnalysisTaskJets *jetana;
   AliJetReader *er = CreateJetReader(jr);
    // Define jet header and jet finder
   AliJetFinder *jetFinder = CreateJetFinder(jf,radius);

   if (jetFinder){
       if (er) jetFinder->SetJetReader(er);
   }

   char *cRadius = "";
   if(radius>0)cRadius = Form("%02d",(int)(radius*10));

   jetana = new AliAnalysisTaskJets(Form("JetAnalysis%s%s%s",jr,jf,cRadius));
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (type == "AOD") jetana->SetNonStdBranch(Form("jets%s",jf));
   
   TString c_jr(jr);
   c_jr.ToLower();
   TString c_jf(jf);
   c_jf.ToLower();

   AliAnalysisDataContainer *cout_jet = mgr->CreateContainer(Form("jethist%s%s%s",c_jr.Data(),c_jf.Data(),cRadius), TList::Class(),
							     AliAnalysisManager::kOutputContainer, Form("jethist_%s_%s%s.root",c_jr.Data(),c_jf.Data(),cRadius));
   // Connect jet finder to task.
   jetana->SetJetFinder(jetFinder);
   jetana->SetConfigFile("");
   jetana->SetDebugLevel(0);
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

AliJetFinder *CreateJetFinder(Char_t *jf,Float_t radius){
  AliJetFinder *jetFinder = 0;

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


  case "FASTJET":
    AliFastJetHeaderV1 *jh = new AliFastJetHeaderV1();
    jh->SetRparam(0.4); // setup parameters                                  
    if(radius>0)jh->SetRparam(radius);
    jh->SetAlgorithm(2); // antikt from fastjet/JetDefinition.hh
    jetFinder = new AliFastJetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "UA1":
    AliUA1JetHeaderV1 *jh=new AliUA1JetHeaderV1();
    jh->SetComment("UA1 jet code with default parameters");
    jh->BackgMode(0);
    jh->SetRadius(0.4);
    if(radius>0)jh->SetRadius(radius);
    jh->SetEtSeed(4.);
    jh->SetEtSeed(4.);
    jh->SetNAcceptJets(6);
    jh->SetLegoNbinPhi(432);
    jh->SetLegoNbinEta(274);
    jh->SetLegoEtaMin(-2);
    jh->SetLegoEtaMax(+2);
    jh->SetMinJetEt(5.);
    jh->SetJetEtaMax(1.5);
    jh->SetJetEtaMin(-1.5);

    jetFinder = new AliUA1JetFinderV1();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "UA1MC":
    AliUA1JetHeaderV1 *jh=new AliUA1JetHeaderV1();
    jh->SetComment("UA1 jet code with default MC parameters");
    jh->BackgMode(0);
    jh->SetRadius(1.0);
    if(radius>0)jh->SetRadius(radius);
    jh->SetEtSeed(4.);
    jh->SetNAcceptJets(6);
    jh->SetLegoNbinPhi(432);
    jh->SetLegoNbinEta(274);
    jh->SetLegoEtaMin(-2);
    jh->SetLegoEtaMax(+2);
    jh->SetMinJetEt(5.);
    jh->SetJetEtaMax(1.5);
    jh->SetJetEtaMin(-1.5);
    jetFinder = new AliUA1JetFinderV1();
    if (jh) jetFinder->SetJetHeader(jh);
    break;
  default:
    Printf("\n >>>>>>> AddTaskJets Error Wrong jet finder selected\n");
    break;
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
  case "MC2":
    AliJetKineReaderHeader *jrh = new AliJetKineReaderHeader();
    jrh->SetComment("MC full Kinematics spearate config charged only");
    jrh->SetFastSimTPC(kFALSE);
    jrh->SetFastSimEMCAL(kFALSE);
    jrh->SetChargedOnly(kTRUE);
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
  case "AODMC":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD MC Reader");
    jrh->SetPtCut(0.);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0.9
    jrh->SetReadAODMC(1);// 1 all primary MC , 2 all primary charged
    // Define reader and set its header
    er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODMC2":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD MC Reader");
    jrh->SetPtCut(0.);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0.9
    jrh->SetReadAODMC(2);// 1 all primary MC , 2 all primary charged
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
