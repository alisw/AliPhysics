Bool_t  kBackgroundMode = 0; // 0/1 =  bkg subtraction (off/on)

AliAnalysisDataContainer* cin_jet_cont = 0x0;

void SetJetFinderExchangeContainer(AliAnalysisDataContainer* acont) {
  cin_jet_cont = acont;

}

AliAnalysisDataContainer* GetJetFinderExchangeContainer() {return cin_jet_cont;}

AliJetFinder *CreateJetFinder(AliAnalysisDataContainer* contname,Char_t *jf,Float_t radius = -1);

AliAnalysisTaskJetsFinder *AddTaskJetsFinder(AliAnalysisDataContainer* contname,Char_t *jf,Float_t radius = -1,Int_t iBack = 0); // for the new AF
Int_t AddTaskJetsFinderDelta(AliAnalysisDataContainer* contname,char *nonStdFile = "",Bool_t kUseAODMC = kTRUE,UInt_t runFlag = 1|4|32|128|256);     
AliAnalysisTaskJetsFinder *AddTaskJetsFinder(AliAnalysisDataContainer* contname);

AliAnalysisTaskJetsFinder *AddTaskJetsFinder(const char* contname,Char_t *jf = "", Float_t radius = 0.4,Int_t iBack = 0) // LEGO trains
{
  AliAnalysisDataContainer* cont = ((AliAnalysisDataContainer*)AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject(contname));

  if (jf != "") {
    return AddTaskJetsFinder(cont,jf,radius,iBack); 
  } else {
    return AddTaskJetsFinder(cont);
  }

}

AliAnalysisTaskJetsFinder *AddTaskJetsFinder(AliAnalysisDataContainer* contname){
  // Reads the standard input "jets" branch in the AOD
  // UA1 as standard chosen, since it is the most robust and simple JF
  // R = 0.4 suffficient to provide accurate jet axis for correlation studies
  // energy resolution suffers a little
  // Acceptance of jets not limited by the Jet Finder but should be done
  // by user to abs(eta) < 0.5 


  return AddTaskJetsFinder(contname,"UA1",0.4);


}



Int_t AddTaskJetsFinderDelta(AliAnalysisDataContainer* contname,char *nonStdFile,Bool_t kUseAODMC,UInt_t runFlag){

  // Adds a whole set of jet finders  all to be written
  // to a delta AOD
  
  // this can in principle be done also with on the fly 
  // if we have an ESD input jet fidner task does automatically fetch the ouput aod

  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskJetsFinderDelta", "No analysis manager to connect to.");
    return 0;
  }  

  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetsFinderDelta", "This task requires an input event handler");
    return 0;
  }


  AliAODHandler *aodh = (AliAODHandler*)mgr->GetOutputEventHandler();
  if (!aodh) {
    ::Error("AddTaskJetsFinderDelta", "This task needs an output event handler");
    return 0;
  }   


  TString outFile(nonStdFile);


  AliAnalysisTaskJetsFinder *taskjetsFinder = 0;
  Int_t iCount = 0;

  // Jet Fidners Selected by run flag first bit 2^0 second by 2^1 etc
  const char *cJF[14]      = {"UA1","UA1","UA1", "CDF", "CDF","DA","DA","SISCONE","SISCONE","FASTJET","FASTJET","FASTKT","FASTKT","UA1LO"};
  const Float_t radius[14] = { 0.4,  0.7,  1.0,  0.4,   0.7,  0.4, 0.7,     0.4,     0.7,      0.4,       0.7,     0.4,     0.7,    0.4};
  UInt_t  flag[14]         = {  7,    7,    7,    7,      7,   7,   7,        7,       7,        7,         7,       7,       7,      7};

  // flag[5] = 0; // set siscone to 0 for proof mode...
  // flag first bit AOD, second bit AODMC2 third bit AODMC2 third (8) bit AOODMC2b (limited acceptance)
  // i.e. 7 all, 6 only MC2 and MC
  // this stay at three

  


  //  for(int i = 0; i< 17;i++){
  for(int i = 0; i< 14;i++){
    if(!(runFlag&(1<<i)))continue;
    cout << ", cJF[" << i << "]: " << cJF[i] << ", radius[" << i << "]: " << radius[i] << endl;
    taskjetsFinder = AddTaskJetsFinder(contname,cJF[i],radius[i]);
    //	if(i==4 || i==5) taskjetsFinder = AddTaskJetsFinder(cJF[i],radius[i]);
    //	if(kMergeAODs) 
    // To be commented for no pythia/hijing merging
    //taskjetsFinder->ReadAODFromOutput();
    //	taskjetsFinder->SetDetConfig(1);
    //	if(i==4 || i==5) taskjetsFinder->SetDetConfig(3);
    //	else taskjetsFinder->SetDetConfig(1);
    cout << "taskjetsFinder: " << taskjetsFinder << endl;
    if(taskjetsFinder){

      char *cRadius = "";
      if(radius[i]>0)cRadius = Form("%02d",(int)((radius[i]+0.01)*10.)); // add an offset beacuse of precision
      cout << "cRadius: " << cRadius << endl;
      //	  taskjetsFinder->SetNonStdBranch(Form("jets%s_%s%s",cReader[ib],cJF[i],cRadius)); // done in addtask jets
      if(outFile.Length()>0)taskjetsFinder->SetNonStdOutputFile(outFile.Data());
      iCount++;
    }
  }
  Printf("Added %d JetFinders",iCount);
  return 0;
}

AliAnalysisTaskJetsFinder *AddTaskJetsFinder(AliAnalysisDataContainer* contname,Char_t *jf, Float_t radius,Int_t iBack)
{
  // Creates a jet finder task, configures it and adds it to the analysis manager.
  kBackgroundMode = iBack;

  SetJetFinderExchangeContainer(contname);

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskJetsFinder", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetsFinder", "This task requires an input event handler");
    return NULL;
  }

  AliAODHandler *aodh = (AliAODHandler*)mgr->GetOutputEventHandler();
  if (!aodh) {
    ::Error("AddTaskJetsFinder", "This task needs an output event handler");
    return NULL;
  }   


  // Create the task and configure it.
  //===========================================================================
  
  Bool_t  kIsreaderOptionsValid = 0;
  TString readerOptions =  mgr->GetGlobalStr(cin_jet_cont->GetName(), kIsreaderOptionsValid);
  if(!kIsreaderOptionsValid) {  ::Error("AddTaskJetsFinder", "This task needs a to be associated with a Jet Reader Task."); }

  AliAnalysisTaskJetsFinder *taskjetsFinder;
  // Define jet header and jet finder
  AliJetFinder *jetFinder = CreateJetFinder(contname,jf,radius);

  TString jr = readerOptions.Tokenize("#")->At(0)->GetName();
  TString cAdd = "";
  cAdd += Form("%02d",(int)((radius+0.01)*10.));
  cAdd += Form("_B%d",(int)((kBackgroundMode)));
  cAdd += readerOptions.Tokenize("#")->At(1)->GetName(); //returns Form("_Filter%05d_Cut%05d",filterMask,(int)((1000.*ptTrackMin)))
  Printf("%s",cAdd.Data());

  taskjetsFinder = new AliAnalysisTaskJetsFinder(Form("JetAnalysis%s_%s%s",jr.Data(),jf,cAdd.Data()));

  TString c_jr(jr);
  c_jr.ToLower();
  TString c_jf(jf);
  c_jf.ToLower();
  

  TString bName =  Form("jets%s_%s%s",jr.Data(),jf,cAdd.Data());
  taskjetsFinder->SetNonStdBranch(bName.Data());
  Printf("Set jet branchname \"%s\"",bName.Data());


  // Connect jet finder to task.
  taskjetsFinder->SetJetFinder(jetFinder);
  taskjetsFinder->SetConfigFile("");
  taskjetsFinder->SetDebugLevel(-1);
  taskjetsFinder->SetFilterPt(0.1);
  //taskjetsFinder->SetBookAODBackground(1);
  mgr->AddTask(taskjetsFinder);
 
  cin_jet_cont = GetJetFinderExchangeContainer();
  if (! cin_jet_cont) { ::Error("AddTaskJetFinder", "This task needs an Exchange container and none were explicitly defined, taking the default one"); cin_jet_cont = mgr->GetContainers()->At(2); }
  if (! cin_jet_cont) { ::Error("AddTaskJetFinder", "This task needs an Exchange container"); return NULL;}

  AliAnalysisDataContainer *cout_jet = mgr->CreateContainer(Form("jethist_%s_%s%s",c_jr.Data(),c_jf.Data(),cAdd.Data()), TList::Class(),
							    AliAnalysisManager::kOutputContainer, Form("%s:PWGJE_jethist_%s_%s%s",
                                                            AliAnalysisManager::GetCommonFileName(),c_jr.Data(),c_jf.Data(),cAdd.Data()));
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  mgr->ConnectInput  (taskjetsFinder, 0, mgr->GetCommonInputContainer());
  mgr->ConnectInput  (taskjetsFinder, 1, cin_jet_cont);
  mgr->ConnectOutput (taskjetsFinder, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (taskjetsFinder, 1, cout_jet);
   

  return taskjetsFinder;
}

AliJetFinder *CreateJetFinder(AliAnalysisDataContainer* contname,Char_t *jf,Float_t radius){
  AliJetFinder *jetFinder = 0;

  switch (jf) {
    // CDF Jet finder ------------------------------------------------
  case "CDF":
    AliCdfJetHeader *jh = new AliCdfJetHeader();
    jh->SetRadius(0.7);
    if(radius>0)jh->SetRadius(radius);    
    jh->SetAODwrite(kTRUE);
    jh->SetAODtracksWrite(kTRUE);
    jh->SetDebug(-1);
    jh->SetAnalyseJets(kTRUE);
    jetFinder = new AliCdfJetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

    // DA Jet finder ------------------------------------------------
  case "DA":
    AliDAJetHeader *jh=new AliDAJetHeader();
    jh->SetComment("DA jet code with default parameters");
    jh->SelectJets(kTRUE);
    jh->SetDebug(-1);
    //  jh->SetNeff(200);
    //  jh->SetEtaEff(2.2);
    if(radius>0)jh->SetRadius(radius);
    jh->SetEtMin(5.);
    jh->SetFiducialEtaMin(-0.9);
    jh->SetFiducialEtaMax(0.9);
    jetFinder = new AliDAJetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

    // Anti Kt Jet finder ------------------------------------------------
  case "FASTJET":
    // ANTI KT
    AliFastJetHeaderV1 *jh = new AliFastJetHeaderV1();
    jh->SetRparam(0.4); // setup parameters                                  
    if(radius>0)jh->SetRparam(radius);
    jh->SetJetEtaMax(0.9 - radius);
    jh->SetJetEtaMin(-0.9 + radius);
    jh->SetJetPhiMin( 0.);
    jh->SetJetPhiMax(2.*TMath::Pi());
    jh->SetAlgorithm(2); // antikt from fastjet/JetDefinition.hh
    jh->SetDebug(-1);
    /*
      $FASTJET/include/fastjet/JetDefinition.hh
      enum JetAlgorithm {kt_algorithm, cambridge_algorithm,
      antikt_algorithm, genkt_algorithm,
      ee_kt_algorithm, ee_genkt_algorithm, ...};

    */
    jetFinder = new AliFastJetFinder();
    jh->SetPtMin(5);
    // Background
    jh->SetBGMode(kBackgroundMode);
    Double_t rBkg = 0.4;
    jh->SetBGAlgorithm(1);
    jh->SetRparamBkg(rBkg);
    jh->SetGhostEtaMax(0.9);
    jh->SetGhostArea(0.005);                            // area of a ghost
    jh->SetRapRange( -0.9+rBkg, 0.9-rBkg);              // rapidity range for subtracting background must be < ghostmaxrap-0.95*R
    jh->SetPhiRange(0 , 2*TMath::Pi());                 // phi range for subtracting background

    // Additional background
    jh->SetBkgFastJetb(1);
    jh->SetBkgFastJetWoHardest(1); 
    if (jh) jetFinder->SetJetHeader(jh);
    break;

    // Fast Kt Jet finder ------------------------------------------------
  case "FASTKT":
    AliFastJetHeaderV1 *jh = new AliFastJetHeaderV1();
    jh->SetRparam(0.4); // setup parameters                                  
    if(radius>0)jh->SetRparam(radius);
    jh->SetAlgorithm(0); // kt from fastjet/JetDefinition.hh
    jh->SetDebug(-1);
    jh->SetPtMin(5);
    jh->SetJetEtaMax(0.9 - radius);
    jh->SetJetEtaMin(-0.9 + radius);
    jh->SetJetPhiMin( 0.);
    jh->SetJetPhiMax(2.*TMath::Pi());
    // Background
    jh->SetBGMode(kBackgroundMode);
    jh->SetBGAlgorithm(1);
    Double_t rBkg = 0.4;
    jh->SetRparamBkg(rBkg);
    jh->SetGhostEtaMax(0.9);
    jh->SetGhostArea(0.005);                            // area of a ghost
    jh->SetRapRange( -0.9+rBkg, 0.9-rBkg);              // rapidity range for subtracting background must be < ghostmaxrap-0.95*R
    jh->SetPhiRange(0 , 2*TMath::Pi());                 // phi range for subtracting background

    jetFinder = new AliFastJetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

    // SISCone Jet finder ------------------------------------------------
  case "SISCONE":
    AliSISConeJetHeader * jh = new AliSISConeJetHeader();

    //siscone parameters
    jh->SetConeRadius(0.4);                   // default cone radius
    if(radius>0)jh->SetConeRadius(radius);    // cone radius
    jh->SetDebug(-1);
    jh->SetJetEtaMax(0.9 - radius);
    jh->SetJetEtaMin(-0.9 + radius);
    jh->SetJetPhiMin( 0.);
    jh->SetJetPhiMax(2.*TMath::Pi());

    jh->SetOverlapThreshold(0.75);            // overlap parameter, between 0 and 1 excluded!! 0.75 value is advised
    jh->SetPtProtojetMin(0);                  // pt min of protojets
    jh->SetMinJetPt(5);                       // Ptmin of jets (GeV)

    // Background
    //do you want to subtract BG (0 = no, 1 = yes)
    jh->SetBGMode(kBackgroundMode);

    //for background
    Double_t rBkg = 0.4;
    jh->SetRapRange( -0.9+rBkg, 0.9-rBkg);              // rapidity range for subtracting background must be < ghostmaxrap-0.95*R
    jh->SetPhiRange(0 , 2*TMath::Pi());       // phi range for subtracting background
    
    //to determine jets area
    jh->SetBGAlgorithm(1);                    // algorithm for rho calculation : 0 = kT, 1 = Cambridge
    jh->SetGhostEtaMax(0.9);                  // eta max where ghosts are generated
    jh->SetGhostArea(0.005);                  // area of a ghost 
    jh->SetMeanGhostKt(1e-100);               // average transverse momentum of the ghosts.
    jh->SetAreaTypeNumber(4);                 // from 1 to 5 : 1 = active_area, 2 = active_area_explicit_ghosts, 3 = one_ghost_passive_area, 4 = passive_area, 5 = voronoi_area
    jetFinder = new AliSISConeJetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

    // UA1 Jet finder ------------------------------------------------
  case "UA1":
    AliUA1JetHeader *jh=new AliUA1JetHeader();
    jh->SetComment("UA1 jet code with default parameters");
    jh->SetRadius(0.4);
    if(radius>0)jh->SetRadius(radius);
    jh->SetDebug(-1);
    jh->SetEtSeed(2.); 
    jh->SetNAcceptJets(6);
    jh->SetLegoNbinPhi(432);
    jh->SetLegoNbinEta(274);
    jh->SetLegoEtaMin(-0.9);
    jh->SetLegoEtaMax(+0.9);
    jh->SetMinJetEt(5.);
    jh->SetJetEtaMax(0.9 - radius);
    jh->SetJetEtaMin(-0.9 + radius);
    jh->SetJetPhiMin(0.);
    jh->SetJetPhiMax(2.*TMath::Pi());
    jh->BackgMode(kBackgroundMode);

    jetFinder = new AliUA1JetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "UA1LO":
    AliUA1JetHeader *jh=new AliUA1JetHeader();
    jh->SetComment("UA1 jet code with Lo Pt settings parameters");
    jh->BackgMode(kBackgroundMode);
    jh->SetRadius(0.4);
    if(radius>0)jh->SetRadius(radius);
    jh->SetDebug(-1);
    jh->SetEtSeed(1.);
    jh->SetNAcceptJets(6);
    jh->SetLegoNbinPhi(432);
    jh->SetLegoNbinEta(274);
    jh->SetLegoEtaMin(-2);
    jh->SetLegoEtaMax(+2);
    jh->SetMinJetEt(1.);
    jh->SetJetEtaMax(1.5);
    jh->SetJetEtaMin(-1.5);

    jetFinder = new AliUA1JetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "UA1MC":
    AliUA1JetHeader *jh=new AliUA1JetHeader();
    jh->SetComment("UA1 jet code with default MC parameters");
    jh->BackgMode(kBackgroundMode);
    jh->SetRadius(1.0);
    if(radius>0)jh->SetRadius(radius);
    jh->SetDebug(-1); 
    jh->SetEtSeed(4.);
    jh->SetNAcceptJets(6);
    jh->SetLegoNbinPhi(432);
    jh->SetLegoNbinEta(274);
    jh->SetLegoEtaMin(-2);
    jh->SetLegoEtaMax(+2);
    jh->SetMinJetEt(5.);
    jh->SetJetEtaMax(1.5);
    jh->SetJetEtaMin(-1.5);
    jetFinder = new AliUA1JetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  default:
    ::Error("AddTaskJetsReader",Form("Wrong jet finder selected: %s\n",jf));
    return 0x0; 
    break;
  }

  return jetFinder;

}
