AliJetReader *CreateJetReader(Char_t *jr,UInt_t filterMask); // Common config
AliJetFinder *CreateJetFinder(Char_t *jf,Float_t radius = -1);

AliAnalysisTaskJets *AddTaskJets(Char_t *jr, Char_t *jf,Float_t radius = -1,UInt_t filterMask = 0,Float_t ptTrackMin = 0.15,Int_t iBack = 0); // for the new AF
Int_t AddTaskJetsDelta(char *nonStdFile = "",UInt_t filterMask = 0,Bool_t kUseAODMC = kTRUE,UInt_t runFlag = 1|4|32|128|256);     
AliAnalysisTaskJets *AddTaskJets(UInt_t filterMask = 0);

Float_t kPtTrackMin = 0.15;
Int_t kBackgroundMode = 0;


AliAnalysisTaskJets *AddTaskJets(UInt_t filterMask ){
  // fills the standard "jets" branch in the AOD
  // need the ESDFilter to run before, to access the AODtracks
  // Tracks selected by the first Filter (1<<0)
  // needs to be adapted for different cuts
  
  // UA1 as standard chosen, since it is the most robust and simple JF
  // R = 0.4 suffficient to provide accurate jet axis for correlation studies
  // energy resolution suffers a little
  // Acceptance of jets not limited by the Jet Finder but should be done
  // by user to abs(eta) < 0.5 

  return AddTaskJets("AOD","UA1",0.4,filterMask);

}



Int_t AddTaskJetsDelta(char *nonStdFile,UInt_t filterMask,Bool_t kUseAODMC,UInt_t runFlag){

  // Adds a whole set of jet finders  all to be written
  // to a delta AOD
  
  // this can in principle be done also with on the fly 
  // if we have an ESD input jet fidner task does automatically fetch the ouput aod

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJetsDelta", "No analysis manager to connect to.");
      return 0;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskJetsDelta", "This task requires an input event handler");
      return 0;
   }

  AliAODHandler *aodh = (AliAODHandler*)mgr->GetOutputEventHandler();
  if (!aodh) {
    ::Error("AddTaskJetsDelta", "This task needs an output event handler");
    return 0;
  }   

  TString outFile(nonStdFile);

  AliAnalysisTaskJets *jetana = 0;
  Int_t iCount = 0;

  // Jet Fidners Selected by run flag first bit 2^0 second by 2^1 etc
  const char *cJF[14]        = {"UA1","UA1","UA1", "CDF", "CDF","DA","DA","SISCONE","SISCONE","FASTJET","FASTJET","FASTKT","FASTKT","UA1LO"};
  const Float_t radius[14]   = {  0.4,  0.7,  1.0,  0.4,    0.7, 0.4, 0.7,      0.4,     0.7,      0.4,       0.7,     0.4,     0.7,    0.4};
  UInt_t  flag[14]           = {    6,    7,    7,    7,      7,   7,   7,        7,       7,        7,         7,       7,       7,      7};
  // flag[5] = 0; // set siscone to 0 for proof mode...
  // flag first bit AOD, second bit AODMC2 third bit AODMC2 third (8) bit AOODMC2b (limited acceptance)
  // i.e. 7 all, 6 only MC2 and MC
  // this stay at three
  const char *cReader[4] = {"AOD","AODMC","AODMC2","AODMC2b"};  

  

  for(int i = 0; i< 14;i++){
    if(!(runFlag&(1<<i)))continue;
    if(!kUseAODMC)flag[i]&=1; // switch OFF MC if we do not have it
    for(int ib = 0;ib<4;ib++){      
      if(flag[i]&(1<<ib)){
	jetana = AddTaskJets(cReader[ib],cJF[i],radius[i],filterMask);
	if(jetana){
	  char *cRadius = "";
	  if(radius[i]>0)cRadius = Form("%02d",(int)((radius[i]+0.01)*10.)); // add an offset beacuse of precision
	  //	  jetana->SetNonStdBranch(Form("jets%s_%s%s",cReader[ib],cJF[i],cRadius)); // done in addtask jets
	  if(outFile.Length()>0)jetana->SetNonStdOutputFile(outFile.Data());
	  iCount++;
	}
      }
    }
  }
    
  Printf("Added %d JetFinders",iCount);
}





AliAnalysisTaskJets *AddTaskJets(Char_t *jr, Char_t *jf, Float_t radius,UInt_t filterMask,float ptTrackMin,int iBack)
{
  // Creates a jet finder task, configures it and adds it to the analysis manager.

  kPtTrackMin = ptTrackMin;
  kBackgroundMode = iBack;
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
   AliJetReader *er = CreateJetReader(jr,filterMask);
   // Define jet header and jet finder
   AliJetFinder *jetFinder = CreateJetFinder(jf,radius);

   if (jetFinder){
       if (er) jetFinder->SetJetReader(er);
   }



   TString cAdd = "";
   cAdd += Form("%02d",(int)((radius+0.01)*10.));
   cAdd += Form("_B%d",(int)((kBackgroundMode)));
   cAdd += Form("_Filter%05d",filterMask);
   cAdd += Form("_Cut%05d",(int)((1000.*kPtTrackMin)));
   Printf("%s",cAdd.Data());
   jetana = new AliAnalysisTaskJets(Form("JetAnalysis%s_%s%s",jr,jf,cAdd.Data()));


   TString c_jr(jr);
   c_jr.ToLower();
   TString c_jf(jf);
   c_jf.ToLower();

   /*
   if(c_jf.CompareTo("ua1")==0 && TMath::Abs(radius-0.4) < 0.01 && c_jr.CompareTo("aod") == 0){
     // do nothing, this is the standard jet finder R = 0.4, UA1 on AOD
   }
   else{

     
 }
   */
   TString bName =  Form("jets%s_%s%s",jr,jf,cAdd.Data());
   jetana->SetNonStdBranch(bName.Data());
   Printf("Set jet branchname \"%s\"",bName.Data());


   AliAnalysisDataContainer *cout_jet = mgr->CreateContainer(
							     Form("jethist_%s_%s%s",
								  c_jr.Data(),
								  c_jf.Data(),
								  cAdd.Data()), 
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer, 
							     Form("%s:PWGJE_jethist_%s_%s%s",AliAnalysisManager::GetCommonFileName(),
								  c_jr.Data(),
								  c_jf.Data(),
								  cAdd.Data()));

   // Connect jet finder to task.
   jetana->SetJetFinder(jetFinder);
   jetana->SetConfigFile("");
   jetana->SetDebugLevel(2);
   if(TMath::Abs((radius-0.4))< 0.02&&c_jf.Contains("fastjet")){
     //     jetana->SetFilterPt(20.);
   }


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
    if(radius>0)jh->SetRadius(radius);    
    jh->SetAODwrite(kTRUE);
    jh->SetAODtracksWrite(kTRUE);
    //    jh->SetDebugCDF(1);
    jetFinder = new AliCdfJetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "DA":
    AliDAJetHeader *jh=new AliDAJetHeader();
    jh->SetComment("DA jet code with default parameters");
    jh->SelectJets(kTRUE);
//    jh->SetNeff(200);
//    jh->SetEtaEff(2.2);
    if(radius>0)jh->SetRadius(radius);
    jh->SetEtMin(5.);
    jetFinder = new AliDAJetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "FASTJET":
    // DEFAULT is ANTI KT
    AliFastJetHeaderV1 *jh = new AliFastJetHeaderV1();
    jh->SetRparam(0.4); // setup parameters                                  
    Double_t rBkg = 0.6;
    if(radius>0)jh->SetRparam(radius);
    jh->SetAlgorithm(2); // antikt from fastjet/JetDefinition.hh
    /*
      $FASTJET/include/fastjet/JetDefinition.hh
      enum JetAlgorithm {kt_algorithm, cambridge_algorithm,
      antikt_algorithm, genkt_algorithm,
      ee_kt_algorithm, ee_genkt_algorithm, ...};

     */
    jetFinder = new AliFastJetFinder();
    jh->SetPtMin(1);
    jh->SetRparamBkg(rBkg);
    jh->SetGhostEtaMax(0.9+rBkg);
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "FASTKT":
    AliFastJetHeaderV1 *jh = new AliFastJetHeaderV1();
    jh->SetRparam(0.4); // setup parameters                                  
    Double_t rBkg = 0.6;
    if(radius>0)jh->SetRparam(radius);
    jh->SetAlgorithm(0); // kt from fastjet/JetDefinition.hh
    jh->SetPtMin(1);
    jh->SetRparamBkg(rBkg);
    jh->SetGhostEtaMax(0.9+rBkg);
    jetFinder = new AliFastJetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "SISCONE":
    AliSISConeJetHeader * jh = new AliSISConeJetHeader();

    jh->SetJetEtaMax(1.5);
    jh->SetJetEtaMin(-1.5);

    //siscone parameters
    jh->SetConeRadius(0.4);                   // default cone radius
    if(radius>0)jh->SetConeRadius(radius);   // cone radius

    jh->SetOverlapThreshold(0.75);            // overlap parameter, between 0 and 1 excluded!! 0.75 value is advised
    jh->SetPtProtojetMin(0);                  // pt min of protojets
    jh->SetMinJetPt(5);                      // Ptmin of jets (GeV)

    //do you want to subtract BG (0 = no, 1 = yes)
    jh->SetBGMode(0); // if 1 set also the radius for the background determination..

    //for background
    jh->SetRapRange( -0.9, 0.9);              // rapidity range for subtracting background must be < ghostmaxrap-0.95*R
    jh->SetPhiRange(0 , 2*TMath::Pi());       // phi range for subtracting background
    
    //to determine jets area
    jh->SetBGAlgorithm(1);                    // algorithm for rho calculation : 0 = kT, 1 = Cambridge
    jh->SetGhostEtaMax(4);                    // eta max where ghosts are generated 
    jh->SetGhostArea(0.05);                   // area of a ghost 
    jh->SetMeanGhostKt(1e-100);               // average transverse momentum of the ghosts.
    jh->SetAreaTypeNumber(4);                 // from 1 to 5 : 1 = active_area, 2 = active_area_explicit_ghosts, 3 = one_ghost_passive_area, 4 = passive_area, 5 = voronoi_area
    jetFinder = new AliSISConeJetFinder();
    if (jh) jetFinder->SetJetHeader(jh);
    break;

  case "UA1":
    AliUA1JetHeaderV1 *jh=new AliUA1JetHeaderV1();
    jh->SetComment("UA1 jet code with default parameters");
    jh->BackgMode(0);
    jh->SetRadius(0.4);
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
    jh->BackgMode(kBackgroundMode);

    jetFinder = new AliUA1JetFinderV1();
    if (jh) jetFinder->SetJetHeader(jh);
    break;
  case "UA1LO":
    AliUA1JetHeaderV1 *jh=new AliUA1JetHeaderV1();
    jh->SetComment("UA1 jet code with Lo Pt settings parameters");
    jh->BackgMode(0);
    jh->SetRadius(0.4);
    if(radius>0)jh->SetRadius(radius);
    jh->SetEtSeed(1.);
    jh->SetNAcceptJets(6);
    jh->SetLegoNbinPhi(432);
    jh->SetLegoNbinEta(274);
    jh->SetLegoEtaMin(-2);
    jh->SetLegoEtaMax(+2);
    jh->SetMinJetEt(1.);
    jh->SetJetEtaMax(1.5);
    jh->SetJetEtaMin(-1.5);
    jh->BackgMode(kBackgroundMode);

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
    jh->BackgMode(kBackgroundMode);
    jetFinder = new AliUA1JetFinderV1();
    if (jh) jetFinder->SetJetHeader(jh);
    break;
  default:
    Printf("\n >>>>>>> AddTaskJets Error Wrong jet finder selected\n");
    break;
  }

  return jetFinder;

}

AliJetReader *CreateJetReader(Char_t *jr,UInt_t filterMask){
  AliJetReader *er = 0;
  Float_t ptCut  = kPtTrackMin ; // cut on track p_T
  switch (jr) {
  case "MC":
    AliJetKineReaderHeader *jrh = new AliJetKineReaderHeader();
    jrh->SetComment("MC full Kinematics");
    jrh->SetFastSimTPC(kFALSE);
    jrh->SetFastSimEMCAL(kFALSE);
    jrh->SetPtCut(ptCut);
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
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0 .9                                                                             
    // Define reader and set its header                                     
    er = new AliJetKineReader();
    er->SetReaderHeader(jrh);
    break;
  case "ESD":
    AliJetESDReaderHeader *jrh = new AliJetESDReaderHeader();
    jrh->SetComment("Testing");
    jrh->SetFirstEvent(0.);
    jrh->SetLastEvent(1000);
    jrh->SetPtCut(ptCut);
    jrh->SetReadSignalOnly(kFALSE);
    // Define reader and set its header                                     
    er = new AliJetESDReader();
    er->SetReaderHeader(jrh);
    break;

  case "AOD":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD Reader");
    jrh->SetPtCut(ptCut); // set low p_T cut of to 150 MeV
    jrh->SetTestFilterMask(32); // Change this one for a different set of cuts
    if(filterMask>0)jrh->SetTestFilterMask(filterMask); 
    // Define reader and set its header
    er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODMC":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD MC Reader");
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0.9
    jrh->SetReadAODMC(1);// 1 all primary MC , 2 all primary charged
    // Define reader and set its header
    er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODMC2":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD MC Reader");
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0.9
    jrh->SetReadAODMC(2);// 1 all primary MC , 2 all primary charged
    // Define reader and set its header
    er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODMC2b":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD MC Reader");
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-0.9,0.9); // to take all MC particles default is 0.9
    jrh->SetReadAODMC(2);// 1 all primary MC , 2 all primary charged
    // Define reader and set its header
    er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODextra":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD Reader with extra branch");
    jrh->SetPtCut(ptCut); // set low p_T cut of to ptCut
    jrh->SetTestFilterMask(32); // Change this one for a different set of cuts
    jrh->SetReadBranches(AliJetAODReaderHeader::kReadStdBranch | AliJetAODReaderHeader::kReadNonStdBranch);
    jrh->SetNonStdBranch("aodExtraTracks");
    if(filterMask>0)jrh->SetTestFilterMask(filterMask); 
    // Define reader and set its header
    er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODextraonly":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD Reader with extra branch");
    jrh->SetPtCut(ptCut); // set low p_T cut of to 150 MeV
    jrh->SetTestFilterMask(32); // Change this one for a different set of cuts
    jrh->SetReadBranches(AliJetAODReaderHeader::kReadNonStdBranch);
    jrh->SetNonStdBranch("aodExtraTracks");
    if(filterMask>0)jrh->SetTestFilterMask(filterMask); 
    // Define reader and set its header
    er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODMCextra":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD MC Reader with extra branch");
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0.9
    jrh->SetReadAODMC(1);// 1 all primary MC , 2 all primary charged
    jrh->SetReadBranches(AliJetAODReaderHeader::kReadStdBranch | AliJetAODReaderHeader::kReadNonStdBranch);
    jrh->SetNonStdBranch("aodExtraMCparticles");
    // Define reader and set its header
    er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
    break;
  case "AODMC2extra":
    AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("AOD MC Reader with extra branch");
    jrh->SetPtCut(ptCut);
    jrh->SetFiducialEta(-2.1,2.1); // to take all MC particles default is 0.9
    jrh->SetReadAODMC(2);// 1 all primary MC , 2 all primary charged
    jrh->SetReadBranches(AliJetAODReaderHeader::kReadStdBranch | AliJetAODReaderHeader::kReadNonStdBranch);
    jrh->SetNonStdBranch("aodExtraMCparticles");
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
