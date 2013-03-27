//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville

//by default this runs locally
//With the argument true this submits jobs to the grid
//As written this requires an xml script tag.xml in the ~/et directory on the grid to submit jobs
void runCaloEt(bool submit = false, // true or false 
	       const char *dataType="simPbPb", // "sim" or "real" etc.
	       const char *pluginRunMode="full", // "test" or "full" or "terminate"
	       const char *det = "EMCal",int production = 1, Bool_t withtender = kTRUE) // "PHOS" or "EMCAL" or EMCalDetail
{
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");

  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I. -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");



  if (!submit) { 
    cout << "local - no submitting" << endl;
  }
  else { 
    cout << "submitting to grid" << endl;
  }
   
  gROOT->ProcessLine(".L AliAnalysisEtCuts.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisHadEtCorrections.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtCommon.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtSelector.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtSelectorPhos.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtSelectorEmcal.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtTrackMatchCorrections.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtRecEffCorrection.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEt.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtMonteCarlo.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtMonteCarloPhos.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtMonteCarloEmcal.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtReconstructed.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtReconstructedPhos.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtReconstructedEmcal.cxx+g");  
  //gROOT->ProcessLine(".L AliAnalysisEtSelectionContainer.cxx+g");
  //gROOT->ProcessLine(".L AliAnalysisEtSelectionHandler.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisTaskTransverseEnergy.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEmEtMonteCarlo.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEmEtReconstructed.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisTaskTotEt.cxx+g");

  TString detStr(det);
  TString dataStr(dataType);
  if ( detStr.Contains("PHOS") ) {
    gSystem->CopyFile("calocorrections.PHOS.root","calocorrections.root",kTRUE);
    if ( dataStr.Contains("sim") ) {
      gSystem->CopyFile("ConfigEtMonteCarlo.PHOS.C","ConfigEtMonteCarlo.C",kTRUE);
    }
    else{
      gSystem->CopyFile("ConfigEtMonteCarlo.PHOS.data.C","ConfigEtMonteCarlo.C",kTRUE);
    }
  }
  else{
    gSystem->CopyFile("calocorrections.EMCAL.root","calocorrections.root",kTRUE);
    if ( dataStr.Contains("sim") ) {
      gSystem->CopyFile("ConfigEtMonteCarlo.EMCAL.C","ConfigEtMonteCarlo.C",kTRUE);
    }
    else{
      gSystem->CopyFile("ConfigEtMonteCarlo.EMCAL.data.C","ConfigEtMonteCarlo.C",kTRUE);
    }
  }


  char *kTreeName = "esdTree" ;
  TChain * chain   = new TChain(kTreeName,"myESDTree") ;
  
  if(submit){      
    gSystem->Load("libNetx") ; 
    gSystem->Load("libgapiUI");
    gSystem->Load("libRAliEn"); 
    TGrid::Connect("alien://") ;
  }
  
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TotEtManager");
  
  TString taskName = "TaskTotEt" + detStr;
  TString dataStrName(dataType);
  dataStrName.ReplaceAll("/",".");
  Bool_t isPb = kFALSE;
  if ( dataStr.Contains("PbPb") ) { isPb = kTRUE;}
  TString suffix = "";
  if(!withtender){
    suffix = "WithoutTender";
  }
  if(!isPb){
    suffix = "pp"+suffix;
  }
  TString outputName = "Et.ESD." + dataStrName + "." + detStr + ".root";
  TString outputDir = "totEt" + dataStr + detStr+suffix;


  cout << " taskName " << taskName
       << " outputName " << outputName 
       << " outputDir (alien) " << outputDir << endl;
  mgr->SetCommonFileName(outputName.Data());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("out1", TList::Class(), AliAnalysisManager::kOutputContainer, outputName);
  if (submit) {
    gROOT->LoadMacro("CreateAlienHandlerCaloEtSim.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandlerCaloEtSim(outputDir, outputName, pluginRunMode, production,detStr.Contains("PHOS"),!isPb,dataStr.Contains("real"));  
    if (!alienHandler) return;
    mgr->SetGridHandler(alienHandler);
  }

  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler* handler = new AliMCEventHandler;
  Bool_t isMc = kTRUE;
  if ( dataStr.Contains("sim") ) {
    cout << " MC " << endl;
    if ( dataStr.Contains("PbPb") ) { // a la: simPbPb/LHC10e18a
      cout << " PbPb " << endl;
      TString fileLocation = "/data/LHC10h8/137161/999/AliESDs.root";//"/home/dsilverm/data/E_T/" + dataStr + "/dir/AliESDs.root";
      cout << "fileLocation " << fileLocation.Data() << endl; 
//       chain->Add(fileLocation.Data()); // link to local test file
      chain->Add("/data/LHC10h8/137161/999/AliESDs.root");//Hijing Pb+Pb
       chain->Add("/data/LHC10h8/137161/111/AliESDs.root");//Hijing Pb+Pb
       chain->Add("/data/LHC10h8/137161/222/AliESDs.root");//Hijing Pb+Pb
    }
    else { // pp
      cout<<"adding pp simulation file"<<endl;
      chain->Add("/data/LHC10d15/1821/AliESDs.root");
      //chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data
      //chain->Add("/home/dsilverm/data/E_T/sim/LHC10d1/117222/100/AliESDs.root"); // link to local test file
    }
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }
  else { // real data
    cout<<"Hello there!  I am data."<<endl;
    isMc = kFALSE;
      chain->Add("/data/LHC10h/pass2_rev15/10000137366041.860/AliESDs.root");
      chain->Add("/data/LHC10h/pass2_rev15/10000137366041.870/AliESDs.root");
      chain->Add("/data/LHC10h/pass2_rev15/10000137366041.880/AliESDs.root");
      chain->Add("/data/LHC10h/pass2_rev15/10000137366041.890/AliESDs.root");
      chain->Add("/data/LHC10h/pass2_rev15/10000137366041.900/AliESDs.root");
//     chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data
    //chain->Add("/home/dsilverm/data/E_T/data/2010/LHC10b/000117222/ESDs/pass2/10000117222021.30/AliESDs.root"); // link to local test file
    cout << " not MC " << endl;
  }


  //if(!isMc && detStr.Contains("EMC")){
  if(detStr.Contains("EMC")){
    cout<<"You are running over EMCal data and using the tender supply"<<endl;
    gSystem->Load("libTENDER.so");
    gSystem->Load("libTENDERSupplies.so"); 
    gROOT->ProcessLine(".include $ALICE_ROOT/Tender/"); 
    gSystem->AddIncludePath("-I$ALICE_ROOT/ANALYSIS "); 
    //this macro is downloaded from the EMCal tender supply twiki 
    //hopefully it will be replaced by something checked in to aliroot
    //I have added the function from GetOCDBRecParam.C in Jiri's example to this so that we don't add gobs of macros to the code
    //I set the defaults to the golden run for PbPb because we are focusing on the golden run, however, this should be thought through!!
    gROOT->LoadMacro("AddTaskEMCALTenderForEtAnalysis.C");
    cout<<"WARNING: YOU ARE USING CALIBRATION FACTORS FROM PbPb RUN 137161!!"<<endl;
//  	// get reco params from grid OCDB
//    gROOT->LoadMacro("./GetOCDBRecParam.C");
//  	// run num, data type pp/PbPb, from grid
//Gets calibration factors from grid if jobs are to be submitted to the grid
//   	AliEMCALRecParam* pars = GetOCDBRecParam( 137161, "PbPb", submit);

    AliTender *tender = AddTaskEMCALTender( "EMCAL_COMPLETEV1", 0);
    //this also likely needs modification
//     tender->SelectCollisionCandidates( AliVEvent::kMB | AliVEvent::kEMCEGA | AliVEvent::kEMC1 | AliVEvent::kEMC7 );
//     if(submit){tender->SetDefaultCDBStorage("raw://");} //uncomment if you work on grid
//     else{tender->SetDefaultCDBStorage("local://$ALICE_ROOT/OCDB");} //uncomment if you work local

    if(submit){
      cout<<"Setting tender to run on grid"<<endl;
      tender->SetDefaultCDBStorage("raw://"); //uncomment if you work on grid
    }
    else{
      cout<<"Setting tender to run locally"<<endl;
      tender->SetDefaultCDBStorage("local://$ALICE_ROOT/OCDB"); //uncomment if you work local
    }
    // one can sellect what collision candidates to use
    // triggered sample only: L1 = AliVEvent::kEMCEGA, AliVEvent::kEMCEJE; L0 = AliVEvent::kEMC1, AliVEvent::kEMC7
    tender->SelectCollisionCandidates( AliVEvent::kAny );
    tender->SetDebugLevel(2);

    //AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("histosTrgContam", TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    //mgr->ConnectOutput(tender,1,coutput3);
    cout<<"Output container name "<<AliAnalysisManager::GetCommonFileName()<<endl;
  }

  if(isMc) cout<<"I am a MC"<<endl;
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  
  AliPhysicsSelectionTask *physicsSelectionTask = AddTaskPhysicsSelection(isMc);//isMC is true when processing monte carlo
  if(isPb){	 
    cout<<"Adding centrality selection task"<<endl;
    gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    //gROOT->ProcessLine(".L AliCentralitySelectionTask.cxx++g");
    AliCentralitySelectionTask *centTask = AddTaskCentrality();
    if(isMc){
     cout<<"Setting up centrality for MC"<<endl;
     centTask->SetMCInput();
   }
    else{
     cout<<"Setting up centrality for data"<<endl;
   }
  }



  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();





  AliAnalysisTaskTotEt *task1 = new AliAnalysisTaskTotEt(taskName);
  task1->SetMcData(isMc);//necessary to tell the task to basically accept all MC events.
  mgr->AddTask(task1);

  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,1,coutput1);


  
  mgr->SetDebugLevel(0);
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(submit){
    mgr->StartAnalysis("grid");
  }
  else{
    mgr->StartAnalysis("local",chain);
  }
  
  timer.Stop();
  timer.Print();
}
