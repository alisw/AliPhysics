//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville

//by default this runs locally
//With the argument true this submits jobs to the grid
//As written this requires an xml script tag.xml in the ~/et directory on the grid to submit jobs
void runCaloEt(bool submit = false, // true or false 
	       const char *dataType="simPbPb", // "sim" or "real" etc.
	       //const char *dataType="realPbPb", // "sim" or "real" etc.
	       const char *pluginRunMode="test", // "test" or "full" or "terminate"
	       const char *det = "EMCal",//EMCal",//"EMCal",//"EMCal",
	       int production = 1, Bool_t withtender = kTRUE, Int_t runnum = 0, Bool_t withNonlinearity = kTRUE, Bool_t withReclusterizing = kFALSE, Int_t trackmatchcuts=0, Bool_t is2011 = kFALSE, Bool_t jethad = kFALSE) // "PHOS" or "EMCAL" or EMCalDetail
{
  bool runCompiledVersion = kFALSE;
  class AliAnalysisEtCuts;
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");

  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I. -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS");
//   gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$CMAKE_INSTALL_PREFIX/include");
//   gSystem->AddIncludePath("-I. -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");

// gSystem->Load("libOADB");
// gSystem->Load("libCDB");
// gSystem->Load("libRAWDatabase");
gSystem->Load("libPHOSUtils");
gSystem->Load("libPHOSbase");
gSystem->Load("libPHOSpi0Calib");
gSystem->Load("libPHOSrec");
gSystem->Load("libPHOSshuttle");
gSystem->Load("libPHOSsim");
gSystem->Load("libPWGGAPHOSTasks");


//gSystem->Load("libTENDER.so");
//  gSystem->Load("libTENDERSupplies.so"); 
    gSystem->Load("libTender.so");
    gSystem->Load("libTenderSupplies.so"); 
    gSystem->Load("libPWGTools.so");
    gSystem->Load("libPWGEMCAL.so");
    gROOT->ProcessLine(".include $ALICE_ROOT/Tender/"); 
    //gSystem->AddIncludePath("-I$ALICE_ROOT/ANALYSIS "); 


  if (!submit) { 
    cout << "local - no submitting" << endl;
  }
  else { 
    cout << "submitting to grid" << endl;
  }
   
  if(runCompiledVersion){
    gSystem->Load("libPWGLFtotEt.so");
  }
  else{
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
  }
  TString detStr(det);
  TString dataStr(dataType);
  if ( detStr.Contains("PHOS") ) {
    if ( dataStr.Contains("sim") ) {      
      if(is2011){
	gSystem->CopyFile("calocorrections.2011.PHOS.sim.root","calocorrections.root",kTRUE); 
      }
      else{
	gSystem->CopyFile("calocorrections.PHOS.sim.root","calocorrections.root",kTRUE); 
      }
      gSystem->CopyFile("ConfigEtMonteCarlo.PHOS.C","ConfigEtMonteCarlo.C",kTRUE);
    }
    else{
      if(is2011){
	gSystem->CopyFile("calocorrections.2011.PHOS.root","calocorrections.root",kTRUE); 
      }
      else{
	gSystem->CopyFile("calocorrections.PHOS.root","calocorrections.root",kTRUE); 
      }
      gSystem->CopyFile("ConfigEtMonteCarlo.PHOS.data.C","ConfigEtMonteCarlo.C",kTRUE);
    }
  }
  else{//EMCal
    if ( dataStr.Contains("sim") ) {
      if(is2011){
	gSystem->CopyFile("ConfigEtMonteCarlo.EMCAL.2011.C","ConfigEtMonteCarlo.C",kTRUE);
	gSystem->CopyFile("calocorrections.2011.EMCAL.sim.root","calocorrections.root",kTRUE); 
      }
      else{
	gSystem->CopyFile("ConfigEtMonteCarlo.EMCAL.C","ConfigEtMonteCarlo.C",kTRUE);
	gSystem->CopyFile("calocorrections.EMCAL.sim.root","calocorrections.root",kTRUE); 
      }
    }
    else{
      if(is2011){
	gSystem->CopyFile("ConfigEtMonteCarlo.EMCAL.2011.data.C","ConfigEtMonteCarlo.C",kTRUE);
	gSystem->CopyFile("calocorrections.2011.EMCAL.root","calocorrections.root",kTRUE); 
      }
      else{
	gSystem->CopyFile("ConfigEtMonteCarlo.EMCAL.data.C","ConfigEtMonteCarlo.C",kTRUE);
	gSystem->CopyFile("calocorrections.EMCAL.root","calocorrections.root",kTRUE); 
      }
    }

  }
  

  if(is2011){
    gSystem->CopyFile("ConfigEtReconstructed.2011.C","ConfigEtReconstructed.C",kTRUE);
  }
  else{
      gSystem->CopyFile("ConfigEtReconstructed.2010.C","ConfigEtReconstructed.C",kTRUE);
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
  if(!withNonlinearity){
    outputDir +="NoNonlinearity";
  }
  if(withReclusterizing){
    outputDir +="WithReclusterizing";
  }
  if(trackmatchcuts!=0){
    outputDir +=Form("TrackMatchCut%i",trackmatchcuts);
  }
  if(jethad) outputDir+="WithJetHadronMethod";

  cout << " taskName " << taskName
       << " outputName " << outputName 
       << " outputDir (alien) " << outputDir << endl;
  mgr->SetCommonFileName(outputName.Data());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("out1", TList::Class(), AliAnalysisManager::kOutputContainer, outputName);
  if(!isPb){ cout<<"I am not PbPb!!"<<endl;}
  if (submit) {
    gROOT->LoadMacro("CreateAlienHandlerCaloEtSim.C");
    cout<<"Passing in production number "<<production<<endl;
    AliAnalysisGrid *alienHandler = CreateAlienHandlerCaloEtSim(outputDir, outputName, pluginRunMode, production,detStr.Contains("PHOS"),!isPb,dataStr.Contains("real"),runnum,runCompiledVersion);  
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
//      chain->Add("/data/tmp/3682/AliESDs.root");
//      chain->Add("/data/tmp/2782/AliESDs.root");
      // chain->Add("/data/LHC10h8/137161/999/AliESDs.root");//Hijing Pb+Pb
//       chain->Add("/data/LHC10h8/137161/111/AliESDs.root");//Hijing Pb+Pb
//       chain->Add("/data/LHC10h8/137161/222/AliESDs.root");//Hijing Pb+Pb
//      chain->Add("/data/LHC14a6/168464/605/AliESDs.root");
//chain->Add("/data/LHC11a10a_bis/139465/001/AliESDs.root");
//      chain->Add("/data/LHC14a6/168464/605/AliESDs.root");//HIJING with embedded signals
//      chain->Add("/data/LHC12d3/168464/201/AliESDs.root");//HIJING with embedded signals - works, full acceptance
      //chain->Add("/data/LHC14a6/168464/888/AliESDs.root");//HIJING with embedded signals
      //      //  chain->Add("/data/LHC11a10a_bis/139465/002/AliESDs.root");
      //chain->Add("/data/LHC11a10a_bis/139465/003/AliESDs.root");
      //chain->Add("/data/LHC11a10a_bis/139465/004/AliESDs.root");
//  chain->Add("/data/LHC11a10a_bis/139465/006/AliESDs.root");
//  chain->Add("/data/LHC11a10a_bis/139465/007/AliESDs.root");
//  chain->Add("/data/LHC11a10a_bis/139465/008/AliESDs.root");
//  chain->Add("/data/LHC11a10a_bis/139465/009/AliESDs.root");
 chain->Add("/data/LHC11a10a_bis/139465/010/AliESDs.root");
chain->Add("/data/LHC11a10a_bis/139465/011/AliESDs.root");
chain->Add("/data/LHC11a10a_bis/139465/012/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/013/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/014/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/015/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/016/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/017/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/018/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/019/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/020/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/021/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/022/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/023/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/024/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/025/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/026/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/027/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/028/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/029/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/030/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/031/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/032/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/033/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/034/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/035/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/036/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/037/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/038/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/039/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/040/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/041/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/042/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/043/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/044/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/045/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/046/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/047/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/048/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/049/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/050/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/051/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/052/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/053/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/054/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/055/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/056/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/057/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/058/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/059/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/060/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/061/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/062/AliESDs.root");

    }
    else { // pp
      cout<<"adding pp simulation file"<<endl;
      chain->Add("/data/LHC11b1b/999/AliESDs.root");
      //chain->Add("/data/LHC11b1a/999/AliESDs.root");
      //chain->Add("/data/LHC10d15/1821/AliESDs.root");
      //chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data
      //chain->Add("/home/dsilverm/data/E_T/sim/LHC10d1/117222/100/AliESDs.root"); // link to local test file
    }
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }
  else { // real data
    cout<<"Hello there!  I am data."<<endl;
    isMc = kFALSE;

    //      chain->Add("/data/tmp/10000139465010.600/AliESDs.root");

    //  chain->Add("/data/LHC11h/pass2/000168464/11000168464082.94/AliESDs.root");
    chain->Add("/data/LHC10h/pass2/139465/10000139465065.990/AliESDs.root");
//       chain->Add("/data/LHC10h/pass2/10000137366041.860/AliESDs.root");
//        chain->Add("/data/LHC10h/pass2/10000137366041.870/AliESDs.root");
//        chain->Add("/data/LHC10h/pass2/10000137366041.880/AliESDs.root");
//        chain->Add("/data/LHC10h/pass2/10000137366041.890/AliESDs.root");
//        chain->Add("/data/LHC10h/pass2/10000137366041.900/AliESDs.root");
//     chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data
    //chain->Add("/home/dsilverm/data/E_T/data/2010/LHC10b/000117222/ESDs/pass2/10000117222021.30/AliESDs.root"); // link to local test file
    cout << " not MC " << endl;
  }


  //if(!isMc && detStr.Contains("EMC")){
    if(detStr.Contains("EMC")){
  //if(0){
    cout<<"You are running over EMCal data and using the tender supply"<<endl;
    //this macro is downloaded from the EMCal tender supply twiki 
    //hopefully it will be replaced by something checked in to aliroot
    //I have added the function from GetOCDBRecParam.C in Jiri's example to this so that we don't add gobs of macros to the code
    //I set the defaults to the golden run for PbPb because we are focusing on the golden run, however, this should be thought through!!
    //AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance(geoname);

    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
    AliEmcalSetupTask *setupTask = AddTaskEmcalSetup();
    setupTask->SetGeoPath("$ALICE_PHYSICS/OADB/EMCAL");
    setupTask->SetOcdbPath(""); 

//     gROOT->LoadMacro("AddTaskEMCALTenderForEtAnalysis.C");
//     //cout<<"WARNING: YOU ARE USING CALIBRATION FACTORS FROM PbPb RUN 137161!!"<<endl;
// //  	// get reco params from grid OCDB
// //    gROOT->LoadMacro("./GetOCDBRecParam.C");
// //  	// run num, data type pp/PbPb, from grid
// //Gets calibration factors from grid if jobs are to be submitted to the grid
// //   	AliEMCALRecParam* pars = GetOCDBRecParam( 137161, "PbPb", submit);
// //EMCAL_FIRSTYEARV1 F-
//     //AliTender *tender = AddTaskEMCALTender( "EMCAL_COMPLETEV1", 0,withNonlinearity,withReclusterizing,trackmatchcuts);
// AliTender *tender = AddTaskEMCALTender( "EMCAL_FIRSTYEARV1", 0,withNonlinearity,withReclusterizing,trackmatchcuts);
//     //this also likely needs modification
// //     tender->SelectCollisionCandidates( AliVEvent::kMB | AliVEvent::kEMCEGA | AliVEvent::kEMC1 | AliVEvent::kEMC7 );
// //     if(submit){tender->SetDefaultCDBStorage("raw://");} //uncomment if you work on grid
// //     else{tender->SetDefaultCDBStorage("local://$ALICE_ROOT/OCDB");} //uncomment if you work local

//     if(submit){
//       cout<<"Setting tender to run on grid"<<endl;
//       tender->SetDefaultCDBStorage("raw://"); //uncomment if you work on grid
//     }
//     else{
//       cout<<"Setting tender to run locally"<<endl;
//       tender->SetDefaultCDBStorage("local://$ALICE_ROOT/OCDB"); //uncomment if you work local
//     }


    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C");//tendertasks
    TString runPeriod = "LHC10h";
  Bool_t distBC         = kTRUE;   //distance to bad channel
  Bool_t recalibClus    = kTRUE;   //recalibrate cluster energy
  Bool_t recalcClusPos  = kTRUE;   //recalculate cluster position
  Bool_t nonLinearCorr  = kTRUE;   //apply non-linearity
  Bool_t remExotic      = kTRUE;   //remove exotic cells
  Bool_t fidRegion      = kTRUE;  //apply fiducial cuts -->  different from defaults
  Bool_t calibEnergy    = kTRUE;   //calibrate energy
  Bool_t calibTime      = kTRUE;   //calibrate timing
  Bool_t remBC          = kTRUE;   //remove bad channels
  UInt_t nonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrected;
  Bool_t reclusterize   = kFALSE;   //reclusterize --> different from defaults
  Float_t seedthresh    = 0.100;   //seed threshold
  Float_t cellthresh    = 0.050;   //cell threshold
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2;
  Bool_t trackMatch     = kTRUE;   //track matching
  Bool_t updateCellOnly = kFALSE;  //only change if you run your own clusterizer task
  Float_t timeMin       = 100e-9;  //minimum time of physical signal in a cell/digit (s)
  Float_t timeMax       = 900e-9;  //maximum time of physical signal in a cell/digit (s)
  Float_t timeCut       = 900e-9;  //maximum time difference between the digits inside EMC cluster (s)
    const char *pass      = 0 ;       //string defining pass (use none if figured out from path)
    //AliAnalysisTaskSE *tender = AddTaskEMCALTender();
    AliAnalysisTaskSE *tender = AddTaskEMCALTender(distBC, recalibClus, recalcClusPos, nonLinearCorr, remExotic, 
						   fidRegion, calibEnergy, calibTime, remBC, nonLinFunct, reclusterize, seedthresh, 
						   cellthresh, clusterizer, trackMatch, updateCellOnly, timeMin, timeMax, timeCut);
    


    // one can sellect what collision candidates to use
    // triggered sample only: L1 = AliVEvent::kEMCEGA, AliVEvent::kEMCEJE; L0 = AliVEvent::kEMC1, AliVEvent::kEMC7
    tender->SelectCollisionCandidates( AliVEvent::kAny );
    //tender->SetDebugLevel(2);

    //AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("histosTrgContam", TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    //mgr->ConnectOutput(tender,1,coutput3);
    cout<<"Output container name "<<AliAnalysisManager::GetCommonFileName()<<endl;
    }
    else{
      if(withtender){
	cout<<"You are running over PHOS data and using the tender supply"<<endl;
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C");//isMc
	TString mcname = "LHC13d2";
	if(is2011) mcname = "LHC13d2";
	//AliAnalysisTaskSE *tender = AddAODPHOSTender("PHOSTenderTask","PHOStender",mcname.Data(),2,isMc);//last argument is pass 2, not set to assume MC
	
	
	// one can sellect what collision candidates to use
	// triggered sample only: L1 = AliVEvent::kEMCEGA, AliVEvent::kEMCEJE; L0 = AliVEvent::kEMC1, AliVEvent::kEMC7
	//tender->SelectCollisionCandidates( AliVEvent::kAny );
      }
    }

  if(isMc) cout<<"I am a MC"<<endl;
  gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  
  AliPhysicsSelectionTask *physicsSelectionTask = AddTaskPhysicsSelection(isMc);//isMC is true when processing monte carlo
  if(isPb){	 
    cout<<"Adding centrality selection task"<<endl;
    gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
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

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  //AliAnalysisTask *AddTaskPIDResponse(Bool_t isMC=kFALSE, Bool_t autoMCesd=kTRUE,
//                                     Bool_t tuneOnData=kFALSE, Int_t recoPass=2,
//                                     Bool_t cachePID=kFALSE, TString detResponse="",
//                                     Bool_t useTPCEtaCorrection = kFALSE);
//  AliAnalysisTask *AddTaskPIDResponse(Bool_t isMC=kFALSE, Bool_t autoMCesd=kTRUE,
//                                  Bool_t tuneOnData=kFALSE, Int_t recoPass=2,
//                                  Bool_t cachePID=kFALSE, TString detResponse="",
//                                  Bool_t useTPCEtaCorrection = kTRUE,
//                                  Bool_t useTPCMultiplicityCorrection = kFALSE
//                                  Int_t  userDataRecoPass = -1)

  AliAnalysisTask *taskPID;
  if(submit){
    taskPID=AddTaskPIDResponse(isMc);//,kTRUE,kTRUE,2,kFALSE,"",kTRUE,kFALSE,2);
  }
  else{
    cout<<"Not submitting so forcing pass number locally so it doesn't crash"<<endl;
    taskPID=AddTaskPIDResponse(isMc,kTRUE,kTRUE,2,kFALSE,"",kTRUE,kFALSE,2);
  }
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  //AddTaskPIDqa();

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();





  AliAnalysisTaskTotEt *task1 = new AliAnalysisTaskTotEt(taskName);
  task1->SetMcData(isMc);//necessary to tell the task to basically accept all MC events.
  task1->SelectCollisionCandidates(AliVEvent::kMB ) ;
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
