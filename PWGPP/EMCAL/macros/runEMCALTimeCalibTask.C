/// \file runEMCALTimeCalibTask.C
/// \ingroup EMCALPerformanceMacros
/// \brief macro to run Time Calibration Task
///
/// For Run1 data run in two iterations. 
/// 1. First iteration over statistics gives calibration constants.
/// Remember to apply energy cuts. We want to take into account only GOOD cells.
/// 2. Run method ProduceCalibConsts(inputFile,outputFile,fFALSE) to obtain reference 
/// file with constants. 
/// 3. Second iteration with applied consatants plus extra timing cut gives corrected calibration constants.
/// 4. Extract final calib constants ProduceCalibConsts(inputFile,outputFile,kTRUE).
///
/// For Run2 data run in three iterations. 
/// 1. The same as 1 point for Run1 but run by run.
/// 2. The same as 2 point for Run1 but run by run.
/// 3. Extract L1 phase via ProduceOffsetForSMsV2(runNumber,inputFile,outputFile)
/// to obtain phases per run; outputFile is the same.
/// 4. Second iteration with info from L1 phase with loose time criteria
/// 5. Run method ProduceCalibConsts(inputFile,outputFile,kTRUE) to obtain reference 
/// file with constants for whole collection.
/// 6. Third iteration with applied consatants plus L1phases, plus extra timing cut gives corrected calibration constants.
///
/// \param type: Int_t, local or grid processing
/// \param isESD: Bool_t, flag to process ESD or AOD
/// \param isPhysicsSelection: Bool_t, flag to use Physics Selection
/// \param isCentralitySelection: Bool_t, flag to use Centrality Selection
/// \author Marie Germain <marie.germain@subatech.in2p3.fr>, SUBATECH
/// \author Adam Matyja <adam.tomasz.matyja@ifj.edu.pl>, INP PAN Cracow
/// \date Jun 3, 2015

void runEMCALTimeCalibTask(Int_t type=0, Bool_t isESD=kTRUE, Bool_t isPhysicsSelection=kFALSE, Bool_t isCentralitySelection=kFALSE)
{
  //type = 0: local 
  //type = 1: grid

  TStopwatch timer;
  timer.Start();
  printf("*** the task is: AliAnalysisTaskEMCALTimeCalib ***\n");
  // Grid connection
  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
  //printf("*** Connected to AliEn? in  ***\n");
  cout<< "Connected to AliEn" <<endl;
  timer.Print();
  
  //libraries
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  //gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  
  // load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSISalice.so");

  gSystem->Load("libEMCALUtils.so");
  gSystem->Load("libPWGPPEMCAL.so");
  
  gSystem->AddIncludePath("-I$ALICE_ROOT");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS");
  gSystem->AddIncludePath("-I./");     

  //ANALYSIS PART
  const char *collectionfile;
  if(isESD==kTRUE) collectionfile = "esd_156889.xml";
  else collectionfile= "aod_156889.xml";
  TChain* chain = NULL;
  if(type==0){//local files
    if(isESD==kTRUE) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      chain = CreateESDChain("files.txt", 1);
    } else {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      chain = CreateAODChain("filesAOD.txt", 1);
    }
  }
  else if(type==1){//grid files
    if(isESD==kTRUE) chain = CreateChainFromCollection(collectionfile, "esdTree");
    else chain = CreateChainFromCollection(collectionfile, "aodTree");
  }

  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$(ALICE_ROOT)/include
  // or in each macro
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");



  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MyEmcalAnalysis");
  AliInputEventHandler *inputHandler=NULL;

  //  AliESDInputHandler* esdH = new AliESDInputHandler();
  //  mgr->SetInputEventHandler(esdH);

  if(isESD==kTRUE) inputHandler = new AliESDInputHandler();
  else inputHandler = new AliAODInputHandler();

  mgr->SetInputEventHandler(inputHandler);

  AliPhysicsSelectionTask* physselTask=NULL;
  if(isPhysicsSelection){
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskPhysicsSelection.C");
    //    AliPhysicsSelectionTask* physselTask = AddTaskPhysicsSelection();
    physselTask = AddTaskPhysicsSelection();
    // physselTask->AddCollisionTriggerClass("+CINT1WU-B-NOPF-ALL");
  }
  AliCentralitySelectionTask* centralityTask=NULL;
  if(isCentralitySelection){
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    //    AliCentralitySelectionTask* centralityTask = AddTaskCentrality(); // Create task
    centralityTask = AddTaskCentrality(); // Create task
  }

  //AliCalorimeterUtils *cu = new AliCalorimeterUtils();
  //   cu->SwitchOnBadChannelsRemoval(); 

  // Create task
  //UInt_t kTriggerInt = AliVEvent::kAnyINT;
  //UInt_t kTriggerEMC   = AliVEvent::kEMC8 || AliVEvent::kEMC7 || AliVEvent::kEMCEJE || AliVEvent::kEMCEGA;
  
  AliAnalysisTaskEMCALTimeCalib *taskmbemcal = new AliAnalysisTaskEMCALTimeCalib("TimeCalibTask");
  //  AliEMCALRecoUtils *reco = taskmbemcal->GetEMCALRecoUtils();
  //  taskmbemcal->PrintInfo();
  //  reco->SwitchOnRejectExoticCluster();

  taskmbemcal->SelectCollisionCandidates(AliVEvent::kEMC1|AliVEvent::kEMC7|AliVEvent::kEMC8|AliVEvent::kEMCEJE|AliVEvent::kEMCEGA);
  //taskmbemcal->SetGeometryName("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  //taskmbemcal->SwithOffFillHeavyHisto();

  // to load reference histograms run-by-run in pass2
  //taskmbemcal->SetReferenceRunByRunFileName("ReferenceSM_LHC15n.root");
  //taskmbemcal->LoadReferenceRunByRunHistos();
  ///taskmbemcal->SwitchOffBadReco(); 
  //taskmbemcal->SetMinTime(200);
  //taskmbemcal->SetMaxTime(800);
  //taskmbemcal->SetPassTimeHisto(800,400.,800.);

  //to apply L1 shift correction for badly reconstructed runs in muon_calo_pass1 LHC15f-m in run2
  //additionally do in pass2:
  //taskmbemcal->SwitchOnBadReco();
  //taskmbemcal->SetMinTime(-200);
  //taskmbemcal->SetMaxTime(200);    
  //taskmbemcal->SetPassTimeHisto(800,-200.,200.);


  // to load reference calibration constants in pass3
  //taskmbemcal->SetReferenceFileName("Reference_LHC15i.root");
  //taskmbemcal->LoadReferenceHistos();
  //taskmbemcal->SetPassTimeHisto(1400,-350.,350.);

  taskmbemcal->SetBadChannelMapSource(0);
  if(taskmbemcal->GetBadChannelMapSource()==2) taskmbemcal->SetBadChannelFileName("badMap.root");


  //taskmbemcal->SelectCollisionCandidates(AliVEvent::kAnyINT);
  //taskmbemcal->SetDebugLevel(10);

  // Add task(s)
  // taskmbemcal->PrintInfo();
  if(isPhysicsSelection) mgr->AddTask(physselTask);
  if(isCentralitySelection) mgr->AddTask(centralityTask);
  mgr->AddTask(taskmbemcal);

  TString outputFileName;
  if(type==0)outputFileName="timeOutput.root";
  else if(type==1)outputFileName="timeResults.root";

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutputList = mgr->CreateContainer("chistolist", TList::Class(),
     AliAnalysisManager::kOutputContainer, outputFileName.Data());


  // Connect input/output
  mgr->ConnectInput(taskmbemcal, 0, cinput);

  // No need to connect to a common AOD output container if the task does not
  // fill AOD info.
  //  mgr->ConnectOutput(task, 0, coutput0);
  mgr->ConnectOutput(taskmbemcal, 1, coutputList);

  // Enable debug printouts
  if(type==0) mgr->SetDebugLevel(12);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain);
  

  timer.Stop();
  timer.Print();
  
}

//______________________________________________________________________________
TChain *CreateChainFromCollection(const char* xmlfile, const char *treeName="esdTree")
{
// Create a chain from an alien collection.
   TAlienCollection * myCollection  = TAlienCollection::Open(xmlfile);

   if (!myCollection) {
      ::Error("CreateChainSingle", "Cannot create an AliEn collection from %s", xmlfile) ;
      return NULL ;
  }

  TChain* chain = new TChain(treeName);
  myCollection->Reset() ;
  while ( myCollection->Next() ) chain->Add(myCollection->GetTURL("")) ;
  chain->ls();
  return chain;
}

// End of macro
