/// \file runEMCALTimeCalibTask.C
/// \brief macro to run Time Calibration Task
///
/// Run in two steps. 
/// 1. First iteration over statistics gives calibration constants.
/// Remember to apply energy cuts. We want to take into account only GOOD cells.
/// 2. Run method ProduceCalibConsts(inputFile,outputFile) to obtain reference 
/// file with constants. 
/// 3. Second iteration with applied conatants give results.
///
/// \author Marie Germain <marie.germain@subatech.in2p3.fr>, SUBATECH
/// \author Adam Matyja <adam.tomasz.matyja@ifj.edu.pl>, INP PAN Cracow
/// \date Jun 3, 2015

void runTimeCalibTask(Int_t type=0)
{
  //type = 0: local 
  //type = 1: grid

  Bool_t isPhysicsSelection=kFALSE;// <<-----FIXME <<-----uncomment it in real data
  Bool_t isCentralitySelection=kFALSE;

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
  //gSystem->Load("libAOD.so");
  
  // load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSISalice.so");

  gSystem->Load("libEMCALUtils.so");
  gSystem->Load("libPWGPPEMCAL.so");
  
  gSystem->AddIncludePath("-I$ALICE_ROOT");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS");
  gSystem->AddIncludePath("-I./");     

  // 	// test with lib from PWGGA train
  



  //ANALYSIS PART
  const char *collectionfile = "wn.xml";
  TChain* chain = NULL;
  if(type==0){//local files
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
    chain = CreateESDChain("files.txt", 1);
  }
  else if(type==1){//grid files
    chain = CreateChainFromCollection(collectionfile, "esdTree");
  }

  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$(ALICE_ROOT)/include
  // or in each macro
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MyEmcalAnalysis");
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);
  
  if(isPhysicsSelection){
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physselTask = AddTaskPhysicsSelection();
    // physselTask->AddCollisionTriggerClass("+CINT1WU-B-NOPF-ALL");
  }
  if(isCentralitySelection){
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask* centralityTask = AddTaskCentrality(); // Create task
  }

  AliCalorimeterUtils *cu = new AliCalorimeterUtils();
  //   cu->SwitchOnBadChannelsRemoval(); 

  // Create task
  UInt_t kTriggerInt = AliVEvent::kAnyINT;
  UInt_t kTriggerEMC   = AliVEvent::kEMC8 || AliVEvent::kEMC7 || AliVEvent::kEMCEJE || AliVEvent::kEMCEGA;
  
  AliAnalysisTaskEMCALTimeCalib *taskmbemcal = new AliAnalysisTaskEMCALTimeCalib("TimeCalibTask");
  //  AliEMCALRecoUtils *reco = taskmbemcal->GetEMCALRecoUtils();
  //  taskmbemcal->PrintInfo();
  //  reco->SwitchOnRejectExoticCluster();

  taskmbemcal->SelectCollisionCandidates(AliVEvent::kEMC1|AliVEvent::kEMC7|AliVEvent::kEMC8|AliVEvent::kEMCEJE|AliVEvent::kEMCEGA);
  //taskmbemcal->SelectCollisionCandidates(AliVEvent::kAnyINT);
  //taskmbemcal->SetDebugLevel(10);

  // Add task(s)
  // taskmbemcal->PrintInfo();
  if(isPhysicsSelection) mgr->AddTask(physselTask);
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
