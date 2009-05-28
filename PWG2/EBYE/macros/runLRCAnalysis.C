void runLRCAnalysis(const char* mode = "Local", 
		    const char* inputName= "ESDs.lst") {
// This macro runs AliAnalysisTaskLRC in three modes : "Local" , "Interactive", "PROOF", "GRID" 
// ESD-only
// inputName refers to :
// "Local" - file with plane-text list of filenames 
// "Intaractive","GRID" - XML collection name
// "PROOF" - dataset name 
//  
//

if(mode!="Local" && mode!="Interactive" && mode!="PROOF" && mode!="GRID")
{
cout<<" ! Mode must be : Local , Interactive, PROOF, GRID \n";
cout<<" ! Unknown mode :"<<mode<< " \n";
return;
}

if(mode=="Local")runLRCLocal(inputName);
if(mode=="PROOF")runLRCProof(inputName);
if(mode=="Interactive")runLRCInteractive(inputName);
if(mode=="GRID")runLRCInteractive(inputName);


}

void runLRCLocal(const char* inputName= "ESDs.lst") {
  printf("  ------------------------------------------\n");
  printf("  # LRC local run manager \n");
  cout<<"  # Task from :"<<gSystem->pwd()<<"\n";
  cout<<"  # Input list :"<<inputName<<"\n";
  cout<<"  # Loadnig libs...\n";
  
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG2ebye.so");
  
  
  //___________Compile analysis task using AClic____________//
 
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  cout<<"  # Compiling analysis task\n";
  //gROOT->LoadMacro("AliAnalysisTaskLRC.cxx+g");
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");

  TChain* chain = CreateESDChain(inputName);
  
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);  
  
  
  AliAnalysisTaskLRC *task1 = new AliAnalysisTaskLRC("TaskLRC");
  mgr->AddTask(task1);
  
  task1->SetETAWindows(-0.7,0.0,0.3,1.0);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TList::Class(),AliAnalysisManager::kOutputContainer,"LRC.ESD.testLocal.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}

void runLRCProof(const char* inputName= "/COMMON/COMMON/tutorial_small")
{
  printf("  ------------------------------------------\n");
  printf(" # LRC PROOF run manager \n");
  cout<<"  # Task from :"<<gSystem->pwd()<<"\n";
  cout<<"  # Dataset :"<<inputName<<"\n";


	TProof::Open("alicecaf");
	//TProof::Open("anivanov@localhost");

	
  cout<<"  # Loadnig libs...\n";
		
	gProof->UploadPackage("AF-v4-16");
	gProof->EnablePackage("AF-v4-16");

	  
  
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");

  cout<<"  # Compiling analysis task\n";
  
  
  //gProof->Load("AliAnalysisTaskLRC.cxx++g");   
  task = new AliAnalysisTaskLRC("TestLRC");
  
  task->SetETAWindows(-0.7,0.0,0.3,1.0);

    // Create the analysis manager
  mgr = new AliAnalysisManager("testAnalysis");
  
  
  mgr->AddTask(task);

  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);


  // Create containers for input/output
  cinput = mgr->CreateContainer("cchain", TChain::Class(), AliAnalysisManager::kInputContainer);
  coutput = mgr->CreateContainer("chist", TList::Class(),    AliAnalysisManager::kOutputContainer, "LRC.ESD.PROOFtest.root");
 
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 0, coutput);
 
 
  mgr->InitAnalysis();
  mgr->PrintStatus();

 mgr->StartAnalysis("proof", inputName);

  

};

void runLRCInteractive(const char* inputName= "wn.xml") {
  
  printf("  ------------------------------------------\n");
  printf(" # LRC local-interactive run manager \n");
  cout<<"  # Task from :"<<gSystem->pwd()<<"\n";
  cout<<"  # Collection :"<<inputName<<"\n";

  cout<<"*** Connect to AliEn ***\n";
  TGrid::Connect("alien://","anivanov");

  cout<<"  # Loadnig libs...\n";
  
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  
  gSystem->Load("libPWG2ebye.so");
  
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  
  
 //___________Compile analysis task using AClic____________//
 cout<<"  # Compiling analysis task\n";
 //gROOT->LoadMacro("AliAnalysisTaskLRC.cxx+g");
  
  
   
    TAlienCollection * myCollection = 
	new TAlienCollection("wn.xml",100000) ; 
    if (!myCollection) { 
	cout << "XML collection file: " << xmlFileName << " not found" << endl;
	return;
    }

    TChain* chain = new TChain("esdTree");
    
    cout << "Preparing the file list" << endl; 
    myCollection->Reset() ; 
    while ( myCollection->Next() ) {
	cout << "Adding ESD file: " << myCollection->GetTURL("") << endl; 
	chain->Add(myCollection->GetTURL("")) ; 
    }
    

  
  AliAnalysisManager *mgr = new AliAnalysisManager("Local-Interactive_LRC_Manager");
  AliESDInputHandler* esdH = new AliESDInputHandler();
  
  mgr->SetInputEventHandler(esdH);  
  //____________________________________________//
 
  AliAnalysisTaskLRC *task1 = new AliAnalysisTaskLRC("TaskLRC");
  mgr->AddTask(task1);
  
  task1->SetETAWindows(-0.7,0.0,0.3,1.0);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TList::Class(),AliAnalysisManager::kOutputContainer,"LRC.ESD.Local-Interactive.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}
