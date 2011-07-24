runLRCAnalysis(const char* mode = "GRID", const char* inputName= "wn.xml",Bool_t LoadTaskLocal=kTRUE) {
// This macro runs AliAnalysisTaskLRC in three modes : "Local" , "Interactive", "PROOF", "GRID" 
// ESD-only
// inputName refers to :
// "Local" - file with plane-text list of filenames 
// "Intaractive","GRID" - XML collection name
// "PROOF" - dataset name 
//  ---------------------------------------------------------
// This macro needs AliRoot-v4-17 or later 
// and macro AddTaskLRC.C in workdir

// For PROOF run PARS(ANALYSISalice.par  ANALYSIS.par  AOD.par  ESD.par  PWG2ebye.par  
// STEERBase.par ) are to be in workdir.


// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

// Version line : 3.6
// Version 3.6.6


if(mode!="Local" && mode!="Interactive" && mode!="PROOF" && mode!="GRID")
{
cout<<" ! Mode must be : Local , Interactive, PROOF, GRID \n";
cout<<" ! Unknown mode :"<<mode<< " \n";
return;
}

if(mode=="Local")runLRCLocal(inputName,LoadTaskLocal);
if(mode=="PROOF")runLRCProof(inputName,LoadTaskLocal);
if(mode=="Interactive")runLRCInteractive(inputName,LoadTaskLocal);
if(mode=="GRID")runLRCInteractive(inputName,LoadTaskLocal);


}

void LoadAnalysisLibs(Bool_t LoadTaskLocal=kFALSE)
{
  cout<<"  # Loadnig libs...\n";
 

 gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  if(!LoadTaskLocal)gSystem->Load("libPWG2ebye.so");

  
  //___________Compile analysis task using AClic____________//
 
  if(LoadTaskLocal){
  	gROOT->ProcessLine(".include $ALICE_ROOT/include");
	cout<<"  # Compiling AliLRCBase\n";
	gROOT->LoadMacro("AliLRCBase.cxx+g");	
	cout<<"  # Compiling AliLRCProcess\n";
	gROOT->LoadMacro("AliLRCProcess.cxx+g");
	gROOT->LoadMacro("AliRidgeAnalyser.cxx+g");
	cout<<"  # Compiling LRC analysis task\n";
  	gROOT->LoadMacro("AliAnalysisTaskLRC.cxx+g");
	}

}

Bool_t CreateLRCManager(char* name="LRCmanager",Bool_t runKine=kFALSE,Bool_t runAOD=kFALSE)
{
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager(name);
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);  
  if(runKine){
  	AliMCEventHandler* handler = new AliMCEventHandler;
  	mgr->SetMCtruthEventHandler(handler);
  }
  
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  
  mgr->SetCommonFileName("Test.LRC.root");
  
  

 return kTRUE;

}

void runLRCLocal(const char* inputName= "ESDs.lst",Bool_t LoadTaskLocal=kFALSE,Bool_t runKine=kFALSE,Bool_t runAOD=kFALSE) {
  printf("  ------------------------------------------\n");
  printf("  # LRC local run manager \n");
  cout<<"  # Task from :"<<gSystem->pwd()<<"\n";
  cout<<"  # Input list :"<<inputName<<"\n";
  
  
  TStopwatch timer;
  timer.Start();
  
  LoadAnalysisLibs(LoadTaskLocal);
 
  if (!CreateLRCManager("LocalLRCTest",runKine,runAOD)) return;
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();

  gROOT->LoadMacro("AddTaskLRC.C");
  AddTaskLRC(runKine);
  
 
  
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain(inputName);
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}

void runLRCProof(const char* inputName= "/COMMON/COMMON/tutorial_small",Bool_t LoadTaskLocal=kFALSE,const char* proofLink="anivanov@alice-caf.cern.ch",Bool_t runKine=kFALSE,Bool_t runAOD=kFALSE)
{
  printf("  ------------------------------------------\n");
  printf(" # LRC PROOF run manager \n");
  cout<<"  # Task from :"<<gSystem->pwd()<<"\n";
  cout<<"  # Dataset :"<<inputName<<"\n";


	TProof::Open(proofLink);
	//TProof::Open("anivanov@localhost");.

	
  cout<<"  # Loadnig libs...\n";



gProof->EnablePackage("VO_ALICE@AliRoot::v4-20-13-AN");

  // Use AliRoot includes to compile our task
  if(LoadTaskLocal){
 //	gROOT->ProcessLine(".include $ALICE_ROOT/include");
	cout<<"  # Compiling AliLRCBase\n";
	gProof->Load("AliLRCBase.cxx+g");	
	cout<<"  # Compiling AliLRCProcess\n";
	gProof->Load("AliLRCProcess.cxx+g");
	gProof->Load("AliLRCAnalyser.cxx+g");
  	cout<<"  # Compiling analysis task\n";
  	gProof->Load("AliAnalysisTaskLRC.cxx+g");   
  }

  if (!CreateLRCManager("ProofLRCTest",runKine,runAOD)) return;
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
   
  gROOT->LoadMacro("AddTaskLRC.C");
  AddTaskLRC(runKine);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("proof", inputName);
};

void runLRCInteractive(const char* inputName= "wn.xml",Bool_t LoadTaskLocal=kFALSE,Bool_t runKine=kFALSE,Bool_t runAOD=kFALSE) {
  
  printf("  ------------------------------------------\n");
  printf(" # LRC local-interactive/GRID run manager \n");
  cout<<"  # Task from :"<<gSystem->pwd()<<"\n";
  cout<<"  # Collection :"<<inputName<<"\n";

  TStopwatch timer;
  timer.Start();
  
  cout<<"*** Connect to AliEn ***\n";
  TGrid::Connect("alien://","anivanov");

  LoadAnalysisLibs(LoadTaskLocal);
  
  if (!CreateLRCManager("IntLRCTest",runKine,runAOD)) return;
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
 
  
  gROOT->LoadMacro("AddTaskLRC.C");
  AddTaskLRC(runKine);

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
    
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}

Int_t SetupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename) {
    char processline[1024];
    TString base = gSystem->BaseName(pararchivename);
    TString dir  = gSystem->DirName(pararchivename);
    TString ocwd = gSystem->WorkingDirectory();
    // Move to dir where the par files are and unpack
    gSystem->ChangeDirectory(dir.Data());
    sprintf(processline,".! tar xvzf %s.par",base.Data());
    gROOT->ProcessLine(processline);
    // Move to par folder
    gSystem->ChangeDirectory(base.Data());
                                                                                                                                               
    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");
                                                                                                                                               
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      // If dir not empty, set the full include path
      if (dir.Length()) {
         sprintf(processline, ".include %s", pararchivename);
         gROOT->ProcessLine(processline); 
      }   
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
                                                                                                                                               
    gSystem->ChangeDirectory(ocwd.Data());
  }                                                                                                                                               
  return 1;
}
