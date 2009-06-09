void runLRCAnalysis(const char* mode = "GRID", const char* inputName= "wn.xml") {
// This macro runs AliAnalysisTaskLRC in three modes : "Local" , "Interactive", "PROOF", "GRID" 
// ESD-only
// inputName refers to :
// "Local" - file with plane-text list of filenames 
// "Intaractive","GRID" - XML collection name
// "PROOF" - dataset name 
//  ---------------------------------------------------------
// This macro needs AliRoot-v4-17-03 or later 
// and macro AddTaskLRC.C in workdir

// For PROOF run PARS(ANALYSISalice.par  ANALYSIS.par  AOD.par  ESD.par  PWG2ebye.par  
// STEERBase.par ) from AliRoot-v4-17-03 or later are to be in workdir.

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
 
  //gROOT->ProcessLine(".include $ALICE_ROOT/include");
  //cout<<"  # Compiling analysis task\n";
  //gROOT->LoadMacro("AliAnalysisTaskLRC.cxx+g");
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  gROOT->LoadMacro("AddTaskLRC.C");
  TChain* chain = CreateESDChain(inputName);
  
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestLRCManagerLocal");
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);  
  
  AddLRCTaskSet();
  
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
		

gProof->UploadPackage("STEERBase");
gProof->UploadPackage("ESD");
gProof->UploadPackage("AOD");
gProof->UploadPackage("ANALYSIS");
gProof->UploadPackage("ANALYSISalice");
gProof->UploadPackage("PWG2ebye");


gProof->EnablePackage("STEERBase");
gProof->EnablePackage("ESD");
gProof->EnablePackage("AOD");
gProof->EnablePackage("ANALYSIS");
gProof->EnablePackage("ANALYSISalice");
gProof->EnablePackage("PWG2ebye");

 gROOT->LoadMacro("AddTaskLRC.C");  
  
  // Use AliRoot includes to compile our task
  //gROOT->ProcessLine(".include $ALICE_ROOT/include");

  //cout<<"  # Compiling analysis task\n";
  
  
  //gProof->Load("AliAnalysisTaskLRC.cxx++g");   
  // Create the analysis manager
  mgr = new AliAnalysisManager("TestLRCManagerProof");
   
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  AddLRCTaskSet();
   
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  mgr->StartAnalysis("proof", inputName);

  

};

void runLRCInteractive(const char* inputName= "wn.xml") {
  
  printf("  ------------------------------------------\n");
  printf(" # LRC local-interactive/GRID run manager \n");
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
  
  //gROOT->ProcessLine(".include $ALICE_ROOT/include");
  
  
 //___________Compile analysis task using AClic____________//
 cout<<"  # Compiling analysis task\n";
 // gROOT->LoadMacro("AliAnalysisTaskLRC.cxx+g");
  
  
   
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
 
  gROOT->LoadMacro("AddTaskLRC.C");
   
  AddLRCTaskSet();
  
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
