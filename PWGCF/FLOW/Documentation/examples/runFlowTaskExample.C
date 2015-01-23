// Analysis type can be ESD, AOD, MC, ESDMCkineESD, ESDMCkineMC
// in this simple example only ESD is supported
const TString type = "ESD";
// tracks used to determine the Q vector (Global, Tracklet or FMD at the moment)
const TString rptype = "Global";
// Boolean to fill/not fill the QA histograms
Bool_t QA = kTRUE;   

// Boolean to use/not use weights for the Q vector
Bool_t WEIGHTS[] = {kFALSE,kFALSE,kFALSE}; //Phi, v'(pt), v'(eta)


void runFlowTaskExample(Int_t nRuns = 2, 
              Bool_t DATA = kFALSE, const Char_t* dataDir="/Users/snelling/alice_data/Therminator_midcentral", Int_t offset = 0)
{
  TStopwatch timer;
  timer.Start();
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");


  if (type!="AOD") { TChain* chain = CreateESDChain(dataDir, nRuns, offset);}
  else { TChain* chain = CreateAODChain(dataDir, nRuns, offset);}


  // CORRFW correction framework cuts
  //----------Event cuts----------

  TObjArray* mcEventList = new TObjArray(0); 
  AliCFEventGenCuts* mcEventCuts = new AliCFEventGenCuts("mcEventCuts","MC-level event cuts");
  mcEventList->AddLast(mcEventCuts);

  TObjArray* recEventList = new TObjArray(0);
  AliCFEventRecCuts* recEventCuts = new AliCFEventRecCuts("recEventCuts","rec-level event cuts");
  recEventList->AddLast(recEventCuts);

  //----------Cuts for RP----------
  TObjArray* mcListRP = new TObjArray(0);
  AliCFTrackKineCuts *mcKineCutsRP = new AliCFTrackKineCuts("mcKineCutsRP","MC-level kinematic cuts");
  mcListRP->AddLast(mcKineCutsRP);
  AliCFParticleGenCuts *mcGenCutsRP = new AliCFParticleGenCuts("mcGenCutsRP","MC particle generation cuts for RP");
  mcGenCutsRP->SetRequireIsPrimary();
  mcListRP->AddLast(mcGenCutsRP);

  TObjArray* fPIDCutListRP = new TObjArray(0) ;
  AliCFTrackCutPid* cutPidRP = NULL;
  fPIDCutListRP->AddLast(cutPidRP);

  TObjArray* recListRP = new TObjArray(0) ;
  AliCFTrackKineCuts *recKineCutsRP = new AliCFTrackKineCuts("recKineCutsRP","rec-level kine cuts");
  recListRP->AddLast(recKineCutsRP); 
  AliCFTrackQualityCuts *recQualityCutsRP = new AliCFTrackQualityCuts("recQualityCutsRP","rec-level quality cuts");
  recQualityCutsRP->SetMinNClusterTPC(70);
  recListRP->AddLast(recQualityCutsRP);
  AliCFTrackIsPrimaryCuts *recIsPrimaryCutsRP = new AliCFTrackIsPrimaryCuts("recIsPrimaryCutsRP","rec-level isPrimary cuts");
  recIsPrimaryCutsRP->UseSPDvertex(kTRUE);
  recIsPrimaryCutsRP->UseTPCvertex(kTRUE);
  recListRP->AddLast(recIsPrimaryCutsRP);

  TObjArray* accListRP = new TObjArray(0) ;
  AliCFAcceptanceCuts *mcAccCutsRP = new AliCFAcceptanceCuts("mcAccCutsRP","MC acceptance cuts");
  mcAccCutsRP->SetMinNHitITS(2);
  accListRP->AddLast(mcAccCutsRP);

  //----------Cuts for POI----------
  TObjArray* mcListPOI = new TObjArray(0);
  AliCFTrackKineCuts *mcKineCutsPOI = new AliCFTrackKineCuts("mcKineCutsPOI","MC-level kinematic cuts");
  mcListPOI->AddLast(mcKineCutsPOI);
  AliCFParticleGenCuts *mcGenCutsPOI = new AliCFParticleGenCuts("mcGenCutsPOI","MC particle generation cuts for POI");
  mcGenCutsPOI->SetRequireIsPrimary();
  mcListPOI->AddLast(mcGenCutsPOI);

  TObjArray* fPIDCutListPOI = new TObjArray(0) ;
  AliCFTrackCutPid* cutPidPOI = NULL;
  fPIDCutListPOI->AddLast(cutPidPOI);

  TObjArray* recListPOI = new TObjArray(0) ;
  AliCFTrackKineCuts *recKineCutsPOI = new AliCFTrackKineCuts("recKineCutsPOI","rec-level kine cuts");
  recListPOI->AddLast(recKineCutsPOI); 
  AliCFTrackQualityCuts *recQualityCutsPOI = new AliCFTrackQualityCuts("recQualityCutsPOI","rec-level quality cuts");
  recQualityCutsPOI->SetMinNClusterTPC(70);
  recListPOI->AddLast(recQualityCutsPOI);
  AliCFTrackIsPrimaryCuts *recIsPrimaryCutsPOI = new AliCFTrackIsPrimaryCuts("recIsPrimaryCutsPOI","rec-level isPrimary cuts");
  recIsPrimaryCutsPOI->UseSPDvertex(kTRUE);
  recIsPrimaryCutsPOI->UseTPCvertex(kTRUE);
  recListPOI->AddLast(recIsPrimaryCutsPOI);

  TObjArray* accListPOI = new TObjArray(0) ;
  AliCFAcceptanceCuts *mcAccCutsPOI = new AliCFAcceptanceCuts("mcAccCutsPOI","MC acceptance cuts");
  mcAccCutsPOI->SetMinNHitITS(2);
  accListPOI->AddLast(mcAccCutsPOI);

  //----------Add Cut Lists to the CF Manager----------
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* cfmgrRP = new AliCFManager();
  cfmgrRP->SetNStepEvent(3);
  cfmgrRP->SetEventCutsList(AliCFManager::kEvtGenCuts,mcEventList); 
  cfmgrRP->SetEventCutsList(AliCFManager::kEvtRecCuts,recEventList); 
  cfmgrRP->SetNStepParticle(4); 
  cfmgrRP->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListRP);
  cfmgrRP->SetParticleCutsList(AliCFManager::kPartAccCuts,accListRP);
  cfmgrRP->SetParticleCutsList(AliCFManager::kPartRecCuts,recListRP);
  cfmgrRP->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutListRP);

  AliCFManager* cfmgrPOI = new AliCFManager();
  cfmgrPOI->SetNStepEvent(3);
  cfmgrPOI->SetEventCutsList(AliCFManager::kEvtGenCuts,mcEventList); 
  cfmgrPOI->SetEventCutsList(AliCFManager::kEvtRecCuts,recEventList); 
  cfmgrPOI->SetNStepParticle(4); 
  cfmgrPOI->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListPOI);
  cfmgrPOI->SetParticleCutsList(AliCFManager::kPartAccCuts,accListPOI);
  cfmgrPOI->SetParticleCutsList(AliCFManager::kPartRecCuts,recListPOI);
  cfmgrPOI->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutListPOI);


  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager");
  
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc); 


  //  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"); 
  //  AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection();
  //  if (!DATA) {physicsSelTask->GetPhysicsSelection()->SetAnalyzeMC();}
  
  // Add the task which makes the flowevent to the analysis manager 
  AliAnalysisTaskFlowEvent* taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",rptype,kFALSE,1);
  taskFE->SetAnalysisType(type);  
  //  taskFE->SelectCollisionCandidates();  
  mgr->AddTask(taskFE);

  // Set the cuts for the flowevent
  taskFE->SetCFManager1(cfmgrRP);
  taskFE->SetCFManager2(cfmgrPOI);

  // Instantiate the flow methods and add to the analysis manager
  AliAnalysisTaskMCEventPlane *taskMCEP = new AliAnalysisTaskMCEventPlane("TaskMCEventPlane");
  mgr->AddTask(taskMCEP);
  AliAnalysisTaskQCumulants *taskQC = new AliAnalysisTaskQCumulants("TaskQCumulants",kFALSE);
  mgr->AddTask(taskQC);

  // Create the output container for the data produced by the task
  // Connect to the input and output containers
  //===========================================================================
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *coutputFE = 
    mgr->CreateContainer("cobjFlowEventSimple",  AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput(taskFE,0,cinput1); 
  mgr->ConnectOutput(taskFE,1,coutputFE);


  TString outputMCEP = AliAnalysisManager::GetCommonFileName();
  outputMCEP += ":outputMCEPanalysis";
  outputMCEP+= type;
    
  AliAnalysisDataContainer *coutputMCEP = mgr->CreateContainer("cobjMCEP", TList::Class(),AliAnalysisManager::kOutputContainer,outputMCEP); 
  mgr->ConnectInput(taskMCEP,0,coutputFE); 
  mgr->ConnectOutput(taskMCEP,1,coutputMCEP); 


  TString outputQC = AliAnalysisManager::GetCommonFileName();
  outputQC += ":outputQCanalysis";
  outputQC+= type;

  AliAnalysisDataContainer *coutputQC = mgr->CreateContainer("cobjQC", TList::Class(),AliAnalysisManager::kOutputContainer,outputQC); 
  mgr->ConnectInput(taskQC,0,coutputFE); 
  mgr->ConnectOutput(taskQC,1,coutputQC);

  // Run the analysis
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
  
  timer.Stop();
  timer.Print();
  
}

// Helper macros for creating chains
// from: CreateESDChain.C,v 1.10 jgrosseo Exp

TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...
  
  if (!aDataDir)
    return 0;
  
  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
    {
      printf("%s not found.\n", aDataDir);
      return 0;
    }
  
  TChain* chain = new TChain("esdTree");
  TChain* chaingAlice = 0;
  
  if (flags & 2)
    {
      TString execDir(gSystem->pwd());
      TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
      TList* dirList            = baseDir->GetListOfFiles();
      Int_t nDirs               = dirList->GetEntries();
      gSystem->cd(execDir);
      
      Int_t count = 0;
      
      for (Int_t iDir=0; iDir<nDirs; ++iDir)
	{
	  TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
	  if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
	    continue;
	  
	  if (offset > 0)
	    {
	      --offset;
	      continue;
	    }
	  
	  if (count++ == aRuns)
	    break;
	  
	  TString presentDirName(aDataDir);
	  presentDirName += "/";
	  presentDirName += presentDir->GetName();	  
	  chain->Add(presentDirName + "/AliESDs.root/esdTree");
	  //  cerr<<presentDirName<<endl;
	}
      
    }
  else
    {
      // Open the input stream
      ifstream in;
      in.open(aDataDir);
      
      Int_t count = 0;
      
      // Read the input list of files and add them to the chain
      TString esdfile;
      while(in.good()) {
	in >> esdfile;
	if (!esdfile.Contains("root")) continue; // protection
	
	if (offset > 0)
	  {
	    --offset;
	    continue;
	  }
	
	if (count++ == aRuns)
	  break;
	
	// add esd file
	chain->Add(esdfile);
      }
      
      in.close();
    }
  
  return chain;
}


// Helper macros for creating chains
// from: CreateESDChain.C,v 1.10 jgrosseo Exp

TChain* CreateAODChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliAOD.root
  // <aDataDir>/<dir1>/AliAOD.root
  // ...
  
  if (!aDataDir)
    return 0;
  
  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
    {
      printf("%s not found.\n", aDataDir);
      return 0;
    }
  
  TChain* chain = new TChain("aodTree");
  TChain* chaingAlice = 0;
  
  if (flags & 2)
    {
      TString execDir(gSystem->pwd());
      TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
      TList* dirList            = baseDir->GetListOfFiles();
      Int_t nDirs               = dirList->GetEntries();
      gSystem->cd(execDir);
      
      Int_t count = 0;
      
      for (Int_t iDir=0; iDir<nDirs; ++iDir)
	{
	  TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
	  if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
	    continue;
	  
	  if (offset > 0)
	    {
	      --offset;
	      continue;
	    }
	  
	  if (count++ == aRuns)
	    break;
	  
	  TString presentDirName(aDataDir);
	  presentDirName += "/";
	  presentDirName += presentDir->GetName();	  
	  chain->Add(presentDirName + "/AliAOD.root/aodTree");
	  // cerr<<presentDirName<<endl;
	}
      
    }
  else
    {
      // Open the input stream
      ifstream in;
      in.open(aDataDir);
      
      Int_t count = 0;
      
      // Read the input list of files and add them to the chain
      TString aodfile;
      while(in.good()) {
	in >> aodfile;
	if (!aodfile.Contains("root")) continue; // protection
	
	if (offset > 0)
	  {
	    --offset;
	    continue;
	  }
	
	if (count++ == aRuns)
	  break;
	
	// add aod file
	chain->Add(aodfile);
      }
      
      in.close();
    }
  
  return chain;
}

