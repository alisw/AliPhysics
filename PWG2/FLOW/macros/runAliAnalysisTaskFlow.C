// from CreateESDChain.C - instead of  #include "CreateESDChain.C"
TChain* CreateESDChain(const char* aDataDir = "ESDfiles.txt", Int_t aRuns = 20, Int_t offset = 0) ;
void LookupWrite(TChain* chain, const char* target) ;

//SETTING THE ANALYSIS

//Flow analysis method can be:
// SP    = Scalar Product
// LYZ1  = Lee Yang Zeroes first run
// LYZ2  = Lee Yang Zeroes second run
// LYZEP = Lee Yang Zeroes Event Plane
// CUM   = Cumulants
// MCEP  = Flow calculated from the real MC event plane (only for simulated data)
const TString method = "MCEP";

//Type of analysis can be:
// ESD, AOD, MC, ESDMC0, ESDMC1
const TString type = "ESD";


//SETTING THE CUTS

//for integrated flow
const Double_t ptmin1 = 0.0;
const Double_t ptmax1 = 1000.0;
const Double_t ymin1  = -2.;
const Double_t ymax1  = 2.;
const Int_t mintrackrefsTPC1 = 2;
const Int_t mintrackrefsITS1 = 3;
const Int_t charge1 = 1;
const Int_t PDG1 = 211;
const Int_t minclustersTPC1 = 50;
const Int_t maxnsigmatovertex1 = 3;

//for differential flow
const Double_t ptmin2 = 0.0;
const Double_t ptmax2 = 1000.0;
const Double_t ymin2  = -2.;
const Double_t ymax2  = 2.;
const Int_t mintrackrefsTPC2 = 2;
const Int_t mintrackrefsITS2 = 3;
const Int_t charge2 = 1;
const Int_t PDG2 = 321;
const Int_t minclustersTPC2 = 50;
const Int_t maxnsigmatovertex2 = 3;



void runAliAnalysisTaskFlow(Int_t nRuns = 10, const Char_t* dataDir="/data/alice1/kolk/TherminatorFIX", Int_t offset = 0) 

{
  TStopwatch timer;
  timer.Start();

  // include path (to find the .h files when compiling)
  gSystem->AddIncludePath("-I$ALICE_ROOT/include") ;
  gSystem->AddIncludePath("-I$ROOTSYS/include") ;

  // load needed libraries
  gSystem->Load("libTree.so");
  gSystem->Load("libESD.so");
  cerr<<"libESD loaded..."<<endl;
  gSystem->Load("libANALYSIS.so");
  cerr<<"libANALYSIS.so loaded..."<<endl;
  gSystem->Load("libANALYSISRL.so");
  cerr<<"libANALYSISRL.so loaded..."<<endl;
  gSystem->Load("libCORRFW.so");
  cerr<<"libCORRFW.so loaded..."<<endl;
  gSystem->Load("libPWG2flow.so");
  cerr<<"libPWG2flow.so loaded..."<<endl;

  // create the TChain. CreateESDChain() is defined in CreateESDChain.C
  TChain* chain = CreateESDChain(dataDir, nRuns, offset);
  cout<<"chain ("<<chain<<")"<<endl;
 
  //____________________________________________//
  //Create cuts using correction framework

  //############# cuts on MC
  AliCFTrackKineCuts* mcKineCuts1 = new AliCFTrackKineCuts("mcKineCuts1","MC-level kinematic cuts");
  mcKineCuts1->SetPtRange(ptmin1,ptmax1);
  mcKineCuts1->SetRapidityRange(ymin1,ymax1);
  mcKineCuts1->SetChargeMC(charge1);
  
  AliCFTrackKineCuts* mcKineCuts2 = new AliCFTrackKineCuts("mcKineCuts2","MC-level kinematic cuts");
  mcKineCuts2->SetPtRange(ptmin2,ptmax2);
  mcKineCuts2->SetRapidityRange(ymin2,ymax2);
  mcKineCuts2->SetChargeMC(charge2);

  AliCFParticleGenCuts* mcGenCuts1 = new AliCFParticleGenCuts("mcGenCuts1","MC particle generation cuts");
  mcGenCuts1->SetRequireIsPrimary();
  mcGenCuts1->SetRequirePdgCode(PDG1);

  AliCFParticleGenCuts* mcGenCuts2 = new AliCFParticleGenCuts("mcGenCuts2","MC particle generation cuts");
  mcGenCuts2->SetRequireIsPrimary();
  mcGenCuts2->SetRequirePdgCode(PDG2);

  //############# Acceptance Cuts  ????????
  AliCFAcceptanceCuts *mcAccCuts = new AliCFAcceptanceCuts("mcAccCuts","MC acceptance cuts");
  mcAccCuts->SetMinNHitITS(mintrackrefsITS1);
  mcAccCuts->SetMinNHitTPC(mintrackrefsTPC1);
  
  //############# Rec-Level kinematic cuts
  AliCFTrackKineCuts *recKineCuts1 = new AliCFTrackKineCuts("recKineCuts1","rec-level kine cuts");
  recKineCuts1->SetPtRange(ptmin1,ptmax1);
  recKineCuts1->SetRapidityRange(ymin1,ymax1);
  recKineCuts1->SetChargeRec(charge1);

  AliCFTrackKineCuts *recKineCuts2 = new AliCFTrackKineCuts("recKineCuts2","rec-level kine cuts");
  recKineCuts2->SetPtRange(ptmin2,ptmax2);
  recKineCuts2->SetRapidityRange(ymin2,ymax2);
  recKineCuts2->SetChargeRec(charge2);

  AliCFTrackQualityCuts *recQualityCuts = new AliCFTrackQualityCuts("recQualityCuts","rec-level quality cuts");
  recQualityCuts->SetMinNClusterTPC(minclustersTPC1);
  recQualityCuts->SetRequireITSRefit(kTRUE);

  AliCFTrackIsPrimaryCuts *recIsPrimaryCuts = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts","rec-level isPrimary cuts");
  recIsPrimaryCuts->SetMaxNSigmaToVertex(maxnsigmatovertex1);
  
  AliCFTrackCutPid* cutPID1 = new AliCFTrackCutPid("cutPID1","ESD_PID") ;
  AliCFTrackCutPid* cutPID2 = new AliCFTrackCutPid("cutPID2","ESD_PID") ;
  int n_species = AliPID::kSPECIES ;
  Double_t* prior = new Double_t[n_species];

  prior[0] = 0.0244519 ;
  prior[1] = 0.0143988 ;
  prior[2] = 0.805747  ;
  prior[3] = 0.0928785 ;
  prior[4] = 0.0625243 ;
  
  cutPID1->SetPriors(prior);
  cutPID1->SetProbabilityCut(0.0);
  cutPID1->SetDetectors("TPC TOF");
  switch(TMath::Abs(PDG1)) {
  case 11   : cutPID1->SetParticleType(AliPID::kElectron, kTRUE); break;
  case 13   : cutPID1->SetParticleType(AliPID::kMuon    , kTRUE); break;
  case 211  : cutPID1->SetParticleType(AliPID::kPion    , kTRUE); break;
  case 321  : cutPID1->SetParticleType(AliPID::kKaon    , kTRUE); break;
  case 2212 : cutPID1->SetParticleType(AliPID::kProton  , kTRUE); break;
  default   : printf("UNDEFINED PID\n"); break;
  }

  cutPID2->SetPriors(prior);
  cutPID2->SetProbabilityCut(0.0);
  cutPID2->SetDetectors("TPC TOF");
  switch(TMath::Abs(PDG2)) {
  case 11   : cutPID2->SetParticleType(AliPID::kElectron, kTRUE); break;
  case 13   : cutPID2->SetParticleType(AliPID::kMuon    , kTRUE); break;
  case 211  : cutPID2->SetParticleType(AliPID::kPion    , kTRUE); break;
  case 321  : cutPID2->SetParticleType(AliPID::kKaon    , kTRUE); break;
  case 2212 : cutPID2->SetParticleType(AliPID::kProton  , kTRUE); break;
  default   : printf("UNDEFINED PID\n"); break;
  }


  printf("CREATE MC KINE CUTS\n");
  TObjArray* mcList1 = new TObjArray(0);
  mcList1->AddLast(mcKineCuts1);
  mcList1->AddLast(mcGenCuts1);
  
  TObjArray* mcList2 = new TObjArray(0);
  mcList2->AddLast(mcKineCuts2);
  mcList2->AddLast(mcGenCuts2);

  printf("CREATE ACCEPTANCE CUTS\n");
  TObjArray* accList = new TObjArray(0) ;
  accList->AddLast(mcAccCuts);

  printf("CREATE RECONSTRUCTION CUTS\n");
  TObjArray* recList1 = new TObjArray(0) ;
  recList1->AddLast(recKineCuts1);
  recList1->AddLast(recQualityCuts);
  recList1->AddLast(recIsPrimaryCuts);

  TObjArray* recList2 = new TObjArray(0) ;
  recList2->AddLast(recKineCuts2);
  recList2->AddLast(recQualityCuts);
  recList2->AddLast(recIsPrimaryCuts);

  printf("CREATE PID CUTS\n");
  TObjArray* fPIDCutList1 = new TObjArray(0) ;
  fPIDCutList1->AddLast(cutPID1);

  TObjArray* fPIDCutList2 = new TObjArray(0) ;
  fPIDCutList2->AddLast(cutPID2);

  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* cfmgr1 = new AliCFManager();
  cfmgr1->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList1);
  //cfmgr1->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  cfmgr1->SetParticleCutsList(AliCFManager::kPartRecCuts,recList1);
  cfmgr1->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList1);

  AliCFManager* cfmgr2 = new AliCFManager();
  cfmgr2->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList2);
  //cfmgr2->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  cfmgr2->SetParticleCutsList(AliCFManager::kPartRecCuts,recList2);
  cfmgr2->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList2);

 
  if (method == "LYZ2"){  
    // read the input file from the first run 
    TString inputFileName = "outputLYZ1analysis" ;
    inputFileName += type;
    inputFileName += "_firstrun.root";
    cout<<"The input file is "<<inputFileName.Data()<<endl;
    TFile* fInputFile = new TFile(inputFileName.Data(),"READ");
    if(!fInputFile || fInputFile->IsZombie()) { cerr << " ERROR: NO First Run file... " << endl ; }
    cout<<"input file read..."<<endl;
  }

  if (method == "LYZEP") {
    // read the input file from the second LYZ run
    TString inputFileName = "outputLYZ2analysis" ;
    inputFileName += type;
    inputFileName += "_secondrun.root";
    cout<<"The input file is "<<inputFileName.Data()<<endl;
    TFile* fInputFile = new TFile(inputFileName.Data(),"READ");
    if(!fInputFile || fInputFile->IsZombie()) { cerr << " ERROR: NO First Run file... " << endl ; }
    cout<<"input file read..."<<endl;
  }
    

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager");
  
  if (type == "ESD"){
    AliVEventHandler* esdH = new AliESDInputHandler;
    mgr->SetInputEventHandler(esdH); 
    if (method == "MCEP") {
      AliMCEventHandler *mc = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mc);}  }
   
  if (type == "AOD"){
    AliVEventHandler* aodH = new AliAODInputHandler;
    mgr->SetInputEventHandler(aodH);
    if (method == "MCEP") {
      AliMCEventHandler *mc = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mc);}  }
  
  if (type == "MC" || type == "ESDMC0" || type == "ESDMC1"){
    AliVEventHandler* esdH = new AliESDInputHandler;
    mgr->SetInputEventHandler(esdH);

    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc); }

  //____________________________________________//
  // 1st MC EP task
  if (method == "SP"){
    AliAnalysisTaskScalarProduct *task1 = new AliAnalysisTaskScalarProduct("TaskScalarProduct");
    task1->SetAnalysisType(type);
    task1->SetCFManager1(cfmgr1);
    task1->SetCFManager2(cfmgr2);
    mgr->AddTask(task1);
  }
  else if (method == "LYZ1"){
    AliAnalysisTaskLeeYangZeros *task1 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kTRUE);
    task1->SetAnalysisType(type);
    task1->SetFirstRunLYZ(kTRUE);
    task1->SetUseSumLYZ(kTRUE);
    task1->SetCFManager1(cfmgr1);
    task1->SetCFManager2(cfmgr2);
    mgr->AddTask(task1);
  }
  else if (method == "LYZ2"){
    AliAnalysisTaskLeeYangZeros *task1 = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZeros",kFALSE);
    task1->SetAnalysisType(type);
    task1->SetFirstRunLYZ(kFALSE);
    task1->SetUseSumLYZ(kTRUE);
    task1->SetCFManager1(cfmgr1);
    task1->SetCFManager2(cfmgr2);
    mgr->AddTask(task1);
  }
  else if (method == "LYZEP"){
    AliAnalysisTaskLYZEventPlane *task1 = new AliAnalysisTaskLYZEventPlane("TaskLYZEventPlane");
    task1->SetAnalysisType(type);
    task1->SetCFManager1(cfmgr1);
    task1->SetCFManager2(cfmgr2);
    mgr->AddTask(task1);
  }
  else if (method == "CUM"){
    AliAnalysisTaskCumulants *task1 = new AliAnalysisTaskCumulants("TaskCumulants");
    task1->SetAnalysisType(type);
    //task1->SetCFManager1(cfmgr1);
    //task1->SetCFManager2(cfmgr2);
    mgr->AddTask(task1);
  }
  else if (method == "MCEP"){
    AliAnalysisTaskMCEventPlane *task1 = new AliAnalysisTaskMCEventPlane("TaskMCEventPlane");
    task1->SetAnalysisType(type);
    task1->SetCFManager1(cfmgr1);
    task1->SetCFManager2(cfmgr2);
    mgr->AddTask(task1);
  }


  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = 
    mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);

  if (method == "LYZ2" || method == "LYZEP"){ 
    AliAnalysisDataContainer *cinput2 = 
		    mgr->CreateContainer("cobj2",TList::Class(),AliAnalysisManager::kInputContainer); } 


  TString outputName = "output";
  outputName+= method;
  outputName+= "analysis";
  outputName+= type;
  if (method == "LYZ1") outputName+= "_firstrun";
  if (method == "LYZ2") outputName+= "_secondrun";
  outputName+= ".root";
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("cobj1", TList::Class(),AliAnalysisManager::kOutputContainer,outputName);


  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  if (method == "LYZ2" || method == "LYZEP") { 
    mgr->ConnectInput(task1,1,cinput2); }
  mgr->ConnectOutput(task1,0,coutput1);

  if (method == "LYZ2" || method == "LYZEP"){ 
    cinput2->SetData(fInputFile);}

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
	  cerr<<presentDirName<<endl;
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

void LookupWrite(TChain* chain, const char* target)
{
  // looks up the chain and writes the remaining files to the text file target
  
  chain->Lookup();
  
  TObjArray* list = chain->GetListOfFiles();
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  
  ofstream outfile;
  outfile.open(target);
  
  while ((obj = iter->Next()))
    outfile << obj->GetTitle() << "#AliESDs.root" << endl;
  
  outfile.close();
  
  delete iter;
}
