#include "TStopwatch.h"
#include "TObjArray"
#include "Riostream.h"
#include "TFile.h"

//--------------------------------------------------------------------------------------
// RUN SETTINGS
//flow analysis method can be: (set to kTRUE or kFALSE)
Bool_t SP    = kTRUE;
Bool_t LYZ1  = kTRUE;
Bool_t LYZ2  = kFALSE;  
Bool_t LYZEP = kFALSE; 
Bool_t GFC   = kTRUE;
Bool_t QC    = kTRUE;
Bool_t FQD   = kTRUE;
Bool_t MCEP  = kFALSE; //does not work yet 24/12/08
//--------------------------------------------------------------------------------------

// Weights 
// Use weights for Q vector
Bool_t usePhiWeights = kFALSE; //Phi (correction for non-uniform azimuthal acceptance)
Bool_t usePtWeights  = kFALSE; //v'(pt) (differential flow in pt)
Bool_t useEtaWeights = kFALSE; //v'(eta) (differential flow in eta)

//--------------------------------------------------------------------------------------
// CUT SETTINGS
//integrated selection
Double_t ptMaxInt  = 10.;
Double_t ptMinInt  = 0.;
Double_t etaMaxInt = 1.;
Double_t etaMinInt = -1.;
Double_t phiMaxInt = 7.5;
Double_t phiMinInt = 0.;
Int_t PIDInt       = 211;

//differential selection
Double_t ptMaxDiff  = 10.;
Double_t ptMinDiff  = 0.;
Double_t etaMaxDiff = 1.;
Double_t etaMinDiff = -1.;
Double_t phiMaxDiff = 7.5;
Double_t phiMinDiff = 0.;
Int_t PIDDiff       = 211;
//--------------------------------------------------------------------------------------


Int_t offset = 0;
                                          
int runFlowAnalysis(Int_t aRuns = 100, const char* 
		    // 			  dir="/data/alice1/kolk/KineOnly3/")
		    dir="/Users/snelling/alice_data/KineOnly3/")
{
  TStopwatch timer;
  timer.Start();
  
  if (LYZ1 && LYZ2) {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; exit(); }
  if (LYZ2 && LYZEP) {cout<<"WARNING: you cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }
  if (LYZ1 && LYZEP) {cout<<"WARNING: you cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }


  cout<<endl;
  cout<<" ---- BEGIN ANALYSIS ---- "<<endl;
  cout<<endl;
  

  gSystem->AddIncludePath("-I$ROOTSYS/include");

  // load needed libraries
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");

  
  // in root using pars

  SetupPar("PWG2flowCommon");
  cerr<<"PWG2flowCommon.par loaded..."<<endl;

  //aliroot specific stuff
  SetupPar("STEERBase");
  SetupPar("ESD");
  SetupPar("AOD");

  SetupPar("ANALYSIS");
  SetupPar("ANALYSISalice");
  SetupPar("PWG2AOD");

  SetupPar("CORRFW");

  SetupPar("PWG2flowTasks");
  cerr<<"PWG2flowTasks.par loaded..."<<endl;

  // AliFlow event in root with pars
  AliFlowEventSimpleMaker* fEventMaker = new AliFlowEventSimpleMaker();

  /*
  // In AliRoot
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libPWG2flowCommon.so");
  cerr<<"libPWG2flowCommon.so loaded ..."<<endl;
  gSystem->Load("libPWG2flowTasks.so");
  cerr<<"libPWG2flowTasks.so loaded ..."<<endl;
  cout<<endl; 
  
  // Flow event in AliRoot
  AliFlowEventSimpleMaker* fEventMaker = new AliFlowEventSimpleMaker();
   
  */

  // In root inline compile

  /*    
  // Constants  
  gROOT->LoadMacro("AliFlowCommon/AliFlowCommonConstants.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFlowLYZConstants.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFlowCumuConstants.cxx+");

  // Flow event
  gROOT->LoadMacro("AliFlowCommon/AliFlowVector.cxx+"); 
  gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimple.cxx+");    
  gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimple.cxx+");

  // Cuts
  gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimpleCuts.cxx+");    

  // Output histosgrams
  gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHist.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHistResults.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist1.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist2.cxx+");

  // Functions needed for various methods
  gROOT->LoadMacro("AliFlowCommon/AliCumulantsFunctions.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliQCumulantsFunctions.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFittingFunctionsForQDistribution.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFlowLYZEventPlane.cxx+");

  // Flow Analysis code for various methods
  gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithMCEventPlane.cxx+"); 
  gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithScalarProduct.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithLYZEventPlane.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithLeeYangZeros.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithCumulants.cxx+");
  gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithQCumulants.cxx+"); 
  gROOT->LoadMacro("AliFlowCommon/AliFittingQDistribution.cxx+");

  // Class to fill the FlowEvent without aliroot dependence
  // can be found in the directory FlowEventMakers
  gROOT->LoadMacro("FlowEventMakers/FlowEventSimpleMaker.cxx+");   

  cout << "finished loading macros!" << endl;  

  //flow event
  FlowEventSimpleMaker* fEventMaker = new FlowEventSimpleMaker(); 
  //------------------------------------------------------------------------
  */

  //------------------------------------------------------------------------
  //cuts:
  AliFlowTrackSimpleCuts* cutsInt = new AliFlowTrackSimpleCuts();
  cutsInt->SetPtMax(ptMaxInt);
  cutsInt->SetPtMin(ptMinInt);
  cutsInt->SetEtaMax(etaMaxInt);
  cutsInt->SetEtaMin(etaMinInt);
  cutsInt->SetPhiMax(phiMaxInt);
  cutsInt->SetPhiMin(phiMinInt);
  cutsInt->SetPID(PIDInt);

  AliFlowTrackSimpleCuts* cutsDiff = new AliFlowTrackSimpleCuts();
  cutsDiff->SetPtMax(ptMaxDiff);
  cutsDiff->SetPtMin(ptMinDiff);
  cutsDiff->SetEtaMax(etaMaxDiff);
  cutsDiff->SetEtaMin(etaMinDiff);
  cutsDiff->SetPhiMax(phiMaxDiff);
  cutsDiff->SetPhiMin(phiMinDiff);
  cutsDiff->SetPID(PIDDiff);

  //if the weights are used: 
  TFile *fileWithWeights = NULL;
  TList *listWithWeights = NULL;
  
  if(usePhiWeights||usePtWeights||useEtaWeights) {
    fileWithWeights = TFile::Open("weights.root","READ");
    if(fileWithWeights) {
      listWithWeights = (TList*)fileWithWeights->Get("weights");
    }
    else
      {cout << " WARNING: the file <weights.root> with weights from the previous run was not found."<<endl;
	break;
      }    
  }

  //flow methods:  
  AliFlowAnalysisWithQCumulants    *qc    = NULL;
  AliFlowAnalysisWithCumulants     *gfc   = NULL;
  AliFittingQDistribution          *fqd   = NULL;
  AliFlowAnalysisWithLeeYangZeros  *lyz1  = NULL;
  AliFlowAnalysisWithLeeYangZeros  *lyz2  = NULL;
  AliFlowAnalysisWithLYZEventPlane *lyzep = NULL;
  AliFlowAnalysisWithScalarProduct *sp    = NULL;
  AliFlowAnalysisWithMCEventPlane  *mcep  = NULL;   

  //MCEP = monte carlo event plane
  if (MCEP) {
    AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
    mcep->Init();
  }

  //QC = Q-cumulants  
  if(QC) { 
    AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
    qc->Init();
    if(listWithWeights) qc->SetWeightsList(listWithWeights);
    if(usePhiWeights) qc->SetUsePhiWeights(usePhiWeights);
    if(usePtWeights) qc->SetUsePtWeights(usePtWeights);
    if(useEtaWeights) qc->SetUseEtaWeights(useEtaWeights);
  }
  
  //GFC = Generating Function Cumulants 
  if(GFC) {
    AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
    gfc->CreateOutputObjects();
  }
  
  //FQD = Fitting q-distribution 
  if(FQD) {
    AliFittingQDistribution* fqd = new AliFittingQDistribution();
    fqd->Init();
    if(listWithWeights) fqd->SetWeightsList(listWithWeights);
    if(usePhiWeights) fqd->SetUsePhiWeights(usePhiWeights);
    if(usePtWeights) fqd->SetUsePtWeights(usePtWeights);
    if(useEtaWeights) fqd->SetUseEtaWeights(useEtaWeights);
  }

  //SP = Scalar Product 
  if(SP) {
    AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
    sp->Init();
  }

  //LYZ1 = Lee-Yang Zeroes first run
  if(LYZ1) {
    AliFlowAnalysisWithLeeYangZeros* lyz1 = new AliFlowAnalysisWithLeeYangZeros();
    lyz1->SetFirstRun(kTRUE);
    lyz1->SetUseSum(kTRUE);
    lyz1->Init();
  }

  //LYZ2 = Lee-Yang Zeroes second run
  if(LYZ2) {
    AliFlowAnalysisWithLeeYangZeros* lyz2 = new AliFlowAnalysisWithLeeYangZeros();
    // read the input file from the first run 
    TString inputFileNameLYZ2 = "outputLYZ1analysis.root" ;
    TFile* inputFileLYZ2 = new TFile(inputFileNameLYZ2.Data(),"READ");
    if(!inputFileLYZ2 || inputFileLYZ2->IsZombie()) { 
      cerr << " ERROR: NO First Run file... " << endl ; 
    }
    else { 
      TList* inputListLYZ2 = (TList*)inputFileLYZ2->Get("cobjLYZ1");  
      if (!inputListLYZ2) {cout<<"list is NULL pointer!"<<endl;}
      else {
	cout<<"LYZ2 input file/list read..."<<endl;
	lyz2->SetFirstRunList(inputListLYZ2);
	lyz2->SetFirstRun(kFALSE);
	lyz2->SetUseSum(kTRUE);
	lyz2->Init();
      }
    }
  }
  
 //LYZEP = Lee-Yang Zeroes event plane
  if(LYZEP) {
    AliFlowLYZEventPlane* ep = new AliFlowLYZEventPlane() ;
    AliFlowAnalysisWithLYZEventPlane* lyzep = new AliFlowAnalysisWithLYZEventPlane();
    // read the input file from the second lyz run 
    TString inputFileNameLYZEP = "outputLYZ2analysis.root" ;
    TFile* inputFileLYZEP = new TFile(inputFileNameLYZEP.Data(),"READ");
    if(!inputFileLYZEP || inputFileLYZEP->IsZombie()) { 
      cerr << " ERROR: NO Second Run file... " << endl ; }
    else { 
      TList* inputListLYZEP = (TList*)inputFileLYZEP->Get("cobjLYZ2");  
      if (!inputListLYZEP) {cout<<"list is NULL pointer!"<<endl;}
      else {
	cout<<"LYZEP input file/list read..."<<endl;
	ep   ->SetSecondRunList(inputListLYZEP);
	lyzep->SetSecondRunList(inputListLYZEP);
	ep   ->Init();
	lyzep->Init();
      }
    }
  }
  //------------------------------------------------------------------------
  
  
  //standard code
  Int_t fCount = 0;
  TString execDir(gSystem->pwd());
  TString targetDir(dir);
  TSystemDirectory* baseDir = new TSystemDirectory(".", dir);  
  TList* dirList 	   = baseDir->GetListOfFiles();
  if (!dirList) {
    cout << endl << "No input files in: " << targetDir.Data() << endl;
    break;
  }
  Int_t nDirs		   = dirList->GetEntries();
  cout<<endl;
  cout<<"Int_t nDirs = "<<nDirs<<endl;
  gSystem->cd(execDir);
  
  for(Int_t iDir=0;iDir<nDirs;++iDir)
    {
      TSystemFile* presentDir = (TSystemFile*)dirList->At(iDir);
      if(!presentDir || !presentDir->IsDirectory() || 
	 strcmp(presentDir->GetName(), ".") == 0 || 
	 strcmp(presentDir->GetName(), "..") == 0) 
	{
	  cout << endl; 
	  cout << "Directory (" << iDir << "): " << presentDir->GetName() << 
	    " - Skipping ... " << endl;
	  continue ;   
	}
      
      if(offset > 0) { --offset ; continue ; }
      if((aRuns > 0) && (fCount >= aRuns)) { break ; }
      
      TString presentDirName(dir); // aDataDir
      presentDirName += presentDir->GetName();
      presentDirName += "/";
      //cerr<<" presentDirName = "<<presentDirName<<endl;
      
      TString fileName = presentDirName; 
      fileName += "galice.root";
      Long_t *id, *size, *flags, *modtime;
      if(gSystem->GetPathInfo(fileName.Data(),id,size,flags,modtime)) 
	{ 
	  cout << " File : " << fileName << " does NOT exist ! - Skipping ... " 
	       << endl; 
	  continue; 
	}
      //      cout << endl ; cout << "Directory (" << iDir << "): " << presentDirName << " ... " << endl;
      
      //loop (simulations in the present dir) 
      TSystemDirectory* evtsDir = new TSystemDirectory(".", 
						       presentDirName.Data());
      TList* fileList 	    = evtsDir->GetListOfFiles();
      Int_t nFiles		    = fileList->GetEntries();
      //cout<<" Int_t nFiles = "<<nFiles<<endl;
      gSystem->cd(execDir);      
      for(Int_t iFiles=0; iFiles<nFiles; ++iFiles)
	{
	  TSystemFile* presentFile = (TSystemFile*) fileList->At(iFiles);
	  TString presentFileName(presentDirName);
	  presentFileName += presentFile->GetName();
	  
	  if(!(presentFileName.Contains("Kinematics") && 
	       presentFileName.Contains("root"))) { continue ; }
	  
	  //cout << " found: " << presentFileName.Data() << endl; 
	  
	  TFile* kineFile = new TFile(presentFileName.Data(), "READ"); 
	  // kineFile->ls();
	  Int_t nEvts = kineFile->GetNkeys() ; 
	  //cout << "  . found: " << nEvts << " KineTree(s) in " << presentFileName.Data() << endl;
	  TList* kineEventsList = (TList*)kineFile->GetListOfKeys(); 
	  TTree* kTree;
	  TIter next(kineEventsList); 
	  TKey* key;
	  
	  //loop over the events
	  while( key=(TKey *)next() ) 
	    {
	      TDirectory* tDir = (TDirectory*)key->ReadObj();
	      if(!tDir) break;
	      
	      TString evtDir(tDir->GetName()); 
	      //cout << "  . . found: " << tDir->GetName() << endl;
	      
	      kTree = (TTree *)tDir->Get("TreeK");
	      if(!kTree) break;
	      
	      Int_t nPart = kTree->GetEntries();
	      //cout << "  . . . kTree " << fCount << " has " << nPart << " particles " << endl; 
	      
	      //-----------------------------------------------------------
	      //fill and save the flow event	      
	      AliFlowEventSimple *fEvent = fEventMaker->FillTracks(kTree, cutsInt, cutsDiff); 
	                    
	      //pass the flow event to flow methods for analysis 
	      Double_t fEP = 0.; // temporary number need true value
	      if(MCEP) mcep->Make(fEvent,fEP);  //fix fEP
	      if(QC) qc->Make(fEvent);
	      if(GFC) gfc->Make(fEvent);
	      if(FQD) fqd->Make(fEvent);
	      if(LYZ1) lyz1->Make(fEvent);
	      if(LYZ2) lyz2->Make(fEvent);
	      if(LYZEP) lyzep->Make(fEvent,ep);
	      if(SP) sp->Make(fEvent);
	      //-----------------------------------------------------------
	      
	      fCount++;
	      cout << "# " << fCount << " events processed" << endl;
	      delete kTree;
	      delete fEvent;
	    }
	  delete kineFile ;
	}
      delete evtsDir ;
    }

  //--------------------------------------------------------------
  //calculating and storing the final results of flow analysis
  if(MCEP) {mcep->Finish(); mcep->WriteHistograms("outputMCEPanalysis.root");}
  if(QC) {qc->Finish(); qc->WriteHistograms("outputQCanalysis.root");}
  if(GFC) {gfc->Finish(); gfc->WriteHistograms("outputGFCanalysis.root");}
  if(FQD) {fqd->Finish(); fqd->WriteHistograms("outputFQDanalysis.root");}
  if(LYZ1) {lyz1->Finish(); lyz1->WriteHistograms("outputLYZ1analysis.root");}
  if(LYZ2) {lyz2->Finish(); lyz2->WriteHistograms("outputLYZ2analysis.root");}
  if(LYZEP) {lyzep->Finish(); lyzep->WriteHistograms("outputLYZEPanalysis.root");}
  if(SP) {sp->Finish(); sp->WriteHistograms("outputSPanalysis.root");}
  //--------------------------------------------------------------
  
  cout << endl;
  cout << " Finished ... " << endl;
  cout << endl;
  
  timer.Stop();
  cout << endl;
  timer.Print();
}

void SetupPar(char* pararchivename)
{
  //Load par files, create analysis libraries
  //For testing, if par file already decompressed and modified
  //classes then do not decompress.
 
  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ; 
  TString parpar(Form("%s.par", pararchivename)) ; 
  if ( gSystem->AccessPathName(parpar.Data()) ) {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_ROOT")) ;
    TString processline(Form(".! make %s", parpar.Data())) ; 
    gROOT->ProcessLine(processline.Data()) ;
    gSystem->ChangeDirectory(cdir) ; 
    processline = Form(".! mv /tmp/%s .", parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()) ;
  } 
  if ( gSystem->AccessPathName(pararchivename) ) {  
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }
  
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);
  
  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    
    if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
      Error("runProcess","Cannot Build the PAR Archive! - Abort!");
      return -1;
    }
  }
  // check for SETUP.C and execute
  if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
    printf("*******************************\n");
    printf("*** Setup PAR archive       ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C");
  }
  
  gSystem->ChangeDirectory(ocwd.Data());
  printf("Current dir: %s\n", ocwd.Data());
}

