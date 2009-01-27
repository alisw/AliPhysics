#include "TStopwatch.h"
#include "TObjArray"
#include "Riostream.h"

//RUN SETTINGS
//flow analysis method can be: (set to kTRUE or kFALSE)
Bool_t SP    = kTRUE;
Bool_t LYZ1  = kTRUE;
Bool_t LYZ2  = kFALSE;  
Bool_t LYZEP = kFALSE; 
Bool_t GFC   = kTRUE;
Bool_t QC    = kTRUE;
Bool_t FQD   = kTRUE;
Bool_t MCEP  = kFALSE; //does not work yet 24/12/08

Int_t offset = 0;

//CUT SETTINGS
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

int runFlowAnalysisOnKine(Int_t aRuns = 100, const char* 
			  //			  dir="/data/alice1/kolk/KineOnly3/")
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
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ROOTSYS/include");
    
  /* 
     gROOT->LoadMacro("AliFlowCommonConstants.cxx+");
     gROOT->LoadMacro("AliFlowLYZConstants.cxx+");
     gROOT->LoadMacro("AliFlowCumuConstants.cxx+");
     
     gROOT->LoadMacro("AliFlowEventSimple.cxx+");
     gROOT->LoadMacro("AliFlowTrackSimple.cxx+");
     gROOT->LoadMacro("AliFlowCommonHist.cxx+");
     gROOT->LoadMacro("AliFlowCommonHistResults.cxx+");
     gROOT->LoadMacro("AliFlowLYZHist1.cxx+");
     gROOT->LoadMacro("AliFlowLYZHist2.cxx+");
     gROOT->LoadMacro("AliFlowVector.cxx+");
     gROOT->LoadMacro("AliFlowLYZEventPlane.cxx+");
     gROOT->LoadMacro("AliFlowEventSimpleMaker.cxx+"); 
     gROOT->LoadMacro("AliFlowAnalysisWithMCEventPlane.cxx+"); 
     gROOT->LoadMacro("AliFlowAnalysisWithScalarProduct.cxx+");
     gROOT->LoadMacro("AliFlowAnalysisWithLYZEventPlane.cxx+");
     gROOT->LoadMacro("AliFlowAnalysisWithLeeYangZeros.cxx+");
     gROOT->LoadMacro("AliFlowAnalysisWithCumulants.cxx+");
     gROOT->LoadMacro("AliFlowAnalysisWithQCumulants.cxx+"); 
     gROOT->LoadMacro("AliAnalysisTaskScalarProduct.cxx+");
     gROOT->LoadMacro("AliAnalysisTaskMCEventPlane.cxx+");
     gROOT->LoadMacro("AliAnalysisTaskLYZEventPlane.cxx+");
     gROOT->LoadMacro("AliAnalysisTaskLeeYangZeros.cxx+");
     gROOT->LoadMacro("AliAnalysisTaskCumulants.cxx+"); 
     gROOT->LoadMacro("AliAnalysisTaskQCumulants.cxx+");
     gROOT->LoadMacro("AliCumulantsFunctions.cxx+");
     gROOT->LoadMacro("AliQCumulantsFunctions.cxx+");
     gROOT->LoadMacro("AliAnalysisTaskFittingQDistribution.cxx+");
     gROOT->LoadMacro("AliFittingFunctionsForQDistribution.cxx+");
     gROOT->LoadMacro("AliFittingQDistribution.cxx+");
  */
  
  //load needed libraries
  gSystem->Load("libTree.so");
  
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libPWG2flow.so");
  cerr<<"libPWG2flow.so loaded ..."<<endl;
  cout<<endl; 

  //------------------------------------------------------------------------
  //cuts
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

  //------------------------------------------------------------------------
  //flow event
  AliFlowEventSimpleMaker* fEventMaker = new AliFlowEventSimpleMaker(); 
  
  AliFlowAnalysisWithQCumulants    *qc    = NULL;
  AliFlowAnalysisWithCumulants     *gfc   = NULL;
  AliFittingQDistribution          *fqd   = NULL;
  AliFlowAnalysisWithLeeYangZeros  *lyz1  = NULL;
  AliFlowAnalysisWithLeeYangZeros  *lyz2  = NULL;
  AliFlowAnalysisWithLYZEventPlane *lyzep = NULL;
  AliFlowAnalysisWithScalarProduct *sp    = NULL;
  AliFlowAnalysisWithMCEventPlane  *mcep  = NULL;
   
  //flow methods:
  //MCEP = monte carlo event plane
  AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
  if (MCEP)
    {
      mcep->Init();
    }

  //QC = Q-cumulants  
  if(QC)
    { 
      AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
      qc->CreateOutputObjects();
    }
  
  //GFC = Generating Function Cumulants 
  if(GFC)
    {
      AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
      gfc->CreateOutputObjects();
    }
  
  //FQD = Fitting q-distribution 
  AliFittingQDistribution* fqd = new AliFittingQDistribution();
  if(FQD)
    {
      fqd->CreateOutputObjects();
    }
  
  //LYZ1 = Lee-Yang Zeroes first run
  AliFlowAnalysisWithLeeYangZeros* lyz1 = new AliFlowAnalysisWithLeeYangZeros();
  if(LYZ1)
    {
      lyz1->SetFirstRun(kTRUE);
      lyz1->SetUseSum(kTRUE);
      lyz1->Init();
    }

  //LYZ2 = Lee-Yang Zeroes second run
  AliFlowAnalysisWithLeeYangZeros* lyz2 = new AliFlowAnalysisWithLeeYangZeros();
  if(LYZ2)
    {
      // read the input file from the first run 
      TString inputFileNameLYZ2 = "outputLYZ1analysis.root" ;
      TFile* inputFileLYZ2 = new TFile(inputFileNameLYZ2.Data(),"READ");
      if(!inputFileLYZ2 || inputFileLYZ2->IsZombie()) { 
	cerr << " ERROR: NO First Run file... " << endl ; }
      else { 
	TList* inputListLYZ2 = (TList*)inputFileLYZ2->Get("cobjLYZ1");  
	if (!inputListLYZ2) {cout<<"list is NULL pointer!"<<endl;
	}
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
  AliFlowLYZEventPlane* ep = new AliFlowLYZEventPlane() ;
  AliFlowAnalysisWithLYZEventPlane* lyzep = new AliFlowAnalysisWithLYZEventPlane();
  if(LYZEP)
    {
      // read the input file from the second lyz run 
      TString inputFileNameLYZEP = "outputLYZ2analysis.root" ;
      TFile* inputFileLYZEP = new TFile(inputFileNameLYZEP.Data(),"READ");
      if(!inputFileLYZEP || inputFileLYZEP->IsZombie()) { 
	cerr << " ERROR: NO Second Run file... " << endl ; }
      else { 
	TList* inputListLYZEP = (TList*)inputFileLYZEP->Get("cobjLYZ2");  
	if (!inputListLYZEP) {cout<<"list is NULL pointer!"<<endl;
	}
	else {
	  cout<<"LYZEP input file/list read..."<<endl;
	  ep   ->SetSecondRunList(inputListLYZEP);
	  lyzep->SetSecondRunList(inputListLYZEP);
	  ep   ->Init();
	  lyzep->Init();
	}
      }
    }
   
  //SP = Scalar Product 
  AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
  if(SP)
    {
      sp->Init();
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
      cout << endl ; cout << "Directory (" << iDir << "): " << presentDirName 
	   << " ... " << endl;
      
      //loop (simulations in the present dir) 
      TSystemDirectory* evtsDir = new TSystemDirectory(".", 
						       presentDirName.Data());
      TList* fileList 	    = evtsDir->GetListOfFiles();
      Int_t nFiles		    = fileList->GetEntries();
      cout<<" Int_t nFiles = "<<nFiles<<endl;
      gSystem->cd(execDir);
      
      for(Int_t iFiles=0; iFiles<nFiles; ++iFiles)
	{
	  TSystemFile* presentFile = (TSystemFile*) fileList->At(iFiles);
	  TString presentFileName(presentDirName);
	  presentFileName += presentFile->GetName();
	  
	  if(!(presentFileName.Contains("Kinematics") && 
	       presentFileName.Contains("root"))) { continue ; }
	  
	  cout << " found: " << presentFileName.Data() << endl; 
	  
	  TFile* kineFile = new TFile(presentFileName.Data(), "READ"); 
	  // kineFile->ls();
	  Int_t nEvts = kineFile->GetNkeys() ; 
	  cout << "  . found: " << nEvts << " KineTree(s) in " << 
	    presentFileName.Data() << endl;
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
	      cout << "  . . found: " << tDir->GetName() << endl;
	      
	      kTree = (TTree *)tDir->Get("TreeK");
	      if(!kTree) break;
	      
	      Int_t nPart = kTree->GetEntries();
	      cout << "  . . . kTree " << fCount << " has " << nPart << 
		" particles " << endl; 
	      
	      
	      
	      //-----------------------------------------------------------
	      //fill and save the flow event
	      
	      AliFlowEventSimple* fEvent = fEventMaker->FillTracks(kTree, cutsInt, cutsDiff);
	      
	      //pass the flow event to flow methods for analysis 
	      //MCEP
	      if (MCEP)
		{
		  //mcep->Make(fEvent,fEP);  //fix fEP
		  cout<<"  --> MCEP analysis..."<<endl;
		}
 	      //QC
	      if(QC)
		{  
		  qc->Make(fEvent);
		  cout<<"  --> QC analysis..."<<endl;
		}
	      //GFC
	      if(GFC)
		{
		  gfc->Make(fEvent);
		  cout<<"  --> GFC analysis..."<<endl;
		}
	      //FQD 
	      if(FQD)
		{
		  fqd->Make(fEvent);
		  cout<<"  --> FQD analysis..."<<endl;
		}
	      //LYZ1
	      if(LYZ1)
		{
		  lyz1->Make(fEvent);
		  cout<<"  --> LYZ1 analysis..."<<endl;
		}
	      //LYZ2
	      if(LYZ2)
		{
		  lyz2->Make(fEvent);
		  cout<<"  --> LYZ2 analysis..."<<endl;
		}
	      //LYZEP
	      if(LYZEP)
		{
		  lyzep->Make(fEvent,ep);
		  cout<<"  --> LYZEP analysis..."<<endl;
		}
	      //SP
	      if(SP)
		{
		  sp->Make(fEvent);
		  cout<<"  --> SP analysis..."<<endl;
		}



	      //-----------------------------------------------------------
	      
	      fCount++;
	      delete kTree;
	      delete fEvent;
	    }
	  delete kineFile ;
	}
      delete evtsDir ;
    }
  
  //--------------------------------------------------------------
  //calculating and storing the final results of flow analysis
  //MCEP
  if (MCEP)
    {
      mcep->Finish();
      TString *outputFileNameMCEP = new TString("outputMCEPanalysis.root");
      //mcep->WriteHistograms(outputFileNameMCEP); //add this method to MCEP classes
      delete outputFileNameMCEP;
    }
  //QC
  if(QC)
    {
      qc->Finish();
      TString *outputFileNameQC = new TString("outputQCanalysis.root");
      qc->WriteHistograms(outputFileNameQC);
      delete outputFileNameQC;
    }
  //GFC
  if(GFC)
    {
      gfc->Finish();
      TString *outputFileNameGFC = new TString("outputGFCanalysis.root");
      gfc->WriteHistograms(outputFileNameGFC);
      delete outputFileNameGFC;
    }
  //FQD
  if(FQD)
    {
      fqd->Finish();
      TString *outputFileNameFQD = new TString("outputFQDanalysis.root");
      fqd->WriteHistograms(outputFileNameFQD);
      delete outputFileNameFQD;
    }
  //LYZ1
  if(LYZ1)
    {
      lyz1->Finish();
      TString *outputFileNameLYZ1 = new TString("outputLYZ1analysis.root");
      lyz1->WriteHistograms(outputFileNameLYZ1);
      delete outputFileNameLYZ1;
    }
  //LYZ2
  if(LYZ2)
    {
      lyz2->Finish();
      TString *outputFileNameLYZ2 = new TString("outputLYZ2analysis.root");
      lyz2->WriteHistograms(outputFileNameLYZ2);
      delete outputFileNameLYZ2;
    }
  //LYZEP
  if(LYZEP)
    {
      lyzep->Finish();
      TString *outputFileNameLYZEP = new TString("outputLYZEPanalysis.root");
      lyzep->WriteHistograms(outputFileNameLYZEP);
      delete outputFileNameLYZEP;
    }
  //SP
  if(SP)
    {
      sp->Finish();
      TString *outputFileNameSP = new TString("outputSPanalysis.root");
      sp->WriteHistograms(outputFileNameSP);
      delete outputFileNameSP;
    }



  //--------------------------------------------------------------
  
  cout << endl;
  cout << " Finished ... " << endl;
  cout << endl;
  
  timer.Stop();
  cout << endl;
  timer.Print();
}
