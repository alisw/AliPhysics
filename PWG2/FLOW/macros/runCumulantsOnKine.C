#include "TStopwatch.h"
#include "TObjArray"
#include "Riostream.h"

Int_t offset = 0 ;

//int runCumulantsOnKine(Int_t aRuns = 144, Bool_t fFirstRunLYZ = kTRUE, Bool_t fDouble = kFALSE, const char * dir = "/users/alice/KineOnly3/")
//int runCumulantsOnKine(Int_t aRuns = 10, Bool_t fFirstRunLYZ = kTRUE, Bool_t fDouble = kFALSE, const char * dir = "/data/alice2/LHyquid3_rot/")
//int runCumulantsOnKine(Int_t aRuns = 144, const char * dir = "/Users/snelling/alice_data/KineOnly3/")
int runCumulantsOnKine(Int_t aRuns = 100, const char * dir = "/Users/snelling/alice_data/TherminatorFIX/")
//int runCumulantsOnKine(Int_t aRuns = 200, Bool_t fFirstRunLYZ = kTRUE, Bool_t fDouble = kFALSE, const char * dir = "/data/alice2/abDeleteMeASAP/")
{
  TStopwatch timer;
  timer.Start();

  cout<<" BEGIN ANALYSIS "<<endl;
  gSystem->AddIncludePath("-I$ALICE_ROOT/include") ;
  gSystem->AddIncludePath("-I$ROOTSYS/include") ;
  gROOT->LoadMacro("AliFlowVector.cxx+");
  gROOT->LoadMacro("AliFlowCommonConstants.cxx+");
  gROOT->LoadMacro("AliFlowCumuConstants.cxx+");
  gROOT->LoadMacro("AliFlowTrackSimple.cxx+");
  gROOT->LoadMacro("AliFlowEventSimple.cxx+");
  gROOT->LoadMacro("AliFlowEventSimpleMaker.cxx+");
  gROOT->LoadMacro("AliFlowCommonHist.cxx+");
  gROOT->LoadMacro("AliFlowCommonHistResults.cxx+");
  gROOT->LoadMacro("AliFlowAnalysisWithCumulants.cxx+");

  cout<<" loaded macros "<<endl;


  AliFlowAnalysisWithCumulants* aFlowAnalysisWithCumulants = new AliFlowAnalysisWithCumulants();//ab

  AliFlowEventSimpleMaker* fEventMaker = new AliFlowEventSimpleMaker(); 

  aFlowAnalysisWithCumulants->CreateOutputObjects();

  // standard code
  
  Int_t fCount = 0 ;

  TString execDir(gSystem->pwd());
  TSystemDirectory* baseDir = new TSystemDirectory(".", dir) ;  // aDataDir   // TSystemDirectory(const char* dirname, const char* path)
  TList* dirList 	   = baseDir->GetListOfFiles();
  Int_t nDirs		   = dirList->GetEntries();
  cout<<" Int_t nDirs = "<<nDirs<<endl;
  gSystem->cd(execDir);

  for(Int_t iDir=0; iDir<nDirs; ++iDir)
    {
      TSystemFile* presentDir = (TSystemFile*)dirList->At(iDir) ;
      if(!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0) 
	{
	  cout << endl ; 
	  cout << "Directory (" << iDir << "):  " << presentDir->GetName() << " - Skipping ... " << endl ;
	  continue ;   
	}
      if(offset > 0)  { --offset ; continue ; }
      if((aRuns > 0) && (fCount >= aRuns)) { break ; }
 
      TString presentDirName(dir) ; // aDataDir
      presentDirName += presentDir->GetName();
      presentDirName += "/";
      //cerr<<" presentDirName = "<<presentDirName<<endl;

      TString fileName = presentDirName ; 
      fileName += "galice.root" ;
      Long_t *id, *size, *flags, *modtime ;
      if(gSystem->GetPathInfo(fileName.Data(),id,size,flags,modtime)) 
	{ 
	  cout << " File : " << fileName << " does NOT exist ! - Skipping ... " << endl ; 
	  continue ; 
	}
      cout << endl ; cout << "Directory (" << iDir << "):  " << presentDirName << "  ... " << endl ;

      // loop (simulations in the present dir) //
      TSystemDirectory* evtsDir = new TSystemDirectory(".", presentDirName.Data());
      TList* fileList 	    = evtsDir->GetListOfFiles();
      Int_t nFiles		    = fileList->GetEntries();
      cout<<" Int_t nFiles = "<<nFiles<<endl;
      gSystem->cd(execDir);

      for(Int_t iFiles=0; iFiles<nFiles; ++iFiles)
	{
	  TSystemFile* presentFile = (TSystemFile*) fileList->At(iFiles);

	  TString presentFileName(presentDirName);
	  presentFileName += presentFile->GetName();

	  if(!(presentFileName.Contains("Kinematics") && presentFileName.Contains("root"))) { continue ; }

	  cout << " found: " << presentFileName.Data() << endl ; 
  
	  TFile* kineFile = new TFile(presentFileName.Data(), "READ") ; 
	  // kineFile->ls() ;
	  Int_t nEvts = kineFile->GetNkeys() ; 
	  cout << "  . found: " << nEvts << " KineTree(s) in " << presentFileName.Data() << endl ;
	  TList* kineEventsList = (TList*)kineFile->GetListOfKeys() ; 
	  TTree* kTree ;
	  TIter next(kineEventsList); 
	  TKey* key ;

	  //end common code
	 
	  // Loop over the events
	  while( key=(TKey *)next() ) 
	    {
	      TDirectory* tDir = (TDirectory*)key->ReadObj() ;
	      if(!tDir) break;
 
	      TString evtDir(tDir->GetName()) ; 
	      cout << "  . . found: " << tDir->GetName() << endl ;

	      kTree = (TTree *)tDir->Get("TreeK");
	      if(!kTree) break;

	      Int_t nPart = kTree->GetEntries() ;
	      cout << "  . . . kTree " << fCount << " has " << nPart << " particles " << endl ;
	      
	      // fill and save the flow event

	      AliFlowEventSimple* fEvent = fEventMaker->FillTracks(kTree); 
              aFlowAnalysisWithCumulants->Exec(fEvent); cout<<"Cumu analysis"<<endl;
              
	      fCount++ ;
    
	      delete kTree;
	    }
	  delete kineFile ;
	}
      delete evtsDir ;
    }

  aFlowAnalysisWithCumulants->Terminate(fCount);

  cout <<  endl ;
  cout << " Finished ... " << endl ;
  timer.Stop() ;
  cout << endl ;
  timer.Print() ;
}
