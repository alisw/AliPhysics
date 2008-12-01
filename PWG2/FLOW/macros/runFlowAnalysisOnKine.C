#include "TStopwatch.h"
#include "TObjArray"
#include "Riostream.h"

//RUN SETTINGS
//flow analysis method can be: (set to kTRUE or kFALSE)
Bool_t SP    = kFALSE;
Bool_t LYZ1  = kTRUE;
Bool_t LYZ2  = kFALSE;
Bool_t LYZEP = kFALSE;
Bool_t GFC   = kTRUE;
Bool_t QC    = kTRUE;
Bool_t FQD   = kTRUE;
Bool_t MCEP  = kFALSE;

Int_t offset = 0;

int runFlowAnalysisOnKine(Int_t aRuns = 144, const char* 
			  dir="/data/alice1/kolk/KineOnly3/")

{
 TStopwatch timer;
 timer.Start();

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
 //flow event
 AliFlowEventSimpleMaker* fEventMaker = new AliFlowEventSimpleMaker(); 
 
 AliFlowAnalysisWithQCumulants   *qc   = NULL;
 AliFlowAnalysisWithCumulants    *gfc  = NULL;
 AliFittingQDistribution         *fqd  = NULL;
 AliFlowAnalysisWithLeeYangZeros *lyz1 = NULL;

 //flow methods:
 //QC = Q-cumulants  
 if(QC)
 { 
  AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
  qc->CreateOutputObjects();
 }

 //GFC = Generating Function Cumulants 
 AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
 if(GFC)
 {
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
  lyz1->Init();
 }

 //------------------------------------------------------------------------

 //standard code
 Int_t fCount = 0;
 TString execDir(gSystem->pwd());
 TSystemDirectory* baseDir = new TSystemDirectory(".", dir);  
 TList* dirList 	   = baseDir->GetListOfFiles();
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

    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(kTree);
     
    //pass the flow event to flow methods for analysis  
    //QC
    if(QC)
    {  
    cout <<"what the hell" << endl;       
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
    //-----------------------------------------------------------
    
    fCount++;
    delete kTree;
   }
   delete kineFile ;
  }
  delete evtsDir ;
 }
 
 //--------------------------------------------------------------
 //calculating and storing the final results of flow analysis
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
  //TString *outputFileNameFQD = new TString("outputFQDanalysis.root");
  //fqd->WriteHistograms(outputFileNameFQD);
  //delete outputFileNameFQD;
 }
 //--------------------------------------------------------------

 cout << endl;
 cout << " Finished ... " << endl;
 cout << endl;
 
 timer.Stop();
 cout << endl;
 timer.Print();
}
