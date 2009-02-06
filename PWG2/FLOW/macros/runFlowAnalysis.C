#include "TStopwatch.h"
#include "TObjArray"
#include "Riostream.h"
#include "TFile.h"

//----------------------------------------------------------
//RUN SETTINGS
//flow analysis method can be: (set to kTRUE or kFALSE)
Bool_t SP    = kFALSE;
Bool_t LYZ1  = kTRUE;
Bool_t LYZ2  = kFALSE;  
Bool_t LYZEP = kFALSE; 
Bool_t GFC   = kTRUE;
Bool_t QC    = kTRUE;
Bool_t FQD   = kFALSE;
Bool_t MCEP  = kFALSE; //does not work yet 24/12/08
//----------------------------------------------------------

//----------------------------------------------------------
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
//----------------------------------------------------------

//----------------------------------------------------------
//WEIGHTS SETTINGS 
//to use or not to use the weights - that is a question!
Bool_t useWeightsPhi = kFALSE;//Phi
Bool_t useWeightsPt  = kFALSE;//v'(pt)
Bool_t useWeightsEta = kFALSE;//v'(eta)
//----------------------------------------------------------

Int_t offset = 0;
                                          
int runFlowAnalysis(Int_t aRuns = 44, const char* 
			 			  dir="/data/alice1/kolk/KineOnly3/")
			  //dir="/Users/snelling/alice_data/KineOnly3/")
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

  // In AliRoot
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libPWG2flow.so");
  cerr<<"libPWG2flow.so loaded ..."<<endl;
  cout<<endl; 

  //open the file with the weights: 
  TFile *file = NULL;
  if(useWeightsPhi||useWeightsPt||useWeightsEta)
  {
   file = TFile::Open("weightsForTheSecondRun.root", "READ");
  }
  
  // flow event in AliRoot
  AliFlowEventSimpleMaker* fEventMaker = new AliFlowEventSimpleMaker();

  // if weights are being used initialize the weight's histos in AliFlowEvenSimpleMaker:
  fEventMaker->SetUseWeightsPhi(useWeightsPhi);
  fEventMaker->SetUseWeightsPt(useWeightsPt);
  fEventMaker->SetUseWeightsEta(useWeightsEta);
  if(useWeightsPhi||useWeightsPt||useWeightsEta)
  {
   fEventMaker->Init(file);
  }
  
  // In root

  /*    
  // constants  
  gROOT->LoadMacro("code/AliFlowCommonConstants.cxx+");
  gROOT->LoadMacro("code/AliFlowLYZConstants.cxx+");
  gROOT->LoadMacro("code/AliFlowCumuConstants.cxx+");

  // flow event
  gROOT->LoadMacro("code/AliFlowVector.cxx+"); 
  gROOT->LoadMacro("code/AliFlowTrackSimple.cxx+");    
  gROOT->LoadMacro("code/AliFlowEventSimple.cxx+");

  // cuts
  gROOT->LoadMacro("code/AliFlowTrackSimpleCuts.cxx+");    

  // output histosgrams
  gROOT->LoadMacro("code/AliFlowCommonHist.cxx+");
  gROOT->LoadMacro("code/AliFlowCommonHistResults.cxx+");
  gROOT->LoadMacro("code/AliFlowLYZHist1.cxx+");
  gROOT->LoadMacro("code/AliFlowLYZHist2.cxx+");

  // functions needed for various methods
  gROOT->LoadMacro("code/AliCumulantsFunctions.cxx+");
  gROOT->LoadMacro("code/AliQCumulantsFunctions.cxx+");
  gROOT->LoadMacro("code/AliFittingFunctionsForQDistribution.cxx+");
  gROOT->LoadMacro("code/AliFlowLYZEventPlane.cxx+");

  // Flow Analysis code for various methods
  gROOT->LoadMacro("code/AliFlowAnalysisWithMCEventPlane.cxx+"); 
  gROOT->LoadMacro("code/AliFlowAnalysisWithScalarProduct.cxx+");
  gROOT->LoadMacro("code/AliFlowAnalysisWithLYZEventPlane.cxx+");
  gROOT->LoadMacro("code/AliFlowAnalysisWithLeeYangZeros.cxx+");
  gROOT->LoadMacro("code/AliFlowAnalysisWithCumulants.cxx+");
  gROOT->LoadMacro("code/AliFlowAnalysisWithQCumulants.cxx+"); 
  gROOT->LoadMacro("code/AliFittingQDistribution.cxx+");

  // Class to fill the FlowEvent without aliroot dependence
  // can be found in the directory other
  gROOT->LoadMacro("code/FlowEventSimpleMaker.cxx+");   

  cout << "finished loading macros!" << endl;  

  //flow event
  FlowEventSimpleMaker* fEventMaker = new FlowEventSimpleMaker(); 
  //------------------------------------------------------------------------
  */

  //load needed libraries
  gSystem->Load("libTree.so");

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
  if (MCEP) {
    AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
    mcep->Init();
  }

  //QC = Q-cumulants  
  if(QC) { 
    AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
    qc->CreateOutputObjects();
  }
  
  //GFC = Generating Function Cumulants 
  if(GFC) {
    AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
    gfc->CreateOutputObjects();
  }
  
  //FQD = Fitting q-distribution 
  if(FQD) {
    AliFittingQDistribution* fqd = new AliFittingQDistribution();
    fqd->CreateOutputObjects();
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
