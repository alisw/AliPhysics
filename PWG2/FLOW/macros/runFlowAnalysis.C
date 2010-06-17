#include <Riostream.h>
#include <TStopwatch.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TTimeStamp.h>

#include <AliFlowCommon/AliFlowCommonConstants.h>
#include <AliFlowCommon/AliFlowLYZConstants.h>
#include <AliFlowCommon/AliFlowCumuConstants.h>
#include <AliFlowCommon/AliFlowVector.h>
#include <AliFlowCommon/AliFlowTrackSimple.h>
#include <AliFlowCommon/AliFlowEvent.h>
#include <AliFlowCommon/AliFlowEventSimple.h>
#include <AliFlowCommon/AliFlowTrackSimpleCuts.h>
#include <AliFlowCommon/AliFlowCommonHist.h>
#include <AliFlowCommon/AliFlowCommonHistResults.h>
#include <AliFlowCommon/AliFlowLYZHist1.h>
#include <AliFlowCommon/AliFlowLYZHist2.h>
#include <AliFlowCommon/AliCumulantsFunctions.h>
#include <AliFlowCommon/AliFlowLYZEventPlane.h>
#include <AliFlowCommon/AliFlowAnalysisWithMCEventPlane.h>
#include <AliFlowCommon/AliFlowAnalysisWithScalarProduct.h>
#include <AliFlowCommon/AliFlowAnalysisWithLYZEventPlane.h>
#include <AliFlowCommon/AliFlowAnalysisWithLeeYangZeros.h>
#include <AliFlowCommon/AliFlowAnalysisWithCumulants.h>
#include <AliFlowCommon/AliFlowAnalysisWithQCumulants.h>
#include <AliFlowCommon/AliFlowAnalysisWithFittingQDistribution.h>
#include <AliFlowCommon/AliFlowAnalysisWithMixedHarmonics.h>
#include <AliFlowCommon/AliFlowAnalysisWithNestedLoops.h>

//--------------------------------------------------------------------------------------
// Run flow analysis on local data with custom FlowEvent maker
// RUN SETTINGS
//flow analysis method can be: (set to kTRUE or kFALSE)
Bool_t SP       = kTRUE;
Bool_t LYZ1SUM  = kTRUE;
Bool_t LYZ1PROD = kTRUE;
Bool_t LYZ2SUM  = kFALSE; 
Bool_t LYZ2PROD = kFALSE;
Bool_t LYZEP    = kFALSE; 
Bool_t GFC      = kTRUE;
Bool_t QC       = kTRUE;
Bool_t FQD      = kTRUE;
Bool_t MH       = kTRUE; 
Bool_t NL       = kFALSE; 
Bool_t MCEP     = kFALSE; //does not work yet 24/12/08
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
/*
//--------------------------------------------------------------------------------------
// FLOW SETTINGS (R.Rietkerk)
Int_t nLoops=1; 		// Number of times to use the same particle (nonflow).
Double_t xEllipticFlowValue=0.1;// Add Elliptic Flow. Must be in range [0,1].
Int_t nMultiplicityOfEvent=500; // Set Average Multiplicity.
Double_t xSigmaFlow=0.00;	// Add Elliptic Flow. Must be in range [0,1].
Int_t nSigmaMult=50;            // Set Average Multiplicity.
//--------------------------------------------------------------------------------------
*/

enum anaModes {mLocal,mLocalSource,mLocalPAR,};
//mLocal: Analyze data on your computer using aliroot
//mLocalPAR: Analyze data on your computer using root + PAR files
//mLocalSource: Analyze data on your computer using root + source files

void LoadLibraries(const anaModes mode);

Int_t offset = 0;
                                          
int runFlowAnalysis(const anaModes mode=mLocal, Int_t aRuns = 100, const char* 
		    dir="/data/alice1/kolk/KineOnly3/")
		    //		    dir="/Users/snelling/alice_data/KineOnly3/")
		    //		    dir="/Users/snelling/alice_data/stoomboot/5b/")
{
  TStopwatch timer;
  timer.Start();
  
  if (LYZ1SUM && LYZ2SUM) {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; exit(1); }
  if (LYZ1PROD && LYZ2PROD) {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; exit(1); }
  if (LYZ2SUM && LYZEP) {cout<<"WARNING: you cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(1); }
  if (LYZ1SUM && LYZEP) {cout<<"WARNING: you cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(1); }


  cout<<endl;
  cout<<" ---- BEGIN ANALYSIS ---- "<<endl;
  cout<<endl;
  
  LoadLibraries(mode);

  TRandom3 random3Temp; //init for manual settings (R.Rietkerk)
  TTimeStamp dt;
  Int_t sseed = dt.GetNanoSec()/1000;
  random3Temp.SetSeed(sseed);

  if (mode == mLocal || mode == mLocalPAR) {
  }
  else if (mode == mLocalSource) {
  }
  else{
    cout << "No supported running mode selected!" << endl;
    break;
  }


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
  AliFlowAnalysisWithQCumulants *qc = NULL;
  AliFlowAnalysisWithCumulants *gfc = NULL;
  AliFlowAnalysisWithFittingQDistribution *fqd = NULL;
  AliFlowAnalysisWithLeeYangZeros *lyz1sum = NULL;
  AliFlowAnalysisWithLeeYangZeros *lyz1prod = NULL;
  AliFlowAnalysisWithLeeYangZeros *lyz2sum = NULL;
  AliFlowAnalysisWithLeeYangZeros *lyz2prod = NULL;
  AliFlowAnalysisWithLYZEventPlane *lyzep = NULL;
  AliFlowAnalysisWithScalarProduct *sp = NULL;
  AliFlowAnalysisWithMCEventPlane *mcep = NULL;     
  AliFlowAnalysisWithMixedHarmonics *mh = NULL;
  AliFlowAnalysisWithNestedLoops *nl = NULL;

  //MCEP = monte carlo event plane
  if (MCEP) {
    AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
    mcep->Init();
  }

  //QC = Q-cumulants  
  if(QC) { 
    AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
    if(listWithWeights) qc->SetWeightsList(listWithWeights);
    if(usePhiWeights) qc->SetUsePhiWeights(usePhiWeights);
    if(usePtWeights) qc->SetUsePtWeights(usePtWeights);
    if(useEtaWeights) qc->SetUseEtaWeights(useEtaWeights);
    qc->Init();
  }
  
  //GFC = Generating Function Cumulants 
  if(GFC) {
    AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
    if(listWithWeights) gfc->SetWeightsList(listWithWeights);
    if(usePhiWeights) gfc->SetUsePhiWeights(usePhiWeights);
    if(usePtWeights) gfc->SetUsePtWeights(usePtWeights);
    if(useEtaWeights) gfc->SetUseEtaWeights(useEtaWeights);
    gfc->Init();
  }
  
  //FQD = Fitting q-distribution 
  if(FQD) {
    AliFlowAnalysisWithFittingQDistribution* fqd = new AliFlowAnalysisWithFittingQDistribution();
    if(listWithWeights) fqd->SetWeightsList(listWithWeights);
    if(usePhiWeights) fqd->SetUsePhiWeights(usePhiWeights);
    fqd->Init();
  }

  //SP = Scalar Product 
  if(SP) {
    AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
    if(usePhiWeights) sp->SetUsePhiWeights(usePhiWeights);
    sp->Init();
  }

  //LYZ1 = Lee-Yang Zeroes first run
  if(LYZ1SUM) {
    AliFlowAnalysisWithLeeYangZeros* lyz1sum = new AliFlowAnalysisWithLeeYangZeros();
    lyz1sum->SetFirstRun(kTRUE);
    lyz1sum->SetUseSum(kTRUE);
    lyz1sum->Init();
  }
  if(LYZ1PROD) {
    AliFlowAnalysisWithLeeYangZeros* lyz1prod = new AliFlowAnalysisWithLeeYangZeros();
    lyz1prod->SetFirstRun(kTRUE);
    lyz1prod->SetUseSum(kFALSE);
    lyz1prod->Init();
  }
  //LYZ2 = Lee-Yang Zeroes second run
  if(LYZ2SUM) {
    AliFlowAnalysisWithLeeYangZeros* lyz2sum = new AliFlowAnalysisWithLeeYangZeros();
    // read the input file from the first run 
    TString inputFileNameLYZ2SUM = "outputLYZ1SUManalysis.root" ;
    TFile* inputFileLYZ2SUM = new TFile(inputFileNameLYZ2SUM.Data(),"READ");
    if(!inputFileLYZ2SUM || inputFileLYZ2SUM->IsZombie()) { 
      cerr << " ERROR: To run LYZ2SUM you need the output file from LYZ1SUM. This file is not there! Please run LYZ1SUM first." << endl ;
      break; 
    }
    else { 
      TList* inputListLYZ2SUM = (TList*)inputFileLYZ2SUM->Get("cobjLYZ1SUM");  
      if (!inputListLYZ2SUM) {cout<<"SUM Input list is NULL pointer!"<<endl; break;}
      else {
	cout<<"LYZ2SUM input file/list read..."<<endl;
	lyz2sum->SetFirstRunList(inputListLYZ2SUM);
	lyz2sum->SetFirstRun(kFALSE);
	lyz2sum->SetUseSum(kTRUE);
	lyz2sum->Init();
      }
    }
  }
  if(LYZ2PROD) {
    AliFlowAnalysisWithLeeYangZeros* lyz2prod = new AliFlowAnalysisWithLeeYangZeros();
    // read the input file from the first run 
    TString inputFileNameLYZ2PROD = "outputLYZ1PRODanalysis.root" ;
    TFile* inputFileLYZ2PROD = new TFile(inputFileNameLYZ2PROD.Data(),"READ");
    if(!inputFileLYZ2PROD || inputFileLYZ2PROD->IsZombie()) { 
      cerr << " ERROR: To run LYZ2PROD you need the output file from LYZ1PROD. This file is not there! Please run LYZ1PROD first." << endl ;
      break; 
    }
    else { 
      TList* inputListLYZ2PROD = (TList*)inputFileLYZ2PROD->Get("cobjLYZ1PROD");  
      if (!inputListLYZ2PROD) {cout<<"PROD Input list is NULL pointer!"<<endl; break;}
      else {
	cout<<"LYZ2PROD input file/list read..."<<endl;
	lyz2prod->SetFirstRunList(inputListLYZ2PROD);
	lyz2prod->SetFirstRun(kFALSE);
	lyz2prod->SetUseSum(kTRUE);
	lyz2prod->Init();
      }
    }
  }
 //LYZEP = Lee-Yang Zeroes event plane
  if(LYZEP) {
    AliFlowLYZEventPlane* ep = new AliFlowLYZEventPlane() ;
    AliFlowAnalysisWithLYZEventPlane* lyzep = new AliFlowAnalysisWithLYZEventPlane();
    // read the input file from the second lyz run 
    TString inputFileNameLYZEP = "outputLYZ2SUManalysis.root" ;
    TFile* inputFileLYZEP = new TFile(inputFileNameLYZEP.Data(),"READ");
    if(!inputFileLYZEP || inputFileLYZEP->IsZombie()) { 
      cerr << " ERROR: To run LYZEP you need the output file from LYZ2SUM. This file is not there! Please run LYZ2SUM first." << endl ; 
      break;
    }
    else { 
      TList* inputListLYZEP = (TList*)inputFileLYZEP->Get("cobjLYZ2SUM");  
      if (!inputListLYZEP) {cout<<"Input list is NULL pointer!"<<endl; break;}
      else {
	cout<<"LYZEP input file/list read..."<<endl;
	ep   ->SetSecondRunList(inputListLYZEP);
	lyzep->SetSecondRunList(inputListLYZEP);
	ep   ->Init();
	lyzep->Init();
      }
    }
  }
  // MH = Mixed Harmonics:  
  if(MH) { 
    AliFlowAnalysisWithMixedHarmonics* mh = new AliFlowAnalysisWithMixedHarmonics();
    if(listWithWeights) mh->SetWeightsList(listWithWeights);
    //if(usePhiWeights) mh->SetUsePhiWeights(usePhiWeights); // to be improved (enabled)
    //if(usePtWeights) mh->SetUsePtWeights(usePtWeights); // to be improved (enabled)
    //if(useEtaWeights) mh->SetUseEtaWeights(useEtaWeights); // to be improved (enabled)
    mh->Init();
  }
  // NL = Nested Loops:  
  if(NL) { 
    AliFlowAnalysisWithNestedLoops* nl = new AliFlowAnalysisWithNestedLoops();
    if(listWithWeights) nl->SetWeightsList(listWithWeights);
    //if(usePhiWeights) nl->SetUsePhiWeights(usePhiWeights); // to be improved (enabled)
    //if(usePtWeights) nl->SetUsePtWeights(usePtWeights); // to be improved (enabled)
    //if(useEtaWeights) nl->SetUseEtaWeights(useEtaWeights); // to be improved (enabled)
    nl->Init();
  }

  //------------------------------------------------------------------------
  
  
  //standard code to read files in directory
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

	  Double_t xRPAngle;	  
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

	      /*
	      Int_t nNewMultOfEvent = random3Temp.Gaus(nMultiplicityOfEvent,nSigmaMult);
	      cout << "new multiplicity: " << nNewMultOfEvent << endl;
  	      Double_t xNewFlowValue = random3Temp.Gaus(xEllipticFlowValue,xSigmaFlow);
	      if ( (fCount % 100) == 0) {
		cout << "new multiplicity: " << nNewMultOfEvent << endl;
		cout << "new flow value: " << xNewFlowValue << endl;
	      }
	      */

	      AliFlowEventSimple *fEvent = new AliFlowEventSimple(kTree, cutsInt, cutsDiff); 
	      //xRPAngle=random3Temp.Uniform(0.0,TMath::TwoPi());
        //fEvent->SetMCReactionPlaneAngle(xRPAngle);
        //fEvent->SetV2(xNewFlowValue);
        //fEvent->CloneTracks(nLoops);
	                    
	      // do flow analysis for various methods
	      if(MCEP)    mcep->Make(fEvent);
	      if(QC)      qc->Make(fEvent);
	      if(GFC)     gfc->Make(fEvent);
	      if(FQD)     fqd->Make(fEvent);
	      if(LYZ1SUM) lyz1sum->Make(fEvent);
	      if(LYZ1PROD)lyz1prod->Make(fEvent);
	      if(LYZ2SUM) lyz2sum->Make(fEvent);
	      if(LYZ2PROD)lyz2prod->Make(fEvent);
	      if(LYZEP)   lyzep->Make(fEvent,ep);
	      if(SP)      sp->Make(fEvent);	      
	      if(MH)      mh->Make(fEvent);
	      if(NL)      nl->Make(fEvent);
	      //-----------------------------------------------------------
	      fCount++;
	      //cout << "# " << fCount << " events processed" << endl;
	      delete kTree;
	      delete fEvent;
	    }
	  delete kineFile ;
	}
      delete evtsDir ;
    }

 //---------------------------------------------------------------------------------------  
 // create a new file which will hold the final results of all methods:
 TString outputFileName = "AnalysisResults.root";  
 TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
 // create a new file for each method wich will hold list with final results:
 const Int_t nMethods = 12;
 TString method[nMethods] = {"MCEP","SP","GFC","QC","FQD","LYZ1SUM","LYZ1PROD","LYZ2SUM","LYZ2PROD","LYZEP","MH","NL"};
 TDirectoryFile *dirFileFinal[nMethods] = {NULL};
 TString fileNameMethod[nMethods]; 
 for(Int_t i=0;i<nMethods;i++)
 {
  // form a file name for each method:
  fileNameMethod[i]+="output";
  fileNameMethod[i]+=method[i].Data();
  fileNameMethod[i]+="analysis";
  dirFileFinal[i] = new TDirectoryFile(fileNameMethod[i].Data(),fileNameMethod[i].Data());
 } 
 
 // calculating and storing the final results of default flow analysis:
 if(MCEP)    {mcep->Finish();    mcep->WriteHistograms(dirFileFinal[0]);}
 if(SP)      {sp->Finish();      sp->WriteHistograms(dirFileFinal[1]);}
 if(GFC)     {gfc->Finish();     gfc->WriteHistograms(dirFileFinal[2]);}
 if(QC)      {qc->Finish();      qc->WriteHistograms(dirFileFinal[3]);}
 if(FQD)     {fqd->Finish();     fqd->WriteHistograms(dirFileFinal[4]);}
 if(LYZ1SUM) {lyz1sum->Finish(); lyz1sum->WriteHistograms(dirFileFinal[5]);}
 if(LYZ1PROD){lyz1prod->Finish();lyz1prod->WriteHistograms(dirFileFinal[6]);}
 if(LYZ2SUM) {lyz2sum->Finish(); lyz2sum->WriteHistograms(dirFileFinal[7]);}
 if(LYZ2PROD){lyz2prod->Finish();lyz2prod->WriteHistograms(dirFileFinal[8]);}
 if(LYZEP)   {lyzep->Finish();   lyzep->WriteHistograms(dirFileFinal[9]);}
 if(MH)      {mh->Finish();      mh->WriteHistograms(dirFileFinal[10]);}
 if(NL)      {nl->Finish();      nl->WriteHistograms(dirFileFinal[11]);}
 //---------------------------------------------------------------------------------------  
 
 outputFile->Close();
 delete outputFile;
  
  cout << endl;
  cout << " Fini ... " << endl;
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

void LoadLibraries(const anaModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  //gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    cerr<<"libCORRFW loaded..."<<endl;
    gSystem->Load("libPWG2flowCommon");
    cerr<<"libPWG2flowCommon loaded..."<<endl;
    gSystem->Load("libPWG2flowTasks");
    cerr<<"libPWG2flowTasks loaded..."<<endl;
  }
  
  else if (mode == mLocalPAR) {
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
     //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    SetupPar("PWG2AOD");
    SetupPar("CORRFW");
    SetupPar("PWG2flowCommon");
    cerr<<"PWG2flowCommon.par loaded..."<<endl;
    SetupPar("PWG2flowTasks");
    cerr<<"PWG2flowTasks.par loaded..."<<endl;
  }
  
  //---------------------------------------------------------
  // <<<<<<<<<< Source mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==mLocalSource) {
 
    // In root inline compile

   
    // Constants  
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonConstants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZConstants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCumuConstants.cxx+");
    
    // Flow event
    gROOT->LoadMacro("AliFlowCommon/AliFlowVector.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimple.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowEvent.cxx+");
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
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZEventPlane.cxx+");
    
    // Flow Analysis code for various methods
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithMCEventPlane.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithScalarProduct.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithLYZEventPlane.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithLeeYangZeros.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithCumulants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithQCumulants.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithFittingQDistribution.cxx+");    
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithMixedHarmonics.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithNestedLoops.cxx+");
    
    cout << "finished loading macros!" << endl;  
    
  }  
  
}


