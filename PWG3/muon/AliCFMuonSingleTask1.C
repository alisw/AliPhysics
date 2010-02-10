// DEFINITION OF A FEW CONSTANTS

// File scan

const Int_t nsub = 800;

// Event var

const Double_t nevtmin= 1;
const Double_t nevtmax = 500000;

const Double_t trgmin=0;
const Double_t trgmax=10;

const Double_t vmin=-100;
const Double_t vmax=100;

const Double_t vzmin=-3000;
const Double_t vzmax=1000;

// Muons

const Int_t    PDG = 13; 

const Double_t ymin  = -4;
const Double_t ymax  =  -2.5;

const Double_t phimin = -180;
const Double_t phimax = 180;

const Double_t ptmin =  0;
const Double_t ptmax =  20;

const Double_t pmin =  0;
const Double_t pmax =  800;

const Double_t hitmin=0;
const Double_t hitmax=20;

const Double_t chi2min=0;
const Double_t chi2max=20;

const Double_t dcamin=0;
const Double_t dcamax=50;

//----------------------------------------------------

Bool_t AliCFMuonSingleTask1( Int_t runmin = 170119, Int_t runmax = 170119,
			    const Bool_t useGrid = 0,
			    const Bool_t readAOD = 0,
			    const char * kTagXMLFile="wn.xml" // XML file containing tags
			    )
{
  
  TBenchmark benchmark;
  benchmark.Start("AliMuonSingleTask1");
  AliLog::SetGlobalDebugLevel(0);
  Load() ; // load the required libraries
  TChain * analysisChain ;

///// INPUT

  if (useGrid) { // data located on AliEn
    TGrid::Connect("alien://") ;    //  Create an AliRunTagCuts and an AliEventTagCuts Object 
                                    //  and impose some selection criteria
    AliRunTagCuts      *runCuts   = new AliRunTagCuts(); 
    AliEventTagCuts    *eventCuts = new AliEventTagCuts(); 
    AliLHCTagCuts      *lhcCuts   = new AliLHCTagCuts(); 
    AliDetectorTagCuts *detCuts   = new AliDetectorTagCuts(); 
    eventCuts->SetMultiplicityRange(0,2000);

    //  Create an AliTagAnalysis Object and chain the tags
    AliTagAnalysis   *tagAna = new AliTagAnalysis(); 
    if (readAOD) tagAna->SetType("AOD");  // for aliroot > v4-05
    else         tagAna->SetType("ESD");  // for aliroot > v4-05
    TAlienCollection *coll   = TAlienCollection::Open(kTagXMLFile); 
    TGridResult      *tagResult = coll->GetGridResult("",0,0);
    tagResult->Print();
    tagAna->ChainGridTags(tagResult);

    // Create a new esd chain and assign the chain that is returned by querying the tags
    analysisChain = tagAna->QueryTags(runCuts,lhcCuts,detCuts,eventCuts); 
  }

  else {// local data
    // here put your input data path
    printf("\n\nRunning on local file, please check the path\n\n");

    if (readAOD) {
      analysisChain = new TChain("aodTree");
      analysisChain->Add("AliAOD.root");
    }
    else {
	analysisChain = new TChain("esdTree");

// my flat single
        Char_t RunFile[256];
        for(Int_t i=1; i<=160; i++){
	    sprintf(RunFile,"/dalice02/lopez/single/myprod/4D_YCUT/run%d-500/AliESDs.root",i);
	    analysisChain->Add(RunFile);
        }

/*
// input ascii list 
 	Char_t asciFile[256];
	for(Int_t i=runmin; i<=runmax; i++){    
	    sprintf(asciFile,"/dalice05/vulpescu/scratch/PDC08/LHC08t/LIST/ESDlist%d.txt",i);
	    Char_t esdFileName[256];
	    ifstream inf(asciFile, ios::in);
	    Int_t nFile = 0;
	    while (!inf.eof()) {
		inf >> esdFileName;
		nFile++;
		if (nFile > nsub) break;
		analysisChain->Add(esdFileName);
	    }
	}
*/
    }
  }
  
///// END INPUT

  Info("AliCFMuonSingleTask1",Form("CHAIN HAS %d ENTRIES",(Int_t)analysisChain->GetEntries()));

  // CONTAINER DEFINITION
  Info("AliCFMuonSingleTask1","SETUP CONTAINER");
  
  // The sensitive variables (15 in this example), their indices
  UInt_t nevt  = 0;
  UInt_t y  = 1;
  UInt_t phi  = 2;
  UInt_t pt = 3;
  UInt_t p = 4;
  UInt_t hit = 5;
  UInt_t chi2 = 6;
  UInt_t matchtrig = 7;
  UInt_t chi2match = 8;
  UInt_t vx = 9;
  UInt_t vy = 10;
  UInt_t vz = 11;
  UInt_t trg = 12;
  UInt_t dca = 13;
  UInt_t zcoor = 14;
  

  // Setting up the container grid
  UInt_t nstep = 2 ; //number of selection steps : MC and ESD
  const Int_t nvar   = 15 ;     //number of variables on the grid
  const Int_t nbin1  = nevtmax ;  
  const Int_t nbin2  = 6 ;     // y 
  const Int_t nbin3  = 45 ;    // phi  
  const Int_t nbin4  = 50 ;    // pt   
  const Int_t nbin5  = 800 ;   // p
  const Int_t nbin6  = 20 ;    // hit
  const Int_t nbin7  = 20 ;    // chi2
  const Int_t nbin8  = 10 ;    // match trigger
  const Int_t nbin9  = 20 ;    // chi2 match
  const Int_t nbin10  = 100 ;  // vx
  const Int_t nbin11  = 100 ;  // vy
  const Int_t nbin12  = 100 ;  // vz
  const Int_t nbin13  = 10 ;   // trg
  const Int_t nbin14  = 100 ;  // dca
  const Int_t nbin15  = 1000 ; // zcoor

  // arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbin1;
  iBin[1]=nbin2;
  iBin[2]=nbin3;
  iBin[3]=nbin4;
  iBin[4]=nbin5;
  iBin[5]=nbin6;
  iBin[6]=nbin7;
  iBin[7]=nbin8;
  iBin[8]=nbin9;
  iBin[9]=nbin10;
  iBin[10]=nbin11;
  iBin[11]=nbin12;
  iBin[12]=nbin13;
  iBin[13]=nbin14;
  iBin[14]=nbin15;

  // arrays for lower bounds :
  Double_t *binLim1=new Double_t[nbin1+1];
  Double_t *binLim2=new Double_t[nbin2+1];
  Double_t *binLim3=new Double_t[nbin3+1];
  Double_t *binLim4=new Double_t[nbin4+1];
  Double_t *binLim5=new Double_t[nbin5+1];
  Double_t *binLim6=new Double_t[nbin6+1];
  Double_t *binLim7=new Double_t[nbin7+1]; 
  Double_t *binLim8=new Double_t[nbin8+1];
  Double_t *binLim9=new Double_t[nbin9+1];
  Double_t *binLim10=new Double_t[nbin10+1];
  Double_t *binLim11=new Double_t[nbin11+1];
  Double_t *binLim12=new Double_t[nbin12+1];
  Double_t *binLim13=new Double_t[nbin13+1];
  Double_t *binLim14=new Double_t[nbin14+1];
  Double_t *binLim15=new Double_t[nbin15+1];

  // values for bin lower bounds
  for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)nevtmin  + (nevtmax-nevtmin)  /nbin1*(Double_t)i ;
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)ymin  + (ymax-ymin)  /nbin2*(Double_t)i ;
  for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)phimin  + (phimax-phimin)  /nbin3*(Double_t)i ;
  for(Int_t i=0; i<=nbin4; i++) binLim4[i]=(Double_t)ptmin  + (ptmax-ptmin)  /nbin4*(Double_t)i ;
  for(Int_t i=0; i<=nbin5; i++) binLim5[i]=(Double_t)pmin  + (pmax-pmin)  /nbin5*(Double_t)i ;
  for(Int_t i=0; i<=nbin6; i++) binLim6[i]=(Double_t)hitmin  + (hitmax-hitmin)  /nbin6*(Double_t)i ;
  for(Int_t i=0; i<=nbin7; i++) binLim7[i]=(Double_t)chi2min  + (chi2max-chi2min)  /nbin7*(Double_t)i ;
  for(Int_t i=0; i<=nbin8; i++) binLim8[i]=(Double_t)trgmin  + (trgmax-trgmin)  /nbin8*(Double_t)i ;
  for(Int_t i=0; i<=nbin9; i++) binLim9[i]=(Double_t)chi2min  + (chi2max-chi2min)  /nbin9*(Double_t)i ;
  for(Int_t i=0; i<=nbin10; i++) binLim10[i]=(Double_t)vmin  + (vmax-vmin)  /nbin10*(Double_t)i ;
  for(Int_t i=0; i<=nbin11; i++) binLim11[i]=(Double_t)vmin  + (vmax-vmin)  /nbin11*(Double_t)i ;
  for(Int_t i=0; i<=nbin12; i++) binLim12[i]=(Double_t)vmin  + (vmax-vmin)  /nbin12*(Double_t)i ;
  for(Int_t i=0; i<=nbin13; i++) binLim13[i]=(Double_t)trgmin  + (trgmax-trgmin)  /nbin13*(Double_t)i ;
  for(Int_t i=0; i<=nbin14; i++) binLim14[i]=(Double_t)dcamin  + (dcamax-dcamin)  /nbin14*(Double_t)i ;
  for(Int_t i=0; i<=nbin15; i++) binLim15[i]=(Double_t)vzmin  + (vzmax-vzmin)  /nbin15*(Double_t)i ;

  // one container  of 2 steps (MC and ESD) with 15 variables
  AliCFContainer* container = new AliCFContainer("container","container for tracks",nstep,nvar,iBin);
  // setting the bin limits
  container -> SetBinLimits(nevt,binLim1);
  container -> SetBinLimits(y,binLim2);
  container -> SetBinLimits(phi,binLim3);
  container -> SetBinLimits(pt,binLim4);
  container -> SetBinLimits(p,binLim5);
  container -> SetBinLimits(hit,binLim6);
  container -> SetBinLimits(chi2,binLim7);
  container -> SetBinLimits(matchtrig,binLim8);
  container -> SetBinLimits(chi2match,binLim9);
  container -> SetBinLimits(vx,binLim10);
  container -> SetBinLimits(vy,binLim11);
  container -> SetBinLimits(vz,binLim12);
  container -> SetBinLimits(trg,binLim13);
  container -> SetBinLimits(dca,binLim14);
  container -> SetBinLimits(zcoor,binLim15);

  // Set list
  TList* qaList = new TList();

  //CREATE THE CUTS
  // Choice of the particle: here mu-
  AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
  mcGenCuts->SetRequirePdgCode(PDG);
  mcGenCuts->SetQAOn(qaList);
 
  TObjArray* mcList = new TObjArray(0) ;
  mcList->AddLast(mcGenCuts);

  // Set a pt cut of the particle
  AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
  mcKineCuts->SetPtRange(ptmin,ptmax);
  mcKineCuts->SetQAOn(qaList);
  // Set a rapidity cut 
  AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts","rec-level kine cuts");
  recKineCuts->SetRapidityRange(ymin,ymax);
  recKineCuts->SetQAOn(qaList);
  
  TObjArray* recList = new TObjArray(0) ;
  recList->AddLast(recKineCuts);
  recList->AddLast(mcKineCuts);

  // CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* man = new AliCFManager() ;
  man->SetParticleContainer     (container);
  man->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList);
  man->SetParticleCutsList(AliCFManager::kPartAccCuts,recList);

  //CREATE THE TASK
  printf("CREATE TASK\n");
  // create the task
  AliCFMuonSingleTask1 *task = new AliCFMuonSingleTask1("AliMuonSingleTask1");
  task->SetCFManager(man); //here is set the CF manager
  task->SetQAList(qaList);
  if (readAOD)       task->SetReadAODData() ;

  //SETUP THE ANALYSIS MANAGER TO READ INPUT CHAIN AND WRITE DESIRED OUTPUTS
  printf("CREATE ANALYSIS MANAGER\n");
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  if (useGrid) mgr->SetAnalysisType(AliAnalysisManager::kGridAnalysis);
  else mgr->SetAnalysisType(AliAnalysisManager::kLocalAnalysis);

  AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
 
  AliInputEventHandler* dataHandler ;
  
  if   (readAOD) dataHandler = new AliAODInputHandler();
  else           dataHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(dataHandler);

  // Create and connect containers for input/output

  // input data 
  AliAnalysisDataContainer *cinput0  = mgr->CreateContainer("cchain0",TChain::Class(),AliAnalysisManager::kInputContainer);

  // output data
  Char_t file[256];

  sprintf(file,"/dalice02/lopez/single/myprod/4D_YCUT/CFMuonSingleTask1.root");
  printf("Analysis output in %s \n",file);

  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TH1I::Class(),AliAnalysisManager::kOutputContainer,file);
  // output Correction Framework Container (for acceptance & efficiency calculations)
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ccontainer0", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,file);


  cinput0->SetData(analysisChain);
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
 
  printf("READY TO RUN\n");
  //RUN !!!
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",analysisChain);
  }
  benchmark.Stop("AliMuonSingleTask1");
  benchmark.Show("AliMuonSingleTask1");

  return kTRUE ;
}

void Load() {

  //load the required aliroot libraries
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;
  gSystem->Load("$ALICE_ROOT/lib/tgt_linux/libCORRFW.so") ;

  //compile online the task class
  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_ROOT/MUON -I$ALICE_ROOT/STEER -I$ROOTSYS/include");
  gROOT->LoadMacro("./AliCFMuonSingleTask1.cxx+");
}
