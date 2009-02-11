// DEFINITION OF A FEW CONSTANTS

const Double_t nevtmin= 1;
const Double_t nevtmax = 15000;

// Muons 
const Double_t ymin  = -4.0 ;
const Double_t ymax  =  -2.5 ;

const Double_t phimin = -180;
const Double_t phimax = 180;

// Resonance
const Int_t    PDG = 443;

const Double_t ptmin =  0.0 ;
const Double_t ptmax =  30 ;
const Double_t pmin =  0.0 ;
const Double_t pmax =  700 ;
const Int_t    charge  = 0 ;
const Double_t mmin =  0.1 ;
const Double_t mmax =  6 ;
const Double_t mymin =  -5 ;
const Double_t mymax =  -1.5 ;

//----------------------------------------------------

Bool_t AliCFMuonResTask1(
			    const Bool_t useGrid = 0,
			    const Bool_t readAOD = 0,
			    const char * kTagXMLFile="wn.xml" // XML file containing tags
			    )
{
  
  TBenchmark benchmark;
  benchmark.Start("AliMuonResTask1");

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
      analysisChain->Add("/scratch/lopez/PDC08jpsi/run2-300/AliESDs.root");
   }
  }
  
///// END INPUT


  Info("AliCFMuonResTask1",Form("CHAIN HAS %d ENTRIES",(Int_t)analysisChain->GetEntries()));

  // CONTAINER DEFINITION
  Info("AliCFMuonResTask1","SETUP CONTAINER");
  
  // The sensitive variables (9 in this example), their indices
  UInt_t nevt  = 0;
  UInt_t y1  = 1;
  UInt_t phi1  = 2;
  UInt_t y2  = 3;
  UInt_t phi2  = 4;
  UInt_t imass  = 5;
  UInt_t y  = 6;
  UInt_t pt = 7;
  UInt_t p = 8;


  // Setting up the container grid
  UInt_t nstep = 2 ; //number of selection steps : MC and ESD 
  const Int_t nvar   = 9 ;     //number of variables on the grid
  const Int_t nbin1  = nevtmax ;  
  const Int_t nbin2  = 100 ;  
  const Int_t nbin3  = 360 ;  
  const Int_t nbin4  = 100 ;  
  const Int_t nbin5  = 360 ;  
  const Int_t nbin6  = 100 ;  
  const Int_t nbin7  = 100 ;  
  const Int_t nbin8  = 100 ;  
  const Int_t nbin9  = 100 ;  

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

  // values for bin lower bounds
  for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)nevtmin  + (nevtmax-nevtmin)  /nbin1*(Double_t)i ;
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)ymin  + (ymax-ymin)  /nbin2*(Double_t)i ;
  for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)phimin  + (phimax-phimin)  /nbin3*(Double_t)i ;
  for(Int_t i=0; i<=nbin4; i++) binLim4[i]=(Double_t)ymin  + (ymax-ymin)  /nbin4*(Double_t)i ;
  for(Int_t i=0; i<=nbin5; i++) binLim5[i]=(Double_t)phimin  + (phimax-phimin)  /nbin5*(Double_t)i ;
  for(Int_t i=0; i<=nbin6; i++) binLim6[i]=(Double_t)mmin  + (mmax-mmin)  /nbin6*(Double_t)i ;
  for(Int_t i=0; i<=nbin7; i++) binLim7[i]=(Double_t)mymin  + (mymax-mymin)  /nbin7*(Double_t)i ;
  for(Int_t i=0; i<=nbin8; i++) binLim8[i]=(Double_t)ptmin + (ptmax-ptmin)/nbin8*(Double_t)i ; 
  for(Int_t i=0; i<=nbin9; i++) binLim9[i]=(Double_t)pmin + (pmax-pmin)/nbin9*(Double_t)i ; 

  // one container  of 2 steps (MC and ESD) with 9 variables
  AliCFContainer* container = new AliCFContainer("container","container for tracks",nstep,nvar,iBin);
  // setting the bin limits
  container -> SetBinLimits(nevt,binLim1);
  container -> SetBinLimits(y1,binLim2);
  container -> SetBinLimits(phi1,binLim3);
  container -> SetBinLimits(y2,binLim4);
  container -> SetBinLimits(phi2,binLim5);
  container -> SetBinLimits(imass,binLim6);
  container -> SetBinLimits(y,binLim7);
  container -> SetBinLimits(pt,binLim8);
  container -> SetBinLimits(p,binLim9);

  // Set list
  TList* qaList = new TList();

  //CREATE THE CUTS
  // Choice of the Resonance
  AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
  mcGenCuts->SetRequirePdgCode(PDG);
  mcGenCuts->SetQAOn(qaList);

  // Set a pt range of the resonance
  AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
  mcKineCuts->SetChargeMC(charge);
  mcKineCuts->SetPtRange(ptmin,ptmax);
  mcKineCuts->SetQAOn(qaList);

  // Create and fill the list associated 
  TObjArray* mcList = new TObjArray(0) ;
  mcList->AddLast(mcKineCuts);
  mcList->AddLast(mcGenCuts);

  // kinematic cuts on muons rapidity 
  AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts","rec-level kine cuts");
  recKineCuts->SetRapidityRange(ymin,ymax);
  recKineCuts->SetQAOn(qaList);
  TObjArray* recList = new TObjArray(0) ;
  recList->AddLast(recKineCuts);

  // CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* man = new AliCFManager() ;
  man->SetParticleContainer     (container);
  man->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList);
  man->SetParticleCutsList(AliCFManager::kPartRecCuts,recList);


  //CREATE THE TASK
  printf("CREATE TASK\n");
  // create the task
  AliCFMuonResTask1 *task = new AliCFMuonResTask1("AliMuonResTask1");
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
  sprintf(file,"CFMuonResTask1.root");
  printf("Analysis output in %s \n",file);

  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TH1I::Class(),AliAnalysisManager::kOutputContainer,file);
  // output Correction Framework Container (for acceptance & efficiency calculations)
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ccontainer0", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,file);

  cinput0->SetData(analysisChain);
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput0);
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);

  printf("READY TO RUN\n");
  //RUN !!!
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",analysisChain);
  }

  benchmark.Stop("AliMuonResTask1");
  benchmark.Show("AliMuonResTask1");

  return kTRUE ;
}

void Load() {

  //load the required aliroot libraries
  gSystem->Load("libANALYSIS") ;

  gSystem->Load("libANALYSISalice") ;

  //gSystem->Load("libCORRFW.so") ;
  gSystem->Load("$ALICE_ROOT/lib/tgt_linux/libCORRFW.so") ;

  //compile online the task class
  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ROOTSYS/include");
  gROOT->LoadMacro("./AliCFMuonResTask1.cxx+");
}
