//DEFINITION OF A FEW CONSTANTS
const Double_t ymin  = -2.0 ;
const Double_t ymax  =  2.0 ;
const Double_t ptmin =  0.0 ;
const Double_t ptmax =  8.0 ;
const Double_t phimin =   0 ;
const Double_t phimax = 360 ;
const Int_t    mintrackrefsTPC = 2 ;
const Int_t    mintrackrefsITS = 3 ;
const Int_t    PDG = 310; 
const Int_t    minclustersTPC = 10 ;
const Int_t    chargeV0 = 0 ;
//----------------------------------------------------

Bool_t AliCFV0Task(
		   const Bool_t useGrid = 1,
		   const Bool_t readAOD = 0,
		   const char * kTagXMLFile="wn.xml", // XML file containing tags
		   )
{
  
  TBenchmark benchmark;
  benchmark.Start("AliSingleTrackTask");

  AliLog::SetGlobalDebugLevel(0);

  Load(useGrid) ; //load the required libraries

  AliLog::SetGlobalDebugLevel(0);

  TChain * analysisChain ;

  if (useGrid) { //data located on AliEn
    TGrid::Connect("alien://") ;
    //  Create an AliRunTagCuts and an AliEventTagCuts Object and impose some selection criteria
    AliRunTagCuts      *runCuts   = new AliRunTagCuts(); 
    AliEventTagCuts    *eventCuts = new AliEventTagCuts(); 
    AliLHCTagCuts      *lhcCuts   = new AliLHCTagCuts(); 
    AliDetectorTagCuts *detCuts   = new AliDetectorTagCuts(); 
    eventCuts->SetMultiplicityRange(0,20000);
    //  Create an AliTagAnalysis Object and chain the tags
    AliTagAnalysis   *tagAna = new AliTagAnalysis(); 
    tagAna->SetType("ESD");  //for aliroot > v4-05
    TAlienCollection *coll   = TAlienCollection::Open(kTagXMLFile); 
    TGridResult      *tagResult = coll->GetGridResult("",0,0);
    tagResult->Print();
    tagAna->ChainGridTags(tagResult);
    //  Create a new esd chain and assign the chain that is returned by querying the tags
    analysisChain = tagAna->QueryTags(runCuts,lhcCuts,detCuts,eventCuts); 
  }
  else {// local data
    //here put your input data path
    if (readAOD) {
      analysisChain = new TChain("aodTree");
      analysisChain->Add("AliAOD.root");
    }
    else {
      analysisChain = new TChain("esdTree");
      analysisChain->Add("AliESDs.root");
    }
  }
  
  
  Info("AliCFV0Task",Form("CHAIN HAS %d ENTRIES",(Int_t)analysisChain->GetEntries()));

  //CONTAINER DEFINITION
  Info("AliCFV0Task","SETUP CONTAINER");
  //the sensitive variables, their indices
  Int_t ipt = 0;
  Int_t iy  = 1;
  //Setting up the container grid... 
  Int_t nstep = 4 ; //number of selection steps MC 
  const Int_t nvar   = 2 ; //number of variables on the grid:pt,y,phi,vtx
  const Int_t nbin1  = 8 ; //bins in pt
  const Int_t nbin2  = 8 ; //bins in y 
  //arrays for the number of bins in each dimension
  const Int_t iBin[nvar] ={nbin1,nbin2};
  //arrays for lower bounds :
  Double_t binLim1[nbin1+1];
  Double_t binLim2[nbin2+1];
  //values for bin lower bounds
  for(Int_t i=0; i<=nbin1; i++) binLim1[i]=ptmin + (ptmax-ptmin)/nbin1*(Double_t)i ; 
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=ymin  + (ymax-ymin)  /nbin2*(Double_t)i ;
  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("container","container for V0s",nstep,nvar,iBin);
  //setting the bin limits
  container -> SetBinLimits(ipt,binLim1);
  container -> SetBinLimits(iy,binLim2);


  //CREATE THE  CUTS -----------------------------------------------
  
  // Gen-Level kinematic cuts
  AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
  mcKineCuts->SetPtRange(ptmin,ptmax);
  mcKineCuts->SetRapidityRange(ymin,ymax);
  mcKineCuts->SetChargeMC(0);

  AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
  //mcGenCuts->SetRequireIsPrimary(); //problem with some particles...
  mcGenCuts->SetRequirePdgCode(PDG);

  //Acceptance Cuts
  AliCFPairAcceptanceCuts *mcAccCuts = new AliCFPairAcceptanceCuts("mcAccCuts","MC acceptance cuts for V0");
  mcAccCuts->SetMinNHitITS(mintrackrefsITS,mintrackrefsITS);
  mcAccCuts->SetMinNHitTPC(mintrackrefsTPC,mintrackrefsTPC);

  // Rec-Level kinematic cuts
  AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts","V0 rec-level kine cuts");
  recKineCuts->SetPtRange(ptmin,ptmax);
  recKineCuts->SetRapidityRange(ymin,ymax);
  recKineCuts->SetChargeRec(0);

  AliCFPairQualityCuts *recQualityCuts = new AliCFPairQualityCuts("recQualityCuts","V0 rec-level quality cuts");
  if (!readAOD) recQualityCuts->SetMinNClusterTPC(minclustersTPC,minclustersTPC);
  recQualityCuts->SetStatus(AliESDtrack::kTPCrefit & AliESDtrack::kITSrefit, 
			    AliESDtrack::kTPCrefit & AliESDtrack::kITSrefit) ;

  AliCFV0TopoCuts *recTopoCuts = new AliCFV0TopoCuts("recTopoCuts","V0 Topological Cuts");
  recTopoCuts->SetMaxDcaDaughters(0.1);
  recTopoCuts->SetMinCosPointAngle(0.995);
  recTopoCuts->SetMinDcaNeg(0.1);
  recTopoCuts->SetMinDcaPos(0.1);


  AliCFPairPidCut* cutPID = new AliCFPairPidCut("cutPID","ESD_PID") ;
  Double_t prior_pp[AliPID::kSPECIES] = {0.0244519,
					 0.0143988,
					 0.805747 ,
					 0.0928785,
					 0.0625243 };

  Double_t prior_pbpb[AliPID::kSPECIES] = {0.0609,
					   0.1064,
					   0.7152 ,
					   0.0442,
					   0.0733 };
  cutPID->SetPriors(prior_pp);
  cutPID->SetDetectors("ITS TPC TRD","ITS TPC TRD");
  if (readAOD) cutPID->SetAODmode(kTRUE);
  else         cutPID->SetAODmode(kFALSE);
  cutPID->SetProbabilityCut(0,0);
  switch(TMath::Abs(PDG)) {
  case  310    : cutPID->SetParticleType(AliPID::kPion  , kTRUE, AliPID::kPion  , kTRUE); break;
  case  3122   : cutPID->SetParticleType(AliPID::kPion  , kTRUE, AliPID::kProton, kTRUE); break;
  case -3122   : cutPID->SetParticleType(AliPID::kProton, kTRUE, AliPID::kPion  , kTRUE); break;
  default      : printf("UNDEFINED PID\n"); break;
  }

  Info("AliCFV0Task","CREATE MC KINE CUTS");
  TObjArray* mcList = new TObjArray(0) ;
  mcList->AddLast(mcKineCuts);
  mcList->AddLast(mcGenCuts);

  Info("AliCFV0Task","CREATE ACCEPTANCE CUTS");
  TObjArray* accList = new TObjArray(0) ;
  accList->AddLast(mcAccCuts);

  Info("AliCFV0Task","CREATE RECONSTRUCTION CUTS");
  TObjArray* recList = new TObjArray(0) ;
  recList->AddLast(recKineCuts);
  recList->AddLast(recQualityCuts);
  recList->AddLast(recTopoCuts);

  Info("AliCFV0Task","CREATE PID CUTS");
  TObjArray* fPIDCutList = new TObjArray(0) ;
  fPIDCutList->AddLast(cutPID);


  //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  Info("AliCFV0Task","CREATE INTERFACE AND CUTS");
  AliCFManager* man = new AliCFManager() ;
  man->SetParticleContainer(container);
  man->SetParticleCutsList (AliCFManager::kPartGenCuts,mcList);
  man->SetParticleCutsList (AliCFManager::kPartAccCuts,accList);
  man->SetParticleCutsList (AliCFManager::kPartRecCuts,recList);
  man->SetParticleCutsList (AliCFManager::kPartSelCuts,fPIDCutList);

  //CREATE THE TASK
  Info("AliCFV0Task","CREATE TASK");
  // create the task
  AliCFV0Task *task = new AliCFV0Task("AliCFV0Task");
  task->SetCFManager(man); //here is set the CF manager
  if (!readAOD)      task->SetRebuildV0s(kTRUE);
  task->SetV0PDG(PDG);

  //SETUP THE ANALYSIS MANAGER TO READ INPUT CHAIN AND WRITE DESIRED OUTPUTS
  Info("AliCFV0Task","CREATE ANALYSIS MANAGER");
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  if (useGrid) mgr->SetAnalysisType(AliAnalysisManager::kGridAnalysis);
  else mgr->SetAnalysisType(AliAnalysisManager::kLocalAnalysis);

  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);

  AliInputEventHandler* dataHandler ;

  if   (readAOD) dataHandler = new AliAODInputHandler();
  else           dataHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(dataHandler);

  // Create and connect containers for input/output

  //input data
  AliAnalysisDataContainer *cinput0  = 
    mgr->CreateContainer("cchain0",TChain::Class(),AliAnalysisManager::kInputContainer);
  //slot 0 : default output tree (by default handled by AliAnalysisTaskSE)
  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("ctree0", TTree::Class(),AliAnalysisManager::kOutputContainer,"outputV0.root");
  // output histo (number of events processed)
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("chist0", TH1I::Class(),AliAnalysisManager::kOutputContainer,"outputV0.root");
  // output Correction Framework Container (for acceptance & efficiency calculations)
  AliAnalysisDataContainer *coutput2 = 
    mgr->CreateContainer("ccontainer0", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,"outputV0.root");

  cinput0->SetData(analysisChain);

  mgr->AddTask(task);
  mgr->ConnectInput (task,0,cinput0 );
  mgr->ConnectOutput(task,0,coutput0);
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);

 
  Info("AliCFV0Task","READY TO RUN");
  //RUN !!!
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",analysisChain);
  }

  benchmark.Stop("AliCFV0Task");
  benchmark.Show("AliCFV0Task");
  
  return kTRUE ;
}

void Load(Bool_t useGrid) {

  //load the required aliroot libraries
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;
  gSystem->Load("libCORRFW.so") ;

  //compile online the task class
  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ROOTSYS/include");
  gROOT->LoadMacro("./AliCFV0Task.cxx+");
}
