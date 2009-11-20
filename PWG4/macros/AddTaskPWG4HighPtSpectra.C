//DEFINITION OF A FEW CONSTANTS
const Double_t ptmin =  2.0 ;
const Double_t ptmax =  100.0 ;

const Int_t    mintrackrefsTPC = 1;
const Int_t    mintrackrefsITS = 1;
const Int_t    charge  = 1 ;

AliPWG4HighPtSpectra* AddTaskPWG4HighPtSpectra()//<some_parameters>)
{
  // Creates a proton analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPWG4HighPtSpectra", "No analysis manager to connect to.");
    return NULL;
  }  

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHighPtSpectra", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  const char *analysisType = "ESD";//"TPC"

  // C. Create the task, add it to manager.
  //===========================================================================
  //CONTAINER DEFINITION
  Info("AliPWG4HighPtSpectra","SETUP CONTAINER");
  //the sensitive variables, their indices
  UInt_t ipt = 0;
  //Setting up the container grid... 
  UInt_t nstep = 5 ; //Steps/Modes for containers
  Int_t kStepReconstructed = 0;
  Int_t kStepReconstructedTPCOnly = 1;
  Int_t kStepSecondaries = 2;
  Int_t kStepMCtrackable = 3;
  Int_t kStepReconstructedMC = 4;

  const Int_t nvar   = 1; //number of variables on the grid:pt
  const Int_t nbin1  = 98; //bins in pt 98 
 
  //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbin1;

  //arrays for lower bounds :
  Double_t *binLim1=new Double_t[nbin1+1];

  //values for bin lower bounds
  for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)ptmin + (ptmax-ptmin)/nbin1*(Double_t)i ; 
 
  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("container","container for tracks",nstep,nvar,iBin);
  //setting the bin limits
  container -> SetBinLimits(ipt,binLim1);
 
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
  //Standard Cuts
  trackCuts->SetAcceptKinkDaughters(kFALSE);//
  trackCuts->SetRequireTPCRefit(kTRUE);//
  trackCuts->SetEtaRange(-0.9,0.9);//-0.5,0.5);//
  trackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);//
  trackCuts->SetPtRange(0.15, 1e10);//
  trackCuts->SetMinNClustersTPC(50);//
  trackCuts->SetMaxChi2PerClusterTPC(3.5);//
  trackCuts->SetRequireITSRefit(kTRUE);
  trackCuts->SetMaxDCAToVertexXY(2.4);
  trackCuts->SetMaxDCAToVertexZ(3.2);
  trackCuts->SetDCAToVertex2D(kTRUE);
 
  // Gen-Level kinematic cuts
  AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
  mcKineCuts->SetPtRange(0.15,1e10); 
  mcKineCuts->SetRapidityRange(-0.9,0.9);//-0.5,0.5);
  mcKineCuts->SetRequireIsCharged(kTRUE);

  //Acceptance Cuts
  AliCFAcceptanceCuts *mcAccCuts = new AliCFAcceptanceCuts("mcAccCuts","MC acceptance cuts");
  // mcAccCuts->SetMinNHitITS(mintrackrefsITS);
  mcAccCuts->SetMinNHitTPC(mintrackrefsTPC);

  TObjArray* recList = new TObjArray(0);
  TObjArray* recMCList = new TObjArray(0);
  TObjArray* recTPConlyList = new TObjArray(0);
  TObjArray* secList = new TObjArray(0) ;

  printf("CREATE MC KINE CUTS\n");
  TObjArray* mcList = new TObjArray(0) ;
  mcList->AddLast(mcKineCuts);
  mcList->AddLast(mcAccCuts);

  //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* man = new AliCFManager() ;
  man->SetParticleContainer(container);
  man->SetParticleCutsList(0,recList);
  man->SetParticleCutsList(1,recTPConlyList);
  man->SetParticleCutsList(2,secList);
  man->SetParticleCutsList(3,mcList);
  man->SetParticleCutsList(4,recMCList);


  AliPWG4HighPtSpectra *taskPWG4HighPtSpectra = new AliPWG4HighPtSpectra("taskPWG4HighPtSpectra");
  taskPWG4HighPtSpectra->SetCuts(trackCuts);
  taskPWG4HighPtSpectra->SetCFManager(man); //here is set the CF manager
 
  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  //------ input data ------
  //  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
//  char *outputfile = "outputAliPWG4HighPtSpectraTestTrain.root";
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG4_HighPtSpectra"; 
  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("chist0", TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("ccontainer0", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);
 
  mgr->AddTask(taskPWG4HighPtSpectra);

  mgr->ConnectInput(taskPWG4HighPtSpectra,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4HighPtSpectra,0,coutput0);
  mgr->ConnectOutput(taskPWG4HighPtSpectra,1,coutput1);

  // Return task pointer at the end
  return taskPWG4HighPtSpectra;
}
