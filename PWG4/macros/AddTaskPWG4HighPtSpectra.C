//DEFINITION OF A FEW CONSTANTS
const Float_t ptmin =  2.0 ;
const Float_t ptmax =  50.0 ;
const Float_t phimin = 0.;
const Float_t phimax = 2.*TMath::Pi();
const Float_t etamin = -0.9;
const Float_t etamax = 0.9;
const Float_t dcarmin = -0.2;
const Float_t dcarmax = 0.2;
const Float_t chi2TPCmin = 0.0;
const Float_t chi2TPCmax = 3.5;

const Float_t ptmin1 =  ptmin ;
const Float_t ptmax1 =  10.0 ;
const Float_t ptmin2 =  ptmax1 ;
const Float_t ptmax2 =  20.0 ;
const Float_t ptmin3 =  ptmax2 ;
const Float_t ptmax3 =  ptmax ;

const Int_t   mintrackrefsTPC = 1;
const Int_t   mintrackrefsITS = 1;
const Int_t   charge  = 1;

AliPWG4HighPtSpectra* AddTaskPWG4HighPtSpectra()
{
  // Creates HighPtSpectra analysis task and adds it to the analysis manager.
  
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
    ::Error("AddTaskPWG4HighPtSpectra", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  const char *analysisType = "ESD";//"TPC"

  // C. Create the task, add it to manager.
  //===========================================================================
 //CONTAINER DEFINITION
  Info("AliPWG4HighPtSpectra","SETUP CONTAINER");
  //the sensitive variables, their indices
  UInt_t ipt   = 0;
  UInt_t iphi  = 1;
  UInt_t ieta  = 2;
  UInt_t idcar = 3;
  UInt_t ichi2TPC = 4;

  //Setting up the container grid... 
  UInt_t nstep = 5; //Steps/Modes for containers
  Int_t kStepReconstructed = 0;
  Int_t kStepReconstructedTPCOnly = 1;
  Int_t kStepSecondaries = 2;
  Int_t kStepMCtrackable = 3;
  Int_t kStepReconstructedMC = 4;
  
  const Int_t nvar   = 5; //number of variables on the grid:pt
  const Int_t nbin11 = (int)(ptmax1-ptmin1);
  const Int_t nbin12 = (int)((ptmax2-ptmin2)/2.)+nbin11;
  const Int_t nbin13 = (int)((ptmax3-ptmin3)/5.)+nbin12;
  const Int_t nbin1  = nbin13; //bins in pt 98 
  const Int_t nbin2  =  18;//36; //bins in phi
  const Int_t nbin3  =  9; //bins in eta
  const Int_t nbin4  =  40; //bins in DCAR
  const Int_t nbin5  =  35; //bins in Chi2/#NclusTPC

  //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbin1;
  iBin[1]=nbin2;
  iBin[2]=nbin3;
  iBin[3]=nbin4;
  iBin[4]=nbin5;
  
  //arrays for lower bounds :
  Double_t *binLim1=new Double_t[nbin1+1];
  Double_t *binLim2=new Double_t[nbin2+1];
  Double_t *binLim3=new Double_t[nbin3+1];
  Double_t *binLim4=new Double_t[nbin4+1];
  Double_t *binLim5=new Double_t[nbin5+1];

  //values for bin lower bounds 
  for(Int_t i=0; i<=nbin1; i++) {
    if(i<=nbin11) binLim1[i]=(Double_t)ptmin1 + (ptmax1-ptmin1)/nbin11*(Double_t)i ;  
    if(i<=nbin12 && i>nbin11) binLim1[i]=(Double_t)ptmin2 + (ptmax2-ptmin2)/(nbin12-nbin11)*((Double_t)i-(Double_t)nbin11) ;  
    if(i<=nbin13 && i>nbin12) binLim1[i]=(Double_t)ptmin3 + (ptmax3-ptmin3)/(nbin13-nbin12)*((Double_t)i-(Double_t)nbin12) ;  
  }
  //  for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)ptmin + (ptmax-ptmin)/nbin1*(Double_t)i ;  
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)phimin + (phimax-phimin)/nbin2*(Double_t)i ;
  for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)etamin + (etamax-etamin)/nbin3*(Double_t)i ;  
  for(Int_t i=0; i<=nbin4; i++) binLim4[i]=(Double_t)dcarmin + (dcarmax-dcarmin)/nbin4*(Double_t)i ;
  for(Int_t i=0; i<=nbin5; i++) binLim5[i]=(Double_t)chi2TPCmin + (chi2TPCmax-chi2TPCmin)/nbin5*(Double_t)i ;
  

  AliCFContainer* containerPos = new AliCFContainer("containerPos","container for positive tracks",nstep,nvar,iBin);
  //setting the bin limits
  containerPos -> SetBinLimits(ipt,binLim1);
  containerPos -> SetBinLimits(iphi,binLim2);
  containerPos -> SetBinLimits(ieta,binLim3);
  containerPos -> SetBinLimits(idcar,binLim4);
  containerPos -> SetBinLimits(ichi2TPC,binLim5);

  AliCFContainer* containerNeg = new AliCFContainer("containerNeg","container for negative tracks",nstep,nvar,iBin);
  //setting the bin limits
  containerNeg -> SetBinLimits(ipt,binLim1);
  containerNeg -> SetBinLimits(iphi,binLim2);
  containerNeg -> SetBinLimits(ieta,binLim3);
  containerNeg -> SetBinLimits(idcar,binLim4);
  containerNeg -> SetBinLimits(ichi2TPC,binLim5);
  
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
//     //Standard Cuts
//     trackCuts->SetAcceptKinkDaughters(kFALSE);
//     trackCuts->SetRequireTPCStandAlone(kTRUE); 
//     trackCuts->SetRequireTPCRefit(kTRUE);
//     trackCuts->SetMinNClustersTPC(70);
//     trackCuts->SetEtaRange(-0.9,0.9);
//     trackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
//     trackCuts->SetPtRange(0.15, 1e10);
//     trackCuts->SetMaxChi2PerClusterTPC(3.5);
//     trackCuts->SetMaxDCAToVertexXY(2.4);
//     trackCuts->SetMaxDCAToVertexZ(3.2);
//     trackCuts->SetDCAToVertex2D(kTRUE);
//     trackCuts->SetRequireITSRefit(kTRUE);
  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetPtRange(0.15, 1e10);
  trackCuts=trackCuts->GetStandardITSTPCTrackCuts2009(kTRUE);//Primary Track Selection

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
  AliCFManager* manPos = new AliCFManager("manPos","Manager for Positive tracks") ;
  manPos->SetParticleContainer(containerPos);
  manPos->SetParticleCutsList(0,recList);
  manPos->SetParticleCutsList(1,recTPConlyList);
  manPos->SetParticleCutsList(2,secList);
  manPos->SetParticleCutsList(3,mcList);
  manPos->SetParticleCutsList(4,recMCList);

  AliCFManager* manNeg = new AliCFManager("manNeg","Manager for Negative tracks") ;
  manNeg->SetParticleContainer(containerNeg);
  manNeg->SetParticleCutsList(0,recList);
  manNeg->SetParticleCutsList(1,recTPConlyList);
  manNeg->SetParticleCutsList(2,secList);
  manNeg->SetParticleCutsList(3,mcList);
  manNeg->SetParticleCutsList(4,recMCList);


  printf("Create task AliPWG4HighPtSpectra\n");
  AliPWG4HighPtSpectra *taskPWG4HighPtSpectra = new AliPWG4HighPtSpectra("taskPWG4HighPtSpectra");
  taskPWG4HighPtSpectra->SetCuts(trackCuts);
  taskPWG4HighPtSpectra->SetCFManagerPos(manPos); //here is set the CF manager
  taskPWG4HighPtSpectra->SetCFManagerNeg(manNeg); //here is set the CF manager


  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  //------ input data ------
  //  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
//  char *outputfile = "outputAliPWG4HighPtSpectraTestTrain.root";
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG4_HighPtSpectra"; 

  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("chist0HighPtSpectra", TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("ccontainer0HighPtSpectra", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ccontainer1HighPtSpectra", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  
  mgr->AddTask(taskPWG4HighPtSpectra);

  mgr->ConnectInput(taskPWG4HighPtSpectra,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4HighPtSpectra,0,coutput0);
  mgr->ConnectOutput(taskPWG4HighPtSpectra,1,coutput1);
  mgr->ConnectOutput(taskPWG4HighPtSpectra,2,coutput2);


  // Return task pointer at the end
  return taskPWG4HighPtSpectra;
}
