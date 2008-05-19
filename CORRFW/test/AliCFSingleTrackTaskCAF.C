//DEFINITION OF A FEW CONSTANTS
const Double_t ymin  = -1.0 ;
const Double_t ymax  =  1.0 ;
const Double_t ptmin =  0.0 ;
const Double_t ptmax =  8.0 ;
const Int_t    mintrackrefsTPC = 2 ;
const Int_t    mintrackrefsITS = 3 ;
const Int_t    charge  = 1 ;
const Int_t    PDG = 2212; 
const Int_t    minclustersTPC = 50 ;
//----------------------------------------------------

Bool_t AliCFSingleTrackTaskCAF()
{
  
  TBenchmark benchmark;
  benchmark.Start("AliSingleTrackTaskCAF");

  //
  // Connect to proof
  //

  TProof::Reset("proof://user@lxb6046.cern.ch");
  TProof::Open("proof://user@lxb6046.cern.ch");

  //   gProof->ClearPackage("STEERBase");
  //   gProof->ClearPackage("ESD");
  //   gProof->ClearPackage("AOD");
  //   gProof->ClearPackage("ANALYSIS");
  //   gProof->ClearPackage("ANALYSISalice");

  // Enable the STEERBase Package
  gProof->UploadPackage("STEERBase.par");
  gProof->EnablePackage("STEERBase");
  // Enable the ESD Package
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  // Enable the AOD Package
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  // Enable the Analysis Package
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");
  gProof->UploadPackage("ANALYSISalice.par");
  gProof->EnablePackage("ANALYSISalice");
  // Enable the CORRFW Package
  // gProof->ClearPackage("CORRFW");
  gProof->UploadPackage("CORRFW.par");
  gProof->EnablePackage("CORRFW");

  gProof->ShowEnabledPackages();
  gProof->Load("./AliCFSingleTrackTask.cxx+g");

  //
  // Create the chain
  //
  gROOT->LoadMacro("CreateESDChain.C");
  TChain* analysisChain = CreateESDChain("ESD1503X_v1.txt", 2);


  //CONTAINER DEFINITION
  Info("AliCFSingleTrackTaskCAF","SETUP CONTAINER");
  //the sensitive variables (2 in this example), their indices
  UInt_t ipt = 0;
  UInt_t iy  = 1;
  //Setting up the container grid... 
  UInt_t nstep = 4 ; //number of selection steps MC 
  const Int_t nvar   = 2 ; //number of variables on the grid:pt,y
  const Int_t nbin1  = 8 ; //bins in pt
  const Int_t nbin2  = 8 ; //bins in y 

  //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbin1;
  iBin[1]=nbin2;

  //arrays for lower bounds :
  Double_t *binLim1=new Double_t[nbin1+1];
  Double_t *binLim2=new Double_t[nbin2+1];

  //values for bin lower bounds
  for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)ptmin + (ptmax-ptmin)/nbin1*(Double_t)i ; 
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)ymin  + (ymax-ymin)  /nbin2*(Double_t)i ;
  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("container","container for tracks",nstep,nvar,iBin);
  //setting the bin limits
  container -> SetBinLimits(ipt,binLim1);
  container -> SetBinLimits(iy,binLim2);


  // SET TLIST FOR QA HISTOS
  TList* qaList = new TList();

  //CREATE THE  CUTS -----------------------------------------------

  // Gen-Level kinematic cuts
  AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
  mcKineCuts->SetPtRange(ptmin,ptmax);
  mcKineCuts->SetRapidityRange(ymin,ymax);
  mcKineCuts->SetChargeMC(charge);
  mcKineCuts->SetQAOn(qaList);

  //Particle-Level cuts:  
  AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
  mcGenCuts->SetRequireIsPrimary();
  mcGenCuts->SetRequirePdgCode(PDG);
  mcGenCuts->SetQAOn(qaList);

  //Acceptance Cuts
  AliCFAcceptanceCuts *mcAccCuts = new AliCFAcceptanceCuts("mcAccCuts","MC acceptance cuts");
  mcAccCuts->SetMinNHitITS(mintrackrefsITS);
  mcAccCuts->SetMinNHitTPC(mintrackrefsTPC);
  mcAccCuts->SetQAOn(qaList);

  // Rec-Level kinematic cuts
  AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts","rec-level kine cuts");
  recKineCuts->SetPtRange(ptmin,ptmax);
  recKineCuts->SetRapidityRange(ymin,ymax);
  recKineCuts->SetChargeRec(charge);
  recKineCuts->SetQAOn(qaList);

  AliCFTrackQualityCuts *recQualityCuts = new AliCFTrackQualityCuts("recQualityCuts","rec-level quality cuts");
  recQualityCuts->SetMinNClusterTPC(minclustersTPC);
  recQualityCuts->SetRequireITSRefit(kTRUE);
  recQualityCuts->SetQAOn(qaList);

  AliCFTrackIsPrimaryCuts *recIsPrimaryCuts = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts","rec-level isPrimary cuts");
  recIsPrimaryCuts->SetMaxNSigmaToVertex(3);
  recIsPrimaryCuts->SetQAOn(qaList);

  AliCFTrackCutPid* cutPID = new AliCFTrackCutPid("cutPID","ESD_PID") ;
  int n_species = AliPID::kSPECIES ;
  Double_t* prior = new Double_t[n_species];
  
  prior[0] = 0.0244519 ;
  prior[1] = 0.0143988 ;
  prior[2] = 0.805747  ;
  prior[3] = 0.0928785 ;
  prior[4] = 0.0625243 ;
  
  cutPID->SetPriors(prior);
  cutPID->SetProbabilityCut(0.0);
  cutPID->SetDetectors("TPC TOF");
  switch(TMath::Abs(PDG)) {
  case 11   : cutPID->SetParticleType(AliPID::kElectron, kTRUE); break;
  case 13   : cutPID->SetParticleType(AliPID::kMuon    , kTRUE); break;
  case 211  : cutPID->SetParticleType(AliPID::kPion    , kTRUE); break;
  case 321  : cutPID->SetParticleType(AliPID::kKaon    , kTRUE); break;
  case 2212 : cutPID->SetParticleType(AliPID::kProton  , kTRUE); break;
  default   : printf("UNDEFINED PID\n"); break;
  }
  cutPID->SetQAOn(qaList);

  printf("CREATE MC KINE CUTS\n");
  TObjArray* mcList = new TObjArray(0) ;
  mcList->AddLast(mcKineCuts);
  mcList->AddLast(mcGenCuts);

  printf("CREATE ACCEPTANCE CUTS\n");
  TObjArray* accList = new TObjArray(0) ;
  accList->AddLast(mcAccCuts);

  printf("CREATE RECONSTRUCTION CUTS\n");
  TObjArray* recList = new TObjArray(0) ;
  recList->AddLast(recKineCuts);
  recList->AddLast(recQualityCuts);
  recList->AddLast(recIsPrimaryCuts);

  printf("CREATE PID CUTS\n");
  TObjArray* fPIDCutList = new TObjArray(0) ;
  fPIDCutList->AddLast(cutPID);

  //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* man = new AliCFManager() ;
  man->SetParticleContainer     (container);
  man->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList);
  man->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  man->SetParticleCutsList(AliCFManager::kPartRecCuts,recList);
  man->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList);


  //CREATE THE TASK
  printf("CREATE TASK\n");
  // create the task
  AliCFSingleTrackTask *task = new AliCFSingleTrackTask("AliSingleTrackTask");
  task->SetCFManager(man); //here is set the CF manager
  task->SetQAList(qaList);

  //SETUP THE ANALYSIS MANAGER TO READ INPUT CHAIN AND WRITE DESIRED OUTPUTS
  printf("CREATE ANALYSIS MANAGER\n");
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  AliESDInputHandler* esdHandler = new AliESDInputHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  mgr->SetInputEventHandler(esdHandler);

  // Create and connect containers for input/output

  //------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->CreateContainer("cchain0",TChain::Class(),AliAnalysisManager::kInputContainer);

  // ----- output data -----
  
  //slot 0 : default output tree (by default handled by AliAnalysisTaskSE)
  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("ctree0", TTree::Class(),AliAnalysisManager::kOutputContainer,"output.root");

  //now comes user's output objects :
  
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TH1I::Class(),AliAnalysisManager::kOutputContainer,"output.root");
  // output Correction Framework Container (for acceptance & efficiency calculations)
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ccontainer0", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,"output.root");
  // output QA histograms 
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("clist0", TList::Class(),AliAnalysisManager::kOutputContainer,"output.root");

  cinput0->SetData(analysisChain);

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput0);
  mgr->ConnectOutput(task,0,coutput0);
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);

  printf("READY TO RUN\n");
  //RUN !!!
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("proof",analysisChain);
  }

  benchmark.Stop("AliSingleTrackTaskCAF");
  benchmark.Show("AliSingleTrackTaskCAF");

  return kTRUE ;
}

