//DEFINITION OF A FEW CONSTANTS
const Double_t   ymin  = -1.0 ;
const Double_t   ymax  =  1.0 ;
const Double_t   ptmin =  0.0 ;
const Double_t   ptmax =  8.0 ;
const Int_t      mintrackrefsTPC = 2 ;
const Int_t      mintrackrefsITS = 3 ;
const Int_t      charge  = 0 ;
const Int_t      PDG = 313; 
const Int_t      minclustersTPC = 50 ;
const Double32_t nsigmavtx = 3. ;       //max track sigma to PVertex
//----------------------------------------------------

Bool_t AliCFRsnTask(
		    const Bool_t useGrid = 1,
		    const Bool_t readAOD = 0,
		    const char * kTagXMLFile="wn.xml", // XML file containing tags
		    )
{
  
  TBenchmark benchmark;
  benchmark.Start("AliRsnTask");

  AliLog::SetGlobalDebugLevel(0);

  Load(useGrid) ; //load the required libraries

  TChain * analysisChain ;

  if (useGrid) { //data located on AliEn
    TGrid::Connect("alien://") ;    //  Create an AliRunTagCuts and an AliEventTagCuts Object and impose some selection criteria
    AliRunTagCuts      *runCuts   = new AliRunTagCuts(); 
    AliEventTagCuts    *eventCuts = new AliEventTagCuts(); 
    AliLHCTagCuts      *lhcCuts   = new AliLHCTagCuts(); 
    AliDetectorTagCuts *detCuts   = new AliDetectorTagCuts(); 
    eventCuts->SetMultiplicityRange(0,20000);
    //  Create an AliTagAnalysis Object and chain the tags
    AliTagAnalysis   *tagAna = new AliTagAnalysis(); 
    tagAna->SetType("ESD");  //for aliroot > v4-05
    TGridCollection *coll   = gGrid->OpenCollection(kTagXMLFile);
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
      analysisChain->Add("your_data_path/001/AliAOD.root");
      analysisChain->Add("your_data_path/002/AliAOD.root");
    }
    else {
      analysisChain = new TChain("esdTree");
      analysisChain->Add("your_data_path/001/AliESDs.root");
      analysisChain->Add("your_data_path/002/AliESDs.root");
    }
  }
  
  
  Info("AliCFRsnTask",Form("CHAIN HAS %d ENTRIES",(Int_t)analysisChain->GetEntries()));

  //CONTAINER DEFINITION
  Info("AliCFRsnTask","SETUP CONTAINER");
  //the sensitive variables (2 in this example), their indices
  Int_t ipt = 0;
  Int_t iy  = 1;
  //Setting up the container grid... 
  Int_t nstep = 4 ; //number of selection steps MC 
  const Int_t nvar   = 2 ; //number of variables on the grid:pt,y
  const Int_t nbin1  = 8 ; //bins in pt
  const Int_t nbin2  = 8 ; //bins in y 
  //arrays for the number of bins in each dimension
  const Int_t iBin[nvar] ={nbin1,nbin2};
  //arrays for lower bounds :
  Double_t binLim1[nbin1+1];
  Double_t binLim2[nbin2+1];
  //values for bin lower bounds
  for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)ptmin + (ptmax-ptmin)/nbin1*(Double_t)i ; 
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)ymin  + (ymax-ymin)  /nbin2*(Double_t)i ;
  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("container","container for tracks",nstep,nvar,iBin);
  //setting the bin limits
  container -> SetBinLimits(ipt,binLim1);
  container -> SetBinLimits(iy,binLim2);


  //CREATE THE  CUTS -----------------------------------------------
  //Particle-Level cuts:  

  // Gen-Level kinematic cuts
  AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
  mcKineCuts->SetPtRange(ptmin,ptmax);
  mcKineCuts->SetRapidityRange(ymin,ymax);
  mcKineCuts->SetChargeMC(charge);

  AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
  //  mcGenCuts->SetRequireIsPrimary();
  mcGenCuts->SetRequirePdgCode(PDG);

  //Acceptance Cuts
  AliCFPairAcceptanceCuts *mcAccCuts = new AliCFPairAcceptanceCuts("mcAccCuts","MC acceptance cuts");
  mcAccCuts->GetNegCut()->SetMinNHitITS(mintrackrefsITS);
  mcAccCuts->GetPosCut()->SetMinNHitITS(mintrackrefsITS);
  mcAccCuts->GetNegCut()->SetMinNHitTPC(mintrackrefsTPC);
  mcAccCuts->GetPosCut()->SetMinNHitTPC(mintrackrefsTPC);

  // Rec-Level kinematic cuts
  AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts","rec-level kine cuts");
  recKineCuts->SetPtRange(ptmin,ptmax);
  recKineCuts->SetRapidityRange(ymin,ymax);
  recKineCuts->SetChargeRec(charge);

  AliCFPairQualityCuts *recQualityCuts = new AliCFPairQualityCuts("recQualityCuts","rec-level quality cuts");
  if (!readAOD) {
    recQualityCuts->GetNegCut()->SetMinNClusterTPC(minclustersTPC);
    recQualityCuts->GetPosCut()->SetMinNClusterTPC(minclustersTPC);
  }
  recQualityCuts->GetNegCut()->SetStatus(AliESDtrack::kTPCrefit & AliESDtrack::kITSrefit);
  recQualityCuts->GetPosCut()->SetStatus(AliESDtrack::kTPCrefit & AliESDtrack::kITSrefit);


  AliCFPairIsPrimaryCuts *recIsPrimaryCuts = new AliCFPairIsPrimaryCuts("recIsPrimaryCuts","rec-level isPrimary cuts");
  if (readAOD) {
    recIsPrimaryCuts->GetNegCut()->SetAODType(AliAODTrack::kPrimary);
    recIsPrimaryCuts->GetPosCut()->SetAODType(AliAODTrack::kPrimary);
  }
  else {
    recIsPrimaryCuts->GetNegCut()->SetMaxNSigmaToVertex(nsigmavtx);
    recIsPrimaryCuts->GetPosCut()->SetMaxNSigmaToVertex(nsigmavtx);
  }

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


  cutPID->GetNegCut()->SetPriors(prior_pbpb);
  cutPID->GetPosCut()->SetPriors(prior_pbpb);
  cutPID->GetNegCut()->SetProbabilityCut(0);
  cutPID->GetPosCut()->SetProbabilityCut(0);
  cutPID->GetNegCut()->SetDetectors("ITS TPC TRD");
  cutPID->GetPosCut()->SetDetectors("ITS TPC TRD");
  if (readAOD) {
    cutPID->GetNegCut()->SetAODmode(kTRUE);
    cutPID->GetPosCut()->SetAODmode(kTRUE);
  }
  else {
    cutPID->GetNegCut()->SetAODmode(kFALSE);
    cutPID->GetPosCut()->SetAODmode(kFALSE);
  }
  
  switch(PDG) {
  case  -313  : 
    cutPID->GetNegCut()->SetParticleType(AliPID::kKaon  ,kTRUE);
    cutPID->GetPosCut()->SetParticleType(AliPID::kPion  ,kTRUE);
    break;
  case   313  : 
    cutPID->GetNegCut()->SetParticleType(AliPID::kPion  ,kTRUE);
    cutPID->GetPosCut()->SetParticleType(AliPID::kKaon  ,kTRUE); 
    break;
  case   333  : 
    cutPID->GetNegCut()->SetParticleType(AliPID::kKaon  ,kTRUE);
    cutPID->GetPosCut()->SetParticleType(AliPID::kKaon  ,kTRUE); 
    break;
  case  3124  : 
    cutPID->GetNegCut()->SetParticleType(AliPID::kKaon  ,kTRUE);
    cutPID->GetPosCut()->SetParticleType(AliPID::kProton,kTRUE); 
    break;
  case -3124  : 
    cutPID->GetNegCut()->SetParticleType(AliPID::kProton,kTRUE);
    cutPID->GetPosCut()->SetParticleType(AliPID::kKaon  ,kTRUE); 
    break;
  default     : printf("UNDEFINED PID\n"); break;
  }
  //cutPID->SetQAOn(kTRUE);

  Info("AliCFRsnTask","CREATE MC KINE CUTS");
  TObjArray* mcList = new TObjArray(0) ;
  mcList->AddLast(mcKineCuts);
  mcList->AddLast(mcGenCuts);

  Info("AliCFRsnTask","CREATE ACCEPTANCE CUTS");
  TObjArray* accList = new TObjArray(0) ;
  accList->AddLast(mcAccCuts);

  Info("AliCFRsnTask","CREATE RECONSTRUCTION CUTS");
  TObjArray* recList = new TObjArray(0) ;
  recList->AddLast(recKineCuts);
  recList->AddLast(recQualityCuts);
  recList->AddLast(recIsPrimaryCuts);

  Info("AliCFRsnTask","CREATE PID CUTS");
  TObjArray* fPIDCutList = new TObjArray(0) ;
  fPIDCutList->AddLast(cutPID);


  //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  Info("AliCFRsnTask","CREATE INTERFACE AND CUTS");
  AliCFManager* man = new AliCFManager() ;
  man->SetParticleContainer     (container);
  man->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList);
  man->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  man->SetParticleCutsList(AliCFManager::kPartRecCuts,recList);
  man->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList);


  //CREATE THE TASK
  Info("AliCFRsnTask","CREATE TASK");
  // create the task
  AliCFRsnTask *task = new AliCFRsnTask("AliRsnTask");
  task->SetCFManager(man); //here is set the CF manager
  task->SetRsnPDG(PDG);


  //SETUP THE ANALYSIS MANAGER TO READ INPUT CHAIN AND WRITE DESIRED OUTPUTS
  Info("AliCFRsnTask","CREATE ANALYSIS MANAGER");
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  if (useGrid) mgr->SetAnalysisType(AliAnalysisManager::kGridAnalysis);
  else mgr->SetAnalysisType(AliAnalysisManager::kLocalAnalysis);

  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);

  AliInputEventHandler* dataHandler ;
  if (readAOD) dataHandler = new AliAODInputHandler();
  else         dataHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(dataHandler);
  


  // Create and connect containers for input/output

  //input data
  AliAnalysisDataContainer *cinput0  = mgr->CreateContainer("cchain0",TChain::Class(),AliAnalysisManager::kInputContainer);

  //slot 0 : default output tree (by default handled by AliAnalysisTaskSE)
  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("ctree0", TTree::Class(),AliAnalysisManager::kOutputContainer,"output_rsn.root");

  // output histo (number of events processed)
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TH1I::Class(),AliAnalysisManager::kOutputContainer,"output_rsn.root");
  // output Correction Framework Container (for acceptance & efficiency calculations)
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ccontainer0", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,"output_rsn.root");

  cinput0->SetData(analysisChain);

  mgr->AddTask(task);
  mgr->ConnectInput (task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,0,coutput0);
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
 

  Info("AliCFRsnTask","READY TO RUN");
  //RUN !!!
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",analysisChain);
  }

  benchmark.Stop("AliRsnTask");
  benchmark.Show("AliRsnTask");
  
  return kTRUE ;
}

void Load(Bool_t useGrid) {
  //load the required aliroot libraries
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;
  gSystem->Load("libPWG2resonances");
  gSystem->Load("libCORRFW");

  gSystem->AddIncludePath("-I$ALICE_ROOT/PWG2/RESONANCES");
  gROOT->LoadMacro("./AliCFRsnTask.cxx+");
}
