//DEFINITION OF A FEW CONSTANTS
const Double_t ptmin  =   0.0 ;
const Double_t ptmax  =   1.0 ;
const Double_t etamin =  -1.5 ;
const Double_t etamax =   1.5 ;
const Int_t    PDG    =    211; 
const Int_t    minclustersTPC = 40 ;
//----------------------------------------------------

/*

Macro to prepare the necessary objects to perform spectrum unfolding 
(to be used by the class AliCFUnfolding)

Note that the efficiency calculation (and therfore the container filling) 
has to be done using MC values (and not reconstructed values)

This macro creates the following objects :
- (AliCFContainer) container    : used to calculate the efficiency (see AliCFEffGrid)
- (THnSparseD)     correlation  : this is the response matrix

Once you have run this macro, you may launch unfolding using 
the example macro CORRFW/test/testUnfolding.C

*/

Bool_t AliCFTaskForUnfolding()
{
  
  TBenchmark benchmark;
  benchmark.Start("AliSingleTrackTask");

  AliLog::SetGlobalDebugLevel(0);

  Load() ; //load the required libraries

  TChain * analysisChain ;
  printf("\n\nRunning on local file, please check the path\n\n");

  analysisChain = new TChain("esdTree");
  analysisChain->Add("your_data_path/001/AliESDs.root");
  analysisChain->Add("your_data_path/002/AliESDs.root");

  
  Info("AliCFTaskForUnfolding",Form("CHAIN HAS %d ENTRIES",(Int_t)analysisChain->GetEntries()));

  //CONTAINER DEFINITION
  Info("AliCFTaskForUnfolding","SETUP CONTAINER");
  //the sensitive variables (2 in this example), their indices
  UInt_t ipt = 0;
  UInt_t ieta  = 1;
  //Setting up the container grid... 
  UInt_t nstep = 3 ; //number of selection steps MC 
  const Int_t nvar   = 2 ; //number of variables on the grid:pt,eta
  const Int_t nbin[nvar] = {20,20} ;

  //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbin[0];
  iBin[1]=nbin[1];

  //arrays for lower bounds :
  Double_t *binLim0=new Double_t[nbin[0]+1];
  Double_t *binLim1=new Double_t[nbin[1]+1];

  //values for bin lower bounds
  for(Int_t i=0; i<=nbin[0]; i++) binLim0[i]=(Double_t)ptmin  + ( ptmax- ptmin)/nbin[0]*(Double_t)i ; 
  for(Int_t i=0; i<=nbin[1]; i++) binLim1[i]=(Double_t)etamin + (etamax-etamin)/nbin[1]*(Double_t)i ; 

  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("container","container for tracks",nstep,nvar,iBin);
  //setting the bin limits
  container -> SetBinLimits(ipt,binLim0);
  container -> SetBinLimits(ieta ,binLim1);
  container -> SetVarTitle(ipt,"pT");
  container -> SetVarTitle(ieta,"#eta");

  // Gen-Level kinematic cuts
  AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
  mcKineCuts->SetPtRange (ptmin ,ptmax );
  mcKineCuts->SetEtaRange(etamin,etamax);

  //Particle-Level cuts:  
  AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
  mcGenCuts->SetRequireIsPrimary();
  mcGenCuts->SetRequirePdgCode(PDG);

  // Rec-Level kinematic cuts
  AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts","rec-level kine cuts");
  recKineCuts->SetPtRange( ptmin, ptmax);
  recKineCuts->SetPtRange(etamin,etamax);

  AliCFTrackQualityCuts *recQualityCuts = new AliCFTrackQualityCuts("recQualityCuts","rec-level quality cuts");
  recQualityCuts->SetMinNClusterTPC(minclustersTPC);

  AliCFTrackIsPrimaryCuts *recIsPrimaryCuts = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts","rec-level isPrimary cuts");
  recIsPrimaryCuts->SetMaxNSigmaToVertex(3);

  printf("CREATE MC KINE CUTS\n");
  TObjArray* mcList = new TObjArray(0) ;
  mcList->AddLast(mcKineCuts);
  mcList->AddLast(mcGenCuts);

  printf("CREATE RECONSTRUCTION CUTS\n");
  TObjArray* recList = new TObjArray(0) ;
  recList->AddLast(recKineCuts);
  recList->AddLast(recQualityCuts);
  recList->AddLast(recIsPrimaryCuts);

  TObjArray* emptyList = new TObjArray(0);

  //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* man = new AliCFManager() ;
  man->SetNStepEvent(0);
  man->SetParticleContainer(container);
  man->SetParticleCutsList(0,mcList);
  man->SetParticleCutsList(1,recList);
  man->SetParticleCutsList(2,emptyList); // this is special step for monte carlo spectrum

  //CREATE THE TASK
  printf("CREATE TASK\n");
  // create the task
  AliCFTaskForUnfolding *task = new AliCFTaskForUnfolding("AliCFTaskForUnfolding");
  task->SetCFManager(man); //here is set the CF manager

  //create correlation matrix for unfolding
  Int_t thnDim[2*nvar];
  for (int k=0; k<nvar; k++) {
    //first half  : reconstructed 
    //second half : MC
    thnDim[k]      = nbin[k];
    thnDim[k+nvar] = nbin[k];
  }
  THnSparseD* correlation = new THnSparseD("correlation","THnSparse with correlations",2*nvar,thnDim);
  Double_t** binEdges = new Double_t[nvar];
  for (int k=0; k<nvar; k++) {
    binEdges[k]=new Double_t[nbin[k]+1];
    container->GetBinLimits(k,binEdges[k]);
    correlation->SetBinEdges(k,binEdges[k]);
    correlation->SetBinEdges(k+nvar,binEdges[k]);
  }
  correlation->Sumw2();
  task->SetCorrelationMatrix(correlation);
  //---

  //SETUP THE ANALYSIS MANAGER TO READ INPUT CHAIN AND WRITE DESIRED OUTPUTS
  printf("CREATE ANALYSIS MANAGER\n");
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  mgr->SetAnalysisType(AliAnalysisManager::kLocalAnalysis);

  AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
 
  AliInputEventHandler* dataHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(dataHandler);

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
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("corr0", THnSparseD::Class(),AliAnalysisManager::kOutputContainer,"output.root");

  cinput0->SetData(analysisChain);

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,0,coutput0);
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
 
  printf("READY TO RUN\n");
  //RUN !!!
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",analysisChain);
  }

  benchmark.Stop("AliSingleTrackTask");
  benchmark.Show("AliSingleTrackTask");

  return kTRUE ;
}

void Load() {

  //load the required aliroot libraries
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;
  gSystem->Load("libCORRFW") ;

  //compile online the task class
  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ROOTSYS/include");
  gROOT->LoadMacro("./AliCFTaskForUnfolding.cxx+");
}
