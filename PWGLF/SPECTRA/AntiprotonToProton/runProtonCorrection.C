Bool_t runProtonCorrection(Int_t stats = 0, const char* dataset = 0x0) {
  //macro used to extract the correction maps
  //using the official correction framework of ALICE
  //for protons and antiprotons
  //Author: Panos.Christakoglou@cern.ch

  //________________________________________//
  //Connect to proof  
  //TProof::Reset("alicecaf.cern.ch");
  TProof::Open("alicecaf.cern.ch");

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
  gProof->UploadPackage("CORRFW.par");
  gProof->EnablePackage("CORRFW");

  gProof->Load("./AliProtonCorrectionTask.cxx+g");

  //________________________________________//
  //Container definition
  //Variables of the GRID
  //For the time being: y-pT
  //Next step: add Vz
  //Setting up the container grid...
  //===============//
  //Global tracking//
  //===============//
  const Double_t ymin  = -1.0;
  const Double_t ymax  =  1.0;
  const Int_t nbin1 = 20; //bins in y
  const Double_t ptmin =  0.4;
  const Double_t ptmax =  3.1;
  const Int_t nbin2 = 27; //bins in pT 
  //===============//
  //  TPC tracking //
  //===============//
  /*const Double_t ymin  = -0.5;
  const Double_t ymax  =  0.5;
  const Int_t nbin1 = 10; //bins in y
  const Double_t ptmin =  0.4;
  const Double_t ptmax =  0.9;
  const Int_t nbin2 = 15;*/ //bins in pT 
  
  //Setting up the container grid... 
  UInt_t nstep = 4; //number of selection steps MC 
  UInt_t iy  = 0;
  UInt_t ipT = 1;
  const Int_t nvar = 2; //number of variables on the grid:y-pT
    //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0]=nbin1;
  iBin[1]=nbin2;
  //arrays for lower bounds :
  Double_t *binLim1=new Double_t[nbin1+1];
  Double_t *binLim2=new Double_t[nbin2+1];
  //values for bin lower bounds
  for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)ymin  + (ymax-ymin)  /nbin1*(Double_t)i;
  for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)ptmin + (ptmax-ptmin)/nbin2*(Double_t)i; 
  //CF container for protons
  AliCFContainer* containerProtons = new AliCFContainer("containerProtons","container for protons",nstep,nvar,iBin);
  //setting the bin limits
  containerProtons->SetBinLimits(iy,binLim1);
  containerProtons->SetBinLimits(ipT,binLim2);
 //CF container for antiprotons
  AliCFContainer* containerAntiProtons = new AliCFContainer("containerAntiProtons","container for antiprotons",nstep,nvar,iBin);
  //setting the bin limits
  containerAntiProtons->SetBinLimits(iy,binLim1);
  containerAntiProtons->SetBinLimits(ipT,binLim2);

  //________________________________________//
  // SET TLIST FOR QA HISTOS
  TList* qaList = new TList();
  //Cuts
  const Int_t    mintrackrefsTPC = 2;
  const Int_t    mintrackrefsITS = 3;
  const Int_t    chargeProtons  = 1;
  const Int_t    PDGProtons = 2212; 
  const Int_t    chargeAntiProtons  = -1;
  const Int_t    PDGAntiProtons = -2212; 

  const Int_t    minclustersTPC = 50;
  const Float_t  maxChi2PerTPCCluster = 3.5;
  const Float_t  maxCov11 = 2.0;
  const Float_t  maxCov22 = 2.0;
  const Float_t  maxCov33 = 0.5;
  const Float_t  maxCov44 = 0.5;
  const Float_t  maxCov55 = 2.0;
  const Float_t  maxSigmaToVertexTPC = 2.5;

  const Int_t    minclustersITS = 5;
  const Float_t  maxSigmaToVertex = 2.5;

  // Gen-Level kinematic cuts
  AliCFTrackKineCuts *mcKineCutsProtons = new AliCFTrackKineCuts("mcKineCutsProtons",
								 "MC-level kinematic cuts");
  mcKineCutsProtons->SetPtRange(ptmin,ptmax);
  mcKineCutsProtons->SetRapidityRange(ymin,ymax);
  mcKineCutsProtons->SetChargeMC(chargeProtons);
  mcKineCutsProtons->SetQAOn(qaList);

  AliCFTrackKineCuts *mcKineCutsAntiProtons = new AliCFTrackKineCuts("mcKineCutsAntiProtons",
								     "MC-level kinematic cuts");
  mcKineCutsAntiProtons->SetPtRange(ptmin,ptmax);
  mcKineCutsAntiProtons->SetRapidityRange(ymin,ymax);
  mcKineCutsAntiProtons->SetChargeMC(chargeAntiProtons);
  mcKineCutsAntiProtons->SetQAOn(qaList);

  //Particle-Level cuts:  
  AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts",
							     "MC particle generation cuts");
  mcGenCuts->SetRequireIsPrimary();
  mcGenCuts->SetRequirePdgCode(PDGProtons);
  mcGenCuts->SetQAOn(qaList);

  //Acceptance Cuts
  AliCFAcceptanceCuts *mcAccCuts = new AliCFAcceptanceCuts("mcAccCuts",
							   "MC acceptance cuts");
  mcAccCuts->SetMinNHitITS(mintrackrefsITS);
  mcAccCuts->SetMinNHitTPC(mintrackrefsTPC);
  mcAccCuts->SetQAOn(qaList);

  // Rec-Level kinematic cuts
  AliCFTrackKineCuts *recKineCutsProtons = new AliCFTrackKineCuts("recKineCutsProtons",
								  "rec-level kine cuts");
  recKineCutsProtons->SetPtRange(ptmin,ptmax);
  recKineCutsProtons->SetRapidityRange(ymin,ymax);
  recKineCutsProtons->SetChargeRec(chargeProtons);
  recKineCutsProtons->SetQAOn(qaList);

  AliCFTrackKineCuts *recKineCutsAntiProtons = new AliCFTrackKineCuts("recKineCutsAntiProtons",
								      "rec-level kine cuts");
  recKineCutsAntiProtons->SetPtRange(ptmin,ptmax);
  recKineCutsAntiProtons->SetRapidityRange(ymin,ymax);
  recKineCutsAntiProtons->SetChargeRec(chargeAntiProtons);
  recKineCutsAntiProtons->SetQAOn(qaList);

  AliCFTrackQualityCuts *recQualityCuts = new AliCFTrackQualityCuts("recQualityCuts",
								    "rec-level quality cuts");
  recQualityCuts->SetMinNClusterTPC(minclustersTPC);
  recQualityCuts->SetMaxChi2PerClusterTPC(maxChi2PerTPCCluster);
  recQualityCuts->SetMaxCovDiagonalElements(maxCov11,maxCov22,maxCov33,maxCov44,maxCov55);
  recQualityCuts->SetRequireTPCRefit(kTRUE);
  
  //recQualityCuts->SetRequireITSRefit(kTRUE);
  recQualityCuts->SetQAOn(qaList);

  AliCFTrackIsPrimaryCuts *recIsPrimaryCuts = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts",
									  "rec-level isPrimary cuts");
  recIsPrimaryCuts->SetMaxNSigmaToVertex(maxSigmaToVertex);
  recIsPrimaryCuts->SetQAOn(qaList);

  AliCFTrackCutPid* cutPID = new AliCFTrackCutPid("cutPID","ESD_PID");
  int n_species = AliPID::kSPECIES;
  Double_t* prior = new Double_t[n_species];
  
  prior[0] = 0.0244519;
  prior[1] = 0.0143988;
  prior[2] = 0.805747 ;
  prior[3] = 0.0928785;
  prior[4] = 0.0625243;
  
  cutPID->SetPriors(prior);
  cutPID->SetProbabilityCut(0.0);
  cutPID->SetDetectors("ITS TPC TOF");
  cutPID->SetParticleType(AliPID::kProton  , kTRUE); 
  cutPID->SetQAOn(qaList);

  //________________________________________// 
  TObjArray* mcListProtons = new TObjArray(0);
  mcListProtons->AddLast(mcKineCutsProtons);
  mcListProtons->AddLast(mcGenCuts);
  TObjArray* mcListAntiProtons = new TObjArray(0);
  mcListAntiProtons->AddLast(mcKineCutsAntiProtons);
  mcListAntiProtons->AddLast(mcGenCuts);

  printf("CREATE ACCEPTANCE CUTS\n");
  TObjArray* accList = new TObjArray(0);
  accList->AddLast(mcAccCuts);

  printf("CREATE RECONSTRUCTION CUTS\n");
  TObjArray* recListProtons = new TObjArray(0);
  recListProtons->AddLast(recKineCutsProtons);
  recListProtons->AddLast(recQualityCuts);
  recListProtons->AddLast(recIsPrimaryCuts);
  TObjArray* recListAntiProtons = new TObjArray(0);
  recListAntiProtons->AddLast(recKineCutsAntiProtons);
  recListAntiProtons->AddLast(recQualityCuts);
  recListAntiProtons->AddLast(recIsPrimaryCuts);

  printf("CREATE PID CUTS\n");
  TObjArray* fPIDCutList = new TObjArray(0);
  fPIDCutList->AddLast(cutPID);

  //________________________________________// 
  //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
  AliCFManager* manProtons = new AliCFManager();
  manProtons->SetParticleContainer(containerProtons);
  manProtons->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListProtons);
  manProtons->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  manProtons->SetParticleCutsList(AliCFManager::kPartRecCuts,recListProtons);
  manProtons->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList);

  AliCFManager* manAntiProtons = new AliCFManager();
  manAntiProtons->SetParticleContainer(containerAntiProtons);
  manAntiProtons->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListAntiProtons);
  manAntiProtons->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  manAntiProtons->SetParticleCutsList(AliCFManager::kPartRecCuts,recListAntiProtons);
  manAntiProtons->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList);

  //________________________________________// 
  //CREATE THE TASK
  AliProtonCorrectionTask *task = new AliProtonCorrectionTask("AliProtonCorrectionTask");
  task->SetCFManagerProtons(manProtons); //here is set the CF manager
  task->SetCFManagerAntiProtons(manAntiProtons); //here is set the CF manager
  task->SetQAList(qaList);

  //SETUP THE ANALYSIS MANAGER TO READ INPUT CHAIN AND WRITE DESIRED OUTPUTS
  printf("CREATE ANALYSIS MANAGER\n");
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  AliESDInputHandler* esdHandler = new AliESDInputHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  mgr->SetInputEventHandler(esdHandler);

  //------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();

  // ----- output data -----
  //slot 0 : default output tree (by default handled by AliAnalysisTaskSE)
  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("ctree0", TTree::Class(),
							    AliAnalysisManager::kOutputContainer,"corrections.root");
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TH1I::Class(),
							    AliAnalysisManager::kOutputContainer,"corrections.root");
  // output Correction Framework Container (for acceptance & efficiency calculations)
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ccontainer0", 
							    AliCFContainer::Class(),
							    AliAnalysisManager::kOutputContainer,"corrections.root");
  // output Correction Framework Container (for acceptance & efficiency calculations)
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("ccontainer1", 
							    AliCFContainer::Class(),
							    AliAnalysisManager::kOutputContainer,"corrections.root");
  // output QA histograms 
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("clist0", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,"corrections.root");
  
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput0);
  mgr->ConnectOutput(task,0,coutput0);
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);

  //________________________________________// 
  if (mgr->InitAnalysis()) {
    if(dataset)
      mgr->StartAnalysis("proof",dataset,stats);
    else {
      // You should get this macro and the txt file from:
      // http://aliceinfo.cern.ch/Offline/Analysis/CAF/
      gROOT->LoadMacro("CreateESDChain.C");
      TChain* chain = 0x0;
      chain = CreateESDChain("ESD82XX_30K.txt",stats);
      chain->SetBranchStatus("*Calo*",0);

      mgr->StartAnalysis("proof",chain);
    }
  }

  return kTRUE;
}

