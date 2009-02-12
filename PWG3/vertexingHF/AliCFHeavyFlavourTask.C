//DEFINITION OF A FEW CONSTANTS
const Double_t ymin  = -2.0 ;
const Double_t ymax  =  2.0 ;
const Double_t ptmin =  0.0 ;
const Double_t ptmax =  10.0 ;
const Int_t    mintrackrefsTPC = 2 ;
const Int_t    mintrackrefsITS = 3 ;
const Int_t    charge  = 1 ;
const Int_t    PDG = 421; 
const Int_t    minclustersTPC = 50 ;
//----------------------------------------------------

Bool_t AliCFHeavyFlavourTask()
{
	
	TBenchmark benchmark;
	benchmark.Start("AliCFHeavyFlavourTask");
	
	AliLog::SetGlobalDebugLevel(0);
	
	//load the required libraries
	Bool_t useParFiles=kFALSE;
	gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
	LoadLibraries(useParFiles);

	

	TChain * analysisChain ;
	TChain * analysisChainFriend ;
	
	//here put your input data path
	printf("\n\nRunning on local file, please check the path\n\n");
	analysisChain = new TChain("aodTree");
	analysisChainFriend = new TChain("aodTree");
	//for (Int_t i = 0; i<99; i++){
	//TString fileName = "/Events/180002/";  // set here your own path!
	//TString aodHFFileName = "/Events/180002/";  // set here your own path!
	//fileName+=Form("%i/AliAOD.root",i);
	//printf(Form("file = %s \n", fileName.Data()));
	//aodHFFileName+=Form("%i/AliAOD.VertexingHF.root",i);
	//analysisChain->Add(fileName);
	//analysisChainFriend->Add(aodHFFileName);
	analysisChain->Add("AliAOD.root");
	analysisChainFriend->Add("AliAOD.VertexingHF.root");
	  //}
			
	analysisChain->AddFriend(analysisChainFriend);
	Info("AliCFHeavyFlavourTask",Form("CHAIN HAS %d ENTRIES",(Int_t)analysisChain->GetEntries()));
	
	//CONTAINER DEFINITION
	Info("AliCFHeavyFlavourTask","SETUP CONTAINER");
	//the sensitive variables (2 in this example), their indices
	UInt_t ipt = 0;
	UInt_t iy  = 1;
	//Setting up the container grid... 
	UInt_t nstep = 2 ; //number of selection steps MC 
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
	
	//CREATE THE  CUTS -----------------------------------------------
	
	// Gen-Level kinematic cuts
	AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
	//   mcKineCuts->SetPtRange(ptmin,ptmax);
	//   mcKineCuts->SetRapidityRange(ymin,ymax);
	//   mcKineCuts->SetChargeMC(charge);
	
	//Particle-Level cuts:  
	AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
	//mcGenCuts->SetRequireIsPrimary();
	mcGenCuts->SetRequirePdgCode(PDG);
	mcGenCuts->SetAODMC(1); //special flag for reading MC in AOD tree (important)
	
	// Rec-Level kinematic cuts
	AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts","rec-level kine cuts");
	//   recKineCuts->SetPtRange(ptmin,ptmax);
	//   recKineCuts->SetRapidityRange(ymin,ymax);
	//   recKineCuts->SetChargeRec(charge);
	
	AliCFTrackQualityCuts *recQualityCuts = new AliCFTrackQualityCuts("recQualityCuts","rec-level quality cuts");
	//recQualityCuts->SetStatus(AliESDtrack::kITSrefit);
	
	AliCFTrackIsPrimaryCuts *recIsPrimaryCuts = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts","rec-level isPrimary cuts");
	//recIsPrimaryCuts->SetAODType(AliAODTrack::kPrimary);
	
	printf("CREATE MC KINE CUTS\n");
	TObjArray* mcList = new TObjArray(0) ;
	mcList->AddLast(mcKineCuts);
	mcList->AddLast(mcGenCuts);
	
	printf("CREATE RECONSTRUCTION CUTS\n");
	TObjArray* recList = new TObjArray(0) ;
	recList->AddLast(recKineCuts);
	recList->AddLast(recQualityCuts);
	recList->AddLast(recIsPrimaryCuts);
	
	//CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
	printf("CREATE INTERFACE AND CUTS\n");
	AliCFManager* man = new AliCFManager() ;
	man->SetParticleContainer     (container);
	man->SetParticleCutsList(0 , mcList); // MC
	man->SetParticleCutsList(1 , recList); // AOD
	
	//CREATE THE TASK
	printf("CREATE TASK\n");
	// create the task
	AliCFHeavyFlavourTask *task = new AliCFHeavyFlavourTask("AliCFHeavyFlavourTask");
	task->SetCFManager(man); //here is set the CF manager
	
	//SETUP THE ANALYSIS MANAGER TO READ INPUT CHAIN AND WRITE DESIRED OUTPUTS
	printf("CREATE ANALYSIS MANAGER\n");
	// Make the analysis manager
	AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
	
	mgr->SetAnalysisType(AliAnalysisManager::kLocalAnalysis);
	
	AliAODInputHandler* dataHandler = new AliAODInputHandler();
	mgr->SetInputEventHandler(dataHandler);
	
	// Create and connect containers for input/output
	
	// ------ input data ------
	AliAnalysisDataContainer *cinput0  = mgr->CreateContainer("cchain0",TChain::Class(),AliAnalysisManager::kInputContainer);
	
	// ----- output data -----
	
	//slot 0 : default output tree (by default handled by AliAnalysisTaskSE)
	AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("ctree0", TTree::Class(),AliAnalysisManager::kOutputContainer,"outputEta.root");
	
	//now comes user's output objects :
	
	// output TH1I for event counting
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TH1I::Class(),AliAnalysisManager::kOutputContainer,"output.root");
	// output Correction Framework Container (for acceptance & efficiency calculations)
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ccontainer0", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,"output.root");
	
	cinput0->SetData(analysisChain);
	
	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput0);
	mgr->ConnectOutput(task,0,coutput0);
	mgr->ConnectOutput(task,1,coutput1);
	mgr->ConnectOutput(task,2,coutput2);
	
	printf("READY TO RUN\n");
	//RUN !!!
	if (mgr->InitAnalysis()) {
		mgr->PrintStatus();
		mgr->StartAnalysis("local",analysisChain);
	}
	
	benchmark.Stop("AliCFHeavyFlavourTask");
	benchmark.Show("AliCFHeavyFlavourTask");
	
	return kTRUE ;
}

