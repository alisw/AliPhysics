// definition of variables ranges
const Double_t nevtmin 	= 0;			
const Double_t nevtmax 	= 5;

const Double_t trigmin 	= 0;
const Double_t trigmax 	= 4;
const Int_t trigsidemin = 0;	
const Int_t trigsidemax = 4;
const Int_t chargemin		= -3;	
const Int_t chargemax		= 3;

const Double_t ymin  		= -4.5; 
const Double_t ymax  		= -2.0;
const Double_t etamin		= -4.5; 
const Double_t etamax		= -2.0;
const Double_t rabsmin 	= 15;	
const Double_t rabsmax 	= 90;
const Double_t ptmin 		= 0;	
const Double_t ptmax 		= 100;
const Double_t thabsmin = 170;
const Double_t thabsmax = 180;
const Double_t vzmin  = -50;	
const Double_t vzmax  = 50;		
const Double_t dcamin  = 0;		
const Double_t dcamax  = 100;		
const Double_t pmin 		= 0;			
const Double_t pmax 		= 100;	
const Double_t mmin 		= 0;	
const Double_t mmax 		= 15;		

// Resonances
const Int_t PDG = 553;	// Upsilon

///// Setting up the container grid

// CONTAINER DEFINITION
UInt_t nevent			= 0;					// number of event
UInt_t y  				= 1;					// dimuon rapidity
UInt_t pt 				= 2;					// dimuon pt
UInt_t imass  		= 3;					// dimuon invariant mass
UInt_t trig  			= 4;					// single muon track-trigger matching
UInt_t ptmu		 		= 5;					// single muon pt
UInt_t pmu 				= 6;				  // single muon p
UInt_t trigside 	= 7;					// event trigger side (AC, B, E)
UInt_t rabsmu		 	= 8;					// single muon Rabs
UInt_t charge		  = 9;					// total charge of dimuon tracks
UInt_t etamu			= 10;					// single muon eta
UInt_t thabsmu		= 11;					// single muon theta_abs
UInt_t vzmu				= 12;					// single muon Vz
UInt_t dcamu			= 13;					// single muon DCA



// Setting up the container grid
UInt_t nstep = 5 ; //number of selection steps : MC and ESD for simulation, CINT1B, CMUS1B for real data

const Int_t nvar   	= 14;     	//number of variables on the grid
const Int_t nbin1		= 5;	
const Int_t nbin2	  = 15;  	
const Int_t nbin3	  = 100;		
const Int_t nbin4	  = 300; 	
const Int_t nbin5	  = 40;	
const Int_t nbin6	  = 100;  	
const Int_t nbin7	  = 100; 	
const Int_t nbin8	  = 4; 	
const Int_t nbin9 	= 300;
const Int_t nbin10	= 6;		
const Int_t nbin11	= 100;	
const Int_t nbin12	= 200; 	
const Int_t nbin13 	= 100;	
const Int_t nbin14	= 100;	

class AliAnalysisGrid;
class TChain;

Bool_t AliCFMuonResUpsilon(
	const char* runtype = "local",						// local, proof or grid 
	const char* taskname = "CFMuonResUpsilon",	// task name
	const bool readAOD = 0,									// 1 = AOD based analysis
	const bool readMC = 0,										// 1 = Read MC
	const bool usePlugin = 0,								// grid mode : use plugin
	const char* gridmode = "test"						// plugin mode ("full","test","offline","submit" or "terminate"
	)
{

  //load the required aliroot libraries
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD") ;
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;

	//load correction framework library
  gSystem->Load("libCORRFW.so") ;

  //compile online the task class
  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ROOTSYS/include -I$ALICE_ROOT/PWG3/muon");
  gROOT->LoadMacro("./AliCFMuonResUpsilon.cxx++g");


	// Setting Container
	// -----------------------------------------------

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

	// values for bin lower bounds
	for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)nevtmin+(nevtmax-nevtmin)/nbin1*(Double_t)i ;
	for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)ymin+(ymax-ymin)/nbin2*(Double_t)i ;
	for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)ptmin+(ptmax-ptmin)/nbin3*(Double_t)i ; 
	for(Int_t i=0; i<=nbin4; i++) binLim4[i]=(Double_t)mmin+(mmax-mmin)/nbin4*(Double_t)i ;
	for(Int_t i=0; i<=nbin5; i++) binLim5[i]=(Double_t)trigmin+(trigmax-trigmin)/nbin5*(Double_t)i ; 
	for(Int_t i=0; i<=nbin6; i++) binLim6[i]=(Double_t)ptmin+(ptmax-ptmin)/nbin6*(Double_t)i ;
	for(Int_t i=0; i<=nbin7; i++) binLim7[i]=(Double_t)pmin+(pmax-pmin)/nbin7*(Double_t)i ;
	for(Int_t i=0; i<=nbin8; i++) binLim8[i]=(Double_t)trigsidemin+(trigsidemax-trigsidemin)/nbin8*(Double_t)i ; 
	for(Int_t i=0; i<=nbin9; i++) binLim9[i]=(Double_t)rabsmin+(rabsmax-rabsmin)/nbin9*(Double_t)i ; 
	for(Int_t i=0; i<=nbin10; i++) binLim10[i]=(Double_t)chargemin+(chargemax-chargemin)/nbin10*(Double_t)i ; 
	for(Int_t i=0; i<=nbin11; i++) binLim11[i]=(Double_t)etamin+(etamax-etamin)/nbin11*(Double_t)i ; 
	for(Int_t i=0; i<=nbin12; i++) binLim12[i]=(Double_t)thabsmin+(thabsmax-thabsmin)/nbin12*(Double_t)i ;
	for(Int_t i=0; i<=nbin13; i++) binLim13[i]=(Double_t)vzmin+(vzmax-vzmin)/nbin13*(Double_t)i ;
	for(Int_t i=0; i<=nbin14; i++) binLim14[i]=(Double_t)dcamin+(dcamax-dcamin)/nbin14*(Double_t)i ;
	
	// one container  of 2 steps (MC and ESD) with 12 variables
	AliCFContainer* container = new AliCFContainer("container","container for tracks",nstep,nvar,iBin);
	// setting the bin limits
	container -> SetBinLimits(nevent,binLim1);
	container -> SetBinLimits(y,binLim2);
	container -> SetBinLimits(pt,binLim3);
	container -> SetBinLimits(imass,binLim4);
	container -> SetBinLimits(trig,binLim5);
	container -> SetBinLimits(ptmu,binLim6);
	container -> SetBinLimits(pmu,binLim7);
	container -> SetBinLimits(trigside,binLim8);
	container -> SetBinLimits(rabsmu,binLim9);
	container -> SetBinLimits(charge,binLim10);
	container -> SetBinLimits(etamu,binLim11);
	container -> SetBinLimits(thabsmu,binLim12);
	container -> SetBinLimits(vzmu,binLim13);
	container -> SetBinLimits(dcamu,binLim14);
	
	// create AliCFManager
  AliCFManager* man = new AliCFManager();
	// set particle container
  man->SetParticleContainer(container);

	// ----------------------------------------------

	// check run type
	if(runtype != "local" && runtype != "proof" && runtype != "grid") {
		printf("Incorrect runtype! choose \"local\", \"prootf\" or \"grid\"\n");
		return;
	}

	printf("runtype : %s\n", runtype);

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager(Form("%sAnalysis",taskname));

	// set run mode
  if (runtype == "grid" && usePlugin) {
		// grid mode
		// set package for grid(API,ROOT,ALIROOT)
		const char *package[3] = {"V1.1x","v5-28-00a","v4-21-18-AN"};

 		AliAnalysisGrid *alienHandler = CreateAlienHandler(taskname,gridmode,package);  
 		if (!alienHandler) return;

 		mgr->SetGridHandler(alienHandler);
  }
  else {
		// local mode
    if (readAOD) {
			// AOD
      analysisChain = new TChain("aodTree");

			if(readMC) {
				analysisChain->Add("/home/sahn/alice/mc_prod/PDC09/AODs-100421-100436v2/000/033/AliAODs.root");		// PDC09 AOD
				//analysisChain->Add("/home/sahn/alice/mc_prod/upsilon/local/run/muonAcc/residual/AODs/AliAODs.root");		// muon acceptance production with residual alignment AOD
			}
			else {
				analysisChain->Add("/home/sahn/alice/p2_prod/LHC10g/AODs/135658/019/AliAOD.Dimuons.root");		// LHC10g-run135658 AOD
				//analysisChain->Add("/home/sahn/alice/p2_prod/LHC10g/AODs/135658/020/AliAOD.Dimuons.root");		// LHC10g-run135658 AOD
			}
	  }
    else {
			// ESD
			analysisChain = new TChain("esdTree");
		
			if(readMC) {
				analysisChain->Add("/home/sahn/alice/mc_prod/upsilon/local/run/muonAcc/residual/005/AliESDs.root");		// muon acceptance production with residual alignment
			}
			else {
				analysisChain->Add("/home/sahn/alice/p2_prod/LHC10g/000135748/ESDs/pass1/10000135748001.150/AliESDs.root");		// LHC10g-run135748
			}
   	}
  }
  
  // create the task
  AliCFMuonResUpsilon *task = new AliCFMuonResUpsilon(taskname);

  // Set list
  TList* qaList = new TList();

  //CREATE THE CUTS
  // Choice of the Resonance
  if(readMC && !readAOD) {	
		AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
  	mcGenCuts->SetRequirePdgCode(PDG);
  	mcGenCuts->SetQAOn(qaList);

  	// Set a pt range of the resonance
  	AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
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

  	man->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList);
  	man->SetParticleCutsList(AliCFManager::kPartAccCuts,recList);
	}
	else if(readMC && readAOD) {	// for MC in AOD
		task->SetRequirePdgCode(PDG);
		task->SetPtRange(ptmin,ptmax);
		task->SetRapidityRange(ymin,ymax);
	}

  task->SetCFManager(man); //here is set the CF manager
  task->SetQAList(qaList);

  if(readAOD) task->SetReadAODData();
	if(readMC) 	task->SetReadMCInfo();

	if(readMC && !readAOD) {
  	AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  	mgr->SetMCtruthEventHandler(mcHandler);
	}

  if(readAOD) { 
		AliAODInputHandler *dataHandler = new AliAODInputHandler(); 
  	mgr->SetInputEventHandler(dataHandler);

		if(!readMC){
		// physics selection
		gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
		AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
		}
	}
  else { 
		AliESDInputHandler *dataHandler = new AliESDInputHandler(); 
  	mgr->SetInputEventHandler(dataHandler);

		if(!readMC) {
		// physics selection
		gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
		AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
		}
	}


  // Create and connect containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->CreateContainer("cchain",TChain::Class(),AliAnalysisManager::kInputContainer);

  // output data
  Char_t file[256];
	sprintf(file,"container.%s.%s.%s.%s.root",taskname,runtype,readAOD ? "AOD":"ESD",readMC ? "MC":"Real");

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TH1D::Class(),AliAnalysisManager::kOutputContainer,file);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("container", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,file);

  mgr->AddTask(task);

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);

  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus();
  if(usePlugin) mgr->StartAnalysis(runtype);
  else mgr->StartAnalysis(runtype,analysisChain);

  return kTRUE ;
}

AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, const char* package[3])
{
  //if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
	
	// AliEn plugin mode
  plugin->SetRunMode(gridmode);

	// SE copy test
	//plugin->SetCheckCopy(kFALSE);

	// JDL Merge
	plugin->SetMergeViaJDL(kTRUE);
	plugin->SetOneStageMerging(kFALSE);
	plugin->SetMaxMergeStages(3);

	// Set package
  plugin->SetAPIVersion(package[0]);
  plugin->SetROOTVersion(package[1]);
  plugin->SetAliROOTVersion(package[2]);
  plugin->SetDataPattern("ESDs/pass1/*AliESDs.root");


	// Set run numbers

	//______________________________________________________
	// LHC10f 
  plugin->SetGridDataDir("/alice/data/2010/LHC10f/");
	const Int_t nruns=24;

	Int_t run_number[nruns] = { 133006, 133010, 133327, 133330, 133414,	133419, 133563, 133800, 133924, 133969,
															133985, 134094, 134198, 134204, 134304,	134497, 134666, 134679, 134685, 134690,
															134841, 134905, 134914, 134919 };

  plugin->SetGridWorkingDir(Form("workdir-%s",taskname));

	//______________________________________________________
	// LHC10g
	/*
  plugin->SetGridDataDir("/alice/data/2010/LHC10g/");

	const Int_t nruns=12;

	Int_t run_number[nruns] = { 135658, 135704, 135709, 135712, 135748, 135761, 135795, 136177, 136180, 136189, 136372, 136376 };

  plugin->SetGridWorkingDir(Form("workdir-%s",taskname));
	*/

	// Set Prefix for real data
	plugin->SetRunPrefix("000");

	// Adding run numbers
  for(int i=0; i<nruns; i++) {
		plugin->AddRunNumber(run_number[i]);
	}

	// output and etc.
  plugin->SetGridOutputDir("output"); 
  plugin->SetAnalysisSource("AliCFMuonResUpsilon.cxx");
  plugin->SetAdditionalLibs("AliCFMuonResUpsilon.h AliCFMuonResUpsilon.cxx");
  //plugin->SetSplitMaxInputFileNumber(20);
  plugin->SetTTL(86400);
  plugin->SetInputFormat("xml-single");
  plugin->SetPrice(1);      
	plugin->SetUseSubmitPolicy();
  return plugin;
}
