//DEFINITION OF A FEW CONSTANTS
const Double_t ymin  = -2.1 ;
const Double_t ymax  =  2.1 ;
const Double_t ptmin_0_4 =  0.0 ;
const Double_t ptmax_0_4 =  4.0 ;
const Double_t ptmin_4_8 =  4.0 ;
const Double_t ptmax_4_8 =  8.0 ;
const Double_t ptmin_8_10 =  8.0 ;
const Double_t ptmax_8_10 =  10.0 ;
const Double_t cosmin = -1.05;
const Double_t cosmax =  1.05;
const Double_t cTmin = 0;  // micron
const Double_t cTmax = 500;  // micron
const Double_t dcamin = 0;  // micron
const Double_t dcamax = 500;  // micron
const Double_t d0min = -1000;  // micron
const Double_t d0max = 1000;  // micron
const Double_t d0xd0min = -100000;  // micron
const Double_t d0xd0max = 100000;  // micron
const Int_t    mintrackrefsTPC = 2 ;
const Int_t    mintrackrefsITS = 3 ;
const Int_t    charge  = 1 ;
const Int_t    PDG = 421; 
const Int_t    minclustersTPC = 50 ;
// cuts
const Double_t ptmin = 0.1;
const Double_t ptmax = 9999.;
const Double_t etamin = -0.9;
const Double_t etamax = 0.9;
const Int_t    minITSClusters = 5;

//----------------------------------------------------

AliCFHeavyFlavourTaskMultiVarMultiStep *AddTaskCFMultiVarMultiStep()
{

	//CONTAINER DEFINITION
	Info("AliCFHeavyFlavourTaskMultiVarMultiStep","SETUP CONTAINER");
	//the sensitive variables (6 in this example), their indices
	UInt_t ipt = 0;
	UInt_t iy  = 1;
	UInt_t icosThetaStar  = 2;
	UInt_t ipTpi  = 3;
	UInt_t ipTk  = 4;
	UInt_t icT  = 5;
	UInt_t idca  = 6;
	UInt_t id0pi  = 7;
	UInt_t id0K  = 8;
	UInt_t id0xd0  = 9;
	UInt_t ipointing  = 10;

	//Setting up the container grid... 
	UInt_t nstep = 6; //number of selection steps: MC, Acceptance, Reco (no cuts), RecoAcceptance, RecoITSClusters (RecoAcceptance included), RecoPPR (RecoAcceptance+RecoITSCluster included) 
	const Int_t nvar   = 11 ; //number of variables on the grid:pt, y, cosThetaStar, pTpi, pTk, cT, dca, d0pi, d0K, d0xd0, cosPointingAngle 
	const Int_t nbin0_0_4  = 8 ; //bins in pt from 0 to 4 GeV
	const Int_t nbin0_4_8  = 4 ; //bins in pt from 4 to 8 GeV
	const Int_t nbin0_8_10  = 1 ; //bins in pt from 8 to 10 GeV
	const Int_t nbin1  = 42 ; //bins in y
	const Int_t nbin2  = 42 ; //bins in cosThetaStar 
	const Int_t nbin3_0_4  = 8 ; //bins in ptPi from 0 to 4 GeV
	const Int_t nbin3_4_8  = 4 ; //bins in ptPi from 4 to 8 GeV
	const Int_t nbin3_8_10  = 1 ; //bins in ptPi from 8 to 10 GeV
	const Int_t nbin4_0_4  = 8 ; //bins in ptKa from 0 to 4 GeV
	const Int_t nbin4_4_8  = 4 ; //bins in ptKa from 4 to 8 GeV
	const Int_t nbin4_8_10  = 1 ; //bins in ptKa from 8 to 10 GeV
	const Int_t nbin5  = 24 ; //bins in cT
	const Int_t nbin6  = 24 ; //bins in dca
	const Int_t nbin7  = 100 ; //bins in d0pi
	const Int_t nbin8  = 100 ; //bins in d0K
	const Int_t nbin9  = 80 ; //bins in d0xd0
	const Int_t nbin10  = 1050 ; //bins in cosPointingAngle
	
	//arrays for the number of bins in each dimension
	Int_t iBin[nvar];
 	iBin[0]=nbin0_0_4+nbin0_4_8+nbin0_8_10;
	iBin[1]=nbin1;
	iBin[2]=nbin2;
 	iBin[3]=nbin3_0_4+nbin3_4_8+nbin3_8_10;
 	iBin[4]=nbin4_0_4+nbin4_4_8+nbin4_8_10;
	iBin[5]=nbin5;
	iBin[6]=nbin6;
	iBin[7]=nbin7;
	iBin[8]=nbin8;
	iBin[9]=nbin9;
	iBin[10]=nbin10;
	
	//arrays for lower bounds :
	Double_t *binLim0=new Double_t[iBin[0]+1];
	Double_t *binLim1=new Double_t[iBin[1]+1];
	Double_t *binLim2=new Double_t[iBin[2]+1];
	Double_t *binLim3=new Double_t[iBin[3]+1];
	Double_t *binLim4=new Double_t[iBin[4]+1];
	Double_t *binLim5=new Double_t[iBin[5]+1];
	Double_t *binLim6=new Double_t[iBin[6]+1];
	Double_t *binLim7=new Double_t[iBin[7]+1];
	Double_t *binLim8=new Double_t[iBin[8]+1];
	Double_t *binLim9=new Double_t[iBin[9]+1];
	Double_t *binLim10=new Double_t[iBin[10]+1];
	
	// checking limits
	if (ptmax_0_4 != ptmin_4_8) {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","max lim 1st range != min lim 2nd range, please check!");
	}
	if (ptmax_4_8 != ptmin_8_10) {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","max lim 2nd range != min lim 3rd range, please check!");
	}

	// values for bin lower bounds
	// pt
	for(Int_t i=0; i<=nbin0_0_4; i++) binLim0[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbin0_0_4*(Double_t)i ; 
	if (binLim0[nbin0_0_4] != ptmin_4_8)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 1st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_4_8; i++) binLim0[i+nbin0_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbin0_4_8*(Double_t)i ; 
	if (binLim0[nbin0_0_4+nbin0_4_8] != ptmin_8_10)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_8_10; i++) binLim0[i+nbin0_0_4+nbin0_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbin0_8_10*(Double_t)i ; 

	// y
	for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)ymin  + (ymax-ymin)  /nbin1*(Double_t)i ;

	// cosThetaStar
	for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)cosmin  + (cosmax-cosmin)  /nbin2*(Double_t)i ;

	// ptPi
	for(Int_t i=0; i<=nbin3_0_4; i++) binLim3[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbin3_0_4*(Double_t)i ; 
	if (binLim3[nbin3_0_4] != ptmin_4_8)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for ptPi - 1st range - differs from expected!");
	}
	for(Int_t i=0; i<=nbin3_4_8; i++) binLim3[i+nbin3_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbin3_4_8*(Double_t)i ; 
	if (binLim3[nbin3_0_4+nbin3_4_8] != ptmin_8_10)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for ptPi - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin3_8_10; i++) binLim3[i+nbin3_0_4+nbin3_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbin3_8_10*(Double_t)i ; 

	// ptKa
	for(Int_t i=0; i<=nbin4_0_4; i++) binLim4[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbin4_0_4*(Double_t)i ; 
	if (binLim4[nbin4_0_4] != ptmin_4_8)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for ptKa - 1st range - differs from expected!");
	}
	for(Int_t i=0; i<=nbin4_4_8; i++) binLim4[i+nbin4_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbin4_4_8*(Double_t)i ; 
	if (binLim4[nbin4_0_4+nbin4_4_8] != ptmin_8_10)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for ptKa - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin4_8_10; i++) binLim4[i+nbin4_0_4+nbin4_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbin4_8_10*(Double_t)i ; 

	// cT
	for(Int_t i=0; i<=nbin5; i++) binLim5[i]=(Double_t)cTmin  + (cTmax-cTmin)  /nbin5*(Double_t)i ;

	// dca
	for(Int_t i=0; i<=nbin6; i++) binLim6[i]=(Double_t)dcamin  + (dcamax-dcamin)  /nbin6*(Double_t)i ;

	// d0pi
	for(Int_t i=0; i<=nbin7; i++) binLim7[i]=(Double_t)d0min  + (d0max-d0min)  /nbin7*(Double_t)i ;

	// d0K
	for(Int_t i=0; i<=nbin8; i++) binLim8[i]=(Double_t)d0min  + (d0max-d0min)  /nbin8*(Double_t)i ;

	// d0xd0
	for(Int_t i=0; i<=nbin9; i++) binLim9[i]=(Double_t)d0xd0min  + (d0xd0max-d0xd0min)  /nbin9*(Double_t)i ;

	// cosPointingAngle
	for(Int_t i=0; i<=nbin10; i++) binLim10[i]=(Double_t)cosmin  + (cosmax-cosmin)  /nbin10*(Double_t)i ;

	// debugging printings
	//Info("AliCFHeavyFlavourTaskMultiVarMultiStep","Printing lower limits for bins in pt");
	//for (Int_t i =0; i<= iBin[0]; i++){
	//	Info("AliCFHeavyFlavourTaskMultiVarMultiStep",Form("i-th bin, lower limit = %f", binLim0[i]));
	//}
	//Info("Printing lower limits for bins in ptPi");
	//for (Int_t i =0; i<= iBin[3]; i++){
	//	Info("AliCFHeavyFlavourTaskMultiVarMultiStep",Form("i-th bin, lower limit = %f", binLim3[i]));
	//}
	//Info("Printing lower limits for bins in ptKa");
	//for (Int_t i =0; i<= iBin[4]; i++){
	//	Info("AliCFHeavyFlavourTaskMultiVarMultiStep",Form("i-th bin, lower limit = %f", binLim4[i]));
	//	}

	//one "container" for MC
	AliCFContainer* container = new AliCFContainer("container","container for tracks",nstep,nvar,iBin);
	//setting the bin limits
	container -> SetBinLimits(ipt,binLim0);
	container -> SetBinLimits(iy,binLim1);
	container -> SetBinLimits(icosThetaStar,binLim2);
	container -> SetBinLimits(ipTpi,binLim3);
	container -> SetBinLimits(ipTk,binLim4);
	container -> SetBinLimits(icT,binLim5);
	container -> SetBinLimits(idca,binLim6);
	container -> SetBinLimits(id0pi,binLim7);
	container -> SetBinLimits(id0K,binLim8);
	container -> SetBinLimits(id0xd0,binLim9);
	container -> SetBinLimits(ipointing,binLim10);
	
	//CREATE THE  CUTS -----------------------------------------------
	
	// Gen-Level kinematic cuts
	AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
	//   mcKineCuts->SetPtRange(ptmin,ptmax);
	//   mcKineCuts->SetRapidityRange(ymin,ymax);
	//   mcKineCuts->SetChargeMC(charge);
	
	//Particle-Level cuts:  
	AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
	//mcGenCuts->SetRequireIsPrimary();
	mcGenCuts->SetRequirePdgCode(PDG, kTRUE);  // kTRUE set in order to include D0_bar
	mcGenCuts->SetAODMC(1); //special flag for reading MC in AOD tree (important)
	
	// Acceptance cuts:
	AliCFAcceptanceCuts* accCuts = new AliCFAcceptanceCuts("accCuts", "Acceptance cuts");
	//accCuts->SetMinNHitITS(3);
	//accCuts->SetMinNHitTPC(2);
	AliCFTrackKineCuts *kineAccCuts = new AliCFTrackKineCuts("kineAccCuts","Kine-Acceptance cuts");
	kineAccCuts->SetPtRange(ptmin,ptmax);
	kineAccCuts->SetEtaRange(etamin,etamax);

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
	
	printf("CREATE ACCEPTANCE CUTS\n");
	TObjArray* accList = new TObjArray(0) ;
	accList->AddLast(kineAccCuts);

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
	man->SetParticleCutsList(1 , accList); // Acceptance 
	man->SetParticleCutsList(2 , recList); // AOD 
	man->SetParticleCutsList(3 , recList); // AOD in Acceptance
	man->SetParticleCutsList(4 , recList); // AOD with required n. of ITS clusters
	man->SetParticleCutsList(5 , recList); // AOD Reco (PPR cuts implemented in Task)
	
	// Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
	  ::Error("AddTaskCompareHF", "No analysis manager to connect to.");
	  return NULL;
	}   
	//CREATE THE TASK
	printf("CREATE TASK\n");
	// create the task
	AliCFHeavyFlavourTaskMultiVarMultiStep *task = new AliCFHeavyFlavourTaskMultiVarMultiStep("AliCFHeavyFlavourTaskMultiVarMultiStep");
	task->SetFillFromGenerated(kFALSE);
	task->SetMinITSClusters(minITSClusters);
	task->SetCFManager(man); //here is set the CF manager
		
	// Create and connect containers for input/output
	
	// ------ input data ------
	AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
	
	// ----- output data -----
	
	//slot 0 : default output tree (by default handled by AliAnalysisTaskSE)
	AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("ctree0", TTree::Class(),AliAnalysisManager::kOutputContainer,"output.root");
	
	//now comes user's output objects :
	
	// output TH1I for event counting
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TH1I::Class(),AliAnalysisManager::kOutputContainer,"output.root");
	// output Correction Framework Container (for acceptance & efficiency calculations)
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ccontainer0", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,"output.root");
	
	mgr->AddTask(task);
	
	mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task,0,coutput0);
	mgr->ConnectOutput(task,1,coutput1);
	mgr->ConnectOutput(task,2,coutput2);
	
	return task;
}

