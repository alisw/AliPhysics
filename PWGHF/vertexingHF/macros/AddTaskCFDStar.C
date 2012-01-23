//DEFINITION OF A FEW CONSTANTS
//
// binning method from C.Zampolli
//
// general
const Double_t ymin  = -2.1 ;
const Double_t ymax  =  2.1 ;
//soft pion
const Double_t ptmin_0_1 =  0.0 ;
const Double_t ptmax_0_1 =  1.0 ;
const Double_t ptmin_1_2 =  1.0 ;
const Double_t ptmax_1_2 =  2.0 ;
const Double_t ptmin_2_10 =  2.0 ;
const Double_t ptmax_2_10 =  15.0 ;
//D0 and D0 prongs
const Double_t ptmin_0_4 =  0.0 ;
const Double_t ptmax_0_4 =  4.0 ;
const Double_t ptmin_4_8 =  4.0 ;
const Double_t ptmax_4_8 =  8.0 ;
const Double_t ptmin_8_10 =  8.0 ;
const Double_t ptmax_8_10 =  20.0 ;
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
const Double_t phimin = 0.0;    
const Int_t    mintrackrefsTPC = 2 ;
const Int_t    mintrackrefsITS = 3 ;
const Int_t    charge  = 1 ; 
const Int_t    minclustersTPC = 50 ;
// cuts
const Double_t ptmin = 0.05;
const Double_t ptmax = 9999.;
const Double_t etamin = -0.9;
const Double_t etamax = 0.9;
const Double_t zmin = -15;
const Double_t zmax = 15;
const Int_t    minITSClusters = 3;
const Int_t    minITSClustersSoft = 2;
//----------------------------------------------------

AliCFTaskForDStarAnalysis *AddTaskCFDStar()
{

	//CONTAINER DEFINITION
	Info("AliCFTaskForDStarAnalysis","SETUP CONTAINER");
	//the sensitive variables, their indices
	UInt_t ipt = 0;
	UInt_t iy  = 1;
	UInt_t icosThetaStar  = 2;
	UInt_t ipTpi  = 3;
	UInt_t ipTD0  = 4;
	UInt_t icT    = 5;
	UInt_t idca   = 6;
	UInt_t id0pi  = 7;
	UInt_t id0K   = 8;
	UInt_t id0xd0     = 9;
	UInt_t ipointing  = 10;
	UInt_t iphi  = 11;
	UInt_t iz    = 12;
        UInt_t ipTD0pi = 13;
        UInt_t ipTD0K  = 14;

	const Double_t phimax = 2*TMath::Pi();

	//Setting up the container grid... 
	UInt_t nstep = 8; //number of selection steps
	const Int_t nvar   = 15 ; //number of variables on the grid:pt, y, cosThetaStar, pTpi, pTk, cT, dca, d0pi, d0K, d0xd0, cosPointingAngle, phi 
	const Int_t nbin0_0_4   = 8 ; //bins in pt from 0 to 4 GeV
	const Int_t nbin0_4_8   = 4 ; //bins in pt from 4 to 8 GeV
	const Int_t nbin0_8_10  = 2 ; //bins in pt from 8 to 10 GeV
	const Int_t nbin1  = 30 ; //bins in y
	const Int_t nbin2  = 30 ; //bins in cosThetaStar 
        // soft pion and D0 from D*
	const Int_t nbin3_0_1  = 8 ; //bins in ptPi from 0 to 4 GeV
	const Int_t nbin3_1_2  = 1 ; //bins in ptPi from 4 to 8 GeV
	const Int_t nbin3_2_10 = 1 ; //bins in ptPi from 8 to 10 GeV
	const Int_t nbin4_0_4  = 8 ; //bins in ptD0 from 0 to 4 GeV
	const Int_t nbin4_4_8  = 3 ; //bins in ptD0 from 4 to 8 GeV
	const Int_t nbin4_8_10 = 1 ; //bins in ptD0 from 8 to 10 GeV
        // D0 prongs - cutting variables
	const Int_t nbin5  = 20 ; //bins in cT
	const Int_t nbin6  = 20 ; //bins in dca
	const Int_t nbin7  = 100 ; //bins in d0pi
	const Int_t nbin8  = 100 ; //bins in d0K
	const Int_t nbin9  = 80 ; //bins in d0xd0
	const Int_t nbin10  = 100 ; //bins in cosPointingAngle
	const Int_t nbin11  = 15 ; //bins in Phi
       	const Int_t nbin12  = 60 ; //bins in z vertex
        // D0 prongs pt and phi
	const Int_t nbin5_0_4   = 8 ; //bins in ptPi from 0 to 4 GeV
	const Int_t nbin5_4_8   = 4 ; //bins in ptPi from 4 to 8 GeV
	const Int_t nbin5_8_10  = 8 ; //bins in ptPi from 8 to 10 GeV
	const Int_t nbin6_0_4   = 8 ; //bins in ptk from 0 to 4 GeV
	const Int_t nbin6_4_8   = 4 ; //bins in ptk from 4 to 8 GeV
	const Int_t nbin6_8_10  = 8 ; //bins in ptk from 8 to 10 GeV
        
	//arrays for the number of bins in each dimension
	Int_t iBin[nvar];

 	iBin[0]=nbin0_0_4+nbin0_4_8+nbin0_8_10;
	iBin[1]=nbin1;
	iBin[2]=nbin2;
 	iBin[3]=nbin3_0_1+nbin3_1_2+nbin3_2_10;
 	iBin[4]=nbin4_0_4+nbin4_4_8+nbin4_8_10;
	iBin[5]=nbin5;
	iBin[6]=nbin6;
	iBin[7]=nbin7;
	iBin[8]=nbin8;
	iBin[9]=nbin9;
	iBin[10]=nbin10;
	iBin[11]=nbin11;
	iBin[12]=nbin12;
        iBin[13]=nbin5_0_4+nbin5_4_8+nbin5_8_10;
	iBin[14]=nbin6_0_4+nbin6_4_8+nbin6_8_10;
	
	//arrays for lower bounds :
	Double_t *binLim0 = new Double_t[iBin[0]+1];
	Double_t *binLim1 = new Double_t[iBin[1]+1];
	Double_t *binLim2 = new Double_t[iBin[2]+1];
	Double_t *binLim3 = new Double_t[iBin[3]+1];
	Double_t *binLim4 = new Double_t[iBin[4]+1];
	Double_t *binLim5 = new Double_t[iBin[5]+1];
	Double_t *binLim6 = new Double_t[iBin[6]+1];
	Double_t *binLim7 = new Double_t[iBin[7]+1];
	Double_t *binLim8 = new Double_t[iBin[8]+1];
	Double_t *binLim9 = new Double_t[iBin[9]+1];
	Double_t *binLim10 = new Double_t[iBin[10]+1];
	Double_t *binLim11 = new Double_t[iBin[11]+1];
	Double_t *binLim12 = new Double_t[iBin[12]+1];
	Double_t *binLim13 = new Double_t[iBin[13]+1];
	Double_t *binLim14 = new Double_t[iBin[14]+1];

	// checking limits
	if (ptmax_0_4 != ptmin_4_8) {
		Error("AliCFTaskForDStarAnalysis","max lim 1st range != min lim 2nd range, please check!");
	}
	if (ptmax_4_8 != ptmin_8_10) {
		Error("AliCFTaskForDStarAnalysis","max lim 2nd range != min lim 3rd range, please check!");
	}

	// values for bin lower bounds
	// pt -----------------------------------------------------------------------------------------
	for(Int_t i=0; i<=nbin0_0_4; i++) binLim0[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbin0_0_4*(Double_t)i ; 
	if (binLim0[nbin0_0_4] != ptmin_4_8)  {
		Error("AliCFDStar","Calculated bin lim for pt - 1st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_4_8; i++) binLim0[i+nbin0_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbin0_4_8*(Double_t)i ; 
	if (binLim0[nbin0_0_4+nbin0_4_8] != ptmin_8_10)  {
		Error("AliCFDStar","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_8_10; i++) binLim0[i+nbin0_0_4+nbin0_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbin0_8_10*(Double_t)i ; 

	// y -----------------------------------------------------------------------------------------
	for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)ymin  + (ymax-ymin)  /nbin1*(Double_t)i ;

	// cosThetaStar -----------------------------------------------------------------------------
	for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)cosmin  + (cosmax-cosmin)  /nbin2*(Double_t)i ;

	// Soft ptPi ---------------------------------------------------------------------------------
	for(Int_t i=0; i<=nbin3_0_1; i++) binLim3[i]=(Double_t)ptmin_0_1 + (ptmax_0_1-ptmin_0_1)/nbin3_0_1*(Double_t)i ; 
	if (binLim3[nbin3_0_1] != ptmin_1_2)  {
		Error("AliCFDStar","Calculated bin lim for ptPi - 1st range - differs from expected!");
	}
	for(Int_t i=0; i<=nbin3_1_2; i++) binLim3[i+nbin3_0_1]=(Double_t)ptmin_1_2 + (ptmax_1_2-ptmin_1_2)/nbin3_1_2*(Double_t)i ; 
	if (binLim3[nbin3_0_1+nbin3_1_2] != ptmin_2_10)  {
		Error("AliCFDStar","Calculated bin lim for ptPi - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin3_2_10; i++) binLim3[i+nbin3_0_1+nbin3_1_2]=(Double_t)ptmin_2_10 + (ptmax_2_10-ptmin_2_10)/nbin3_2_10*(Double_t)i ; 

	// ptD0 --------------------------------------------------------------------------------------
	for(Int_t i=0; i<=nbin4_0_4; i++) binLim4[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbin4_0_4*(Double_t)i ; 
	if (binLim4[nbin4_0_4] != ptmin_4_8)  {
		Error("AliCFDStar","Calculated bin lim for ptKa - 1st range - differs from expected!");
	}
	for(Int_t i=0; i<=nbin4_4_8; i++) binLim4[i+nbin4_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbin4_4_8*(Double_t)i ; 
	if (binLim4[nbin4_0_4+nbin4_4_8] != ptmin_8_10)  {
		Error("AliCFDStar","Calculated bin lim for ptKa - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin4_8_10; i++) binLim4[i+nbin4_0_4+nbin4_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbin4_8_10*(Double_t)i ; 
 
        // D0 ptPi --------------------------------------------------------------------------------------------------------
        for(Int_t i=0; i<=nbin5_0_4; i++) binLim13[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbin5_0_4*(Double_t)i ; 
        if (binLim13[nbin5_0_4] != ptmin_4_8)  {
                Error("AliCFDStar","Calculated bin lim for ptPi - 1st range - differs from expected!");
        }
        for(Int_t i=0; i<=nbin5_4_8; i++) binLim13[i+nbin5_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbin5_4_8*(Double_t)i ; 
        if (binLim13[nbin5_0_4+nbin5_4_8] != ptmin_8_10)  {
                Error("AliCFDStar","Calculated bin lim for ptPi - 2nd range - differs from expected!\n");
        }
        for(Int_t i=0; i<=nbin5_8_10; i++) binLim13[i+nbin5_0_4+nbin5_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbin5_8_10*(Double_t)i ; 

        // D0 ptK ----------------------------------------------------------------------------------------------------------
        for(Int_t i=0; i<=nbin6_0_4; i++) binLim14[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbin6_0_4*(Double_t)i ; 
        if (binLim14[nbin6_0_4] != ptmin_4_8)  {
                Error("AliCFDStar","Calculated bin lim for ptKa - 1st range - differs from expected!");
        }
        for(Int_t i=0; i<=nbin6_4_8; i++) binLim14[i+nbin6_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbin6_4_8*(Double_t)i ; 
        if (binLim14[nbin6_0_4+nbin6_4_8] != ptmin_8_10)  {
                Error("AliCFDStar","Calculated bin lim for ptKa - 2nd range - differs from expected!\n");
        }
        for(Int_t i=0; i<=nbin6_8_10; i++) binLim14[i+nbin6_0_4+nbin6_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbin6_8_10*(Double_t)i ; 

	// cT ---------------------------------------------------------------------------------------------------------------
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

	// Phi
	for(Int_t i=0; i<=nbin11; i++) binLim11[i]=(Double_t)phimin  + (phimax-phimin)  /nbin11*(Double_t)i ;

	// z Primary Vertex
	for(Int_t i=0; i<=nbin12; i++) {
		binLim12[i]=(Double_t)zmin  + (zmax-zmin)  /nbin12*(Double_t)i ;
	}

	//one "container" for MC
	AliCFContainer* container = new AliCFContainer("container","container for tracks",nstep,nvar,iBin);
	//setting the bin limits
	container -> SetBinLimits(ipt,binLim0);
	container -> SetBinLimits(iy,binLim1);
	container -> SetBinLimits(icosThetaStar,binLim2);
	container -> SetBinLimits(ipTpi,binLim3);
	container -> SetBinLimits(ipTD0,binLim4);
	container -> SetBinLimits(icT,binLim5);
	container -> SetBinLimits(idca,binLim6);
	container -> SetBinLimits(id0pi,binLim7);
	container -> SetBinLimits(id0K,binLim8);
	container -> SetBinLimits(id0xd0,binLim9);
	container -> SetBinLimits(ipointing,binLim10);
	container -> SetBinLimits(iphi,binLim11);
	container -> SetBinLimits(iz,binLim12);
	container -> SetBinLimits(ipTD0pi,binLim13);
	container -> SetBinLimits(ipTD0K,binLim14);
	
	//CREATE THE  CUTS -----------------------------------------------
	
	// Gen-Level kinematic cuts
	AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
	
	//Particle-Level cuts:  
	AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
	mcGenCuts->SetRequirePdgCode(413, kTRUE);  // kTRUE set in order to include D*_bar
	mcGenCuts->SetAODMC(1); //special flag for reading MC in AOD tree (important)
	
	// Acceptance cuts:
	AliCFAcceptanceCuts* accCuts = new AliCFAcceptanceCuts("accCuts", "Acceptance cuts");
	AliCFTrackKineCuts *kineAccCuts = new AliCFTrackKineCuts("kineAccCuts","Kine-Acceptance cuts");
	kineAccCuts->SetPtRange(ptmin,ptmax);
	kineAccCuts->SetEtaRange(etamin,etamax);

	// Rec-Level kinematic cuts
	AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts","rec-level kine cuts");
	
	AliCFTrackQualityCuts *recQualityCuts = new AliCFTrackQualityCuts("recQualityCuts","rec-level quality cuts");
	
	AliCFTrackIsPrimaryCuts *recIsPrimaryCuts = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts","rec-level isPrimary cuts");
	
	printf("CREATE MC KINE CUTS\n");
	TObjArray* mcList = new TObjArray(0) ;
	mcList->AddLast(mcKineCuts);
	mcList->AddLast(mcGenCuts);
	
	printf("CREATE ACCEPTANCE CUTS\n");
	TObjArray* accList = new TObjArray(0) ;
	accList->AddLast(kineAccCuts);

	printf("CREATE RECONSTRUCTION CUTS\n");
	TObjArray* recList = new TObjArray(0) ;   // not used!! 
	recList->AddLast(recKineCuts);
	recList->AddLast(recQualityCuts);
	recList->AddLast(recIsPrimaryCuts);
	
	TObjArray* emptyList = new TObjArray(0);

	//CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
	printf("CREATE INTERFACE AND CUTS\n");
	AliCFManager* man = new AliCFManager() ;

	man->SetParticleContainer     (container);
	man->SetParticleCutsList(0 , mcList); // MC
	man->SetParticleCutsList(1 , accList); // Acceptance 
	man->SetParticleCutsList(2 , emptyList); // Vertex 
	man->SetParticleCutsList(3 , emptyList); // Refit 
	man->SetParticleCutsList(4 , emptyList); // AOD
	man->SetParticleCutsList(5 , emptyList); // AOD in Acceptance
	man->SetParticleCutsList(6 , emptyList); // AOD with required n. of ITS clusters
	man->SetParticleCutsList(7 , emptyList); // AOD Reco cuts
	
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
	AliCFTaskForDStarAnalysis *task = new AliCFTaskForDStarAnalysis("AliCFTaskForDStarAnalysis");
	task->SetMinITSClusters(minITSClusters);
	task->SetMinITSClustersSoft(minITSClustersSoft);
	task->SetCFManager(man); //here is set the CF manager
	
        Bool_t AcceptanceUnf = kTRUE; // unfold at acceptance level, otherwise D* cuts
        Int_t thnDim[4];
        
        //first half  : reconstructed 
        //second half : MC
        thnDim[0] = iBin[0];
        thnDim[2] = iBin[0];
        thnDim[1] = iBin[1];
        thnDim[3] = iBin[1];

        THnSparseD* correlation = new THnSparseD("correlation","THnSparse with correlations",4,thnDim);
        Double_t** binEdges = new Double_t[2];

        // set bin limits

        binEdges[0]= binLim0;
        binEdges[1]= binLim1;

        correlation->SetBinEdges(0,binEdges[0]);
        correlation->SetBinEdges(2,binEdges[0]);

        correlation->SetBinEdges(1,binEdges[1]);
        correlation->SetBinEdges(3,binEdges[1]);

        correlation->Sumw2();
  
        // correlation matrix ready
        //------------------------------------------------//

        task->SetCorrelationMatrix(correlation); // correlation matrix for unfolding
	
	// Create and connect containers for input/output
	
	// ------ input data ------
	AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
	
	// ----- output data -----
	
	TString outputfile = AliAnalysisManager::GetCommonFileName();
	outputfile += ":PWG3_D2H_CFtaskDStar";

	//now comes user's output objects :
	// output TH1I for event counting
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("CFDSchist0", TH1I::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	// output Correction Framework Container (for acceptance & efficiency calculations)
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("CFDSccontainer0", AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	// Unfolding - correlation matrix
        AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("CFDScorr0", THnSparseD::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

	mgr->AddTask(task);
	
	mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task,1,coutput1);
	mgr->ConnectOutput(task,2,coutput2);
        mgr->ConnectOutput(task,3,coutput3);

	return task;
}

