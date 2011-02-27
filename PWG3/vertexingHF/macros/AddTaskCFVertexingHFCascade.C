//DEFINITION OF A FEW CONSTANTS
const Double_t ymin  = -2.1 ;
const Double_t ymax  =  2.1 ;
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
const Double_t ptmin = 0.1;
const Double_t ptmax = 9999.;
const Double_t etamin = -0.9;
const Double_t etamax = 0.9;
const Double_t zmin = -15;
const Double_t zmax = 15;
const Int_t    minITSClusters = 1;

//----------------------------------------------------

AliCFTaskVertexingHF *AddTaskCFVertexingHFCascade(const char* cutFile = "./DStartoKpipiCuts.root",Bool_t isKeepDfromB=kFALSE, Bool_t isKeepDfromBOnly=kFALSE, Int_t pdgCode = 413, Char_t isSign = 2)
{
	printf("Addig CF task using cuts from file %s\n",cutFile);

	// isSign = 0 --> D0 only
	// isSign = 1 --> D0bar only
	// isSign = 2 --> D0 + D0bar

	TString expected;
	if (isSign == 0 && pdgCode < 0){
		AliError(Form("Error setting PDG code (%d) and sign (0 --> D0 only): they are not compatible, returning"));
		return 0x0;
	}
	else if (isSign == 1 && pdgCode > 0){
		AliError(Form("Error setting PDG code (%d) and sign (1 --> D0bar only): they are not compatible, returning"));
		return 0x0;
	}
	else if (isSign > 2 || isSign < 0){
		AliError(Form("Sign not valid (%d, possible values are 0, 1, 2), returning"));
		return 0x0;
	}

	TFile* fileCuts = new TFile(cutFile);
	AliRDHFCutsDStartoKpipi *cutsD0toKpi = (AliRDHFCutsDStartoKpipi*)fileCuts->Get("DStartoKpipiCuts");
	
	// check that the fKeepD0fromB flag is set to true when the fKeepD0fromBOnly flag is true
	//  for now the binning is the same than for all D's
	if(isKeepDfromBOnly) isKeepDfromB = true;
	
	Double_t ptmin_0_6;
	Double_t ptmax_0_6;
	Double_t ptmin_6_8;
	Double_t ptmax_6_8;
	Double_t ptmin_8_16;
	Double_t ptmax_8_16;
	Double_t ptmin_16_24;
	Double_t ptmax_16_24;
	
	ptmin_0_6 =  0.0 ;
	ptmax_0_6 =  6.0 ;
	ptmin_6_8 =  6.0 ;
	ptmax_6_8 =  8.0 ;
	ptmin_8_16 =  8.0 ;
	ptmax_8_16 =  16.0 ;
	ptmin_16_24 =  16.0 ;
	ptmax_16_24 =  24.0 ;

	//CONTAINER DEFINITION
	Info("AliCFTaskVertexingHF","SETUP CONTAINER");
	//the sensitive variables, their indices
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
	UInt_t iphi  = 11;
	UInt_t iz  = 12;

	const Double_t phimax = 2*TMath::Pi();

	//Setting up the container grid... 
	UInt_t nstep = 10; 
	const Int_t nvar   = 13;

	//arrays for the number of bins in each dimension
	Int_t iBin[nvar];

	//OPTION 1: defining the pt, ptPi, ptK bins by hand...		
	const Int_t nbin0_0_6  = 6 ; //bins in pt from 0 to 6 GeV
	const Int_t nbin0_6_8  = 1 ; //bins in pt from 6 to 8 GeV
	const Int_t nbin0_8_16  = 2 ; //bins in pt from 8 to 16 GeV
	const Int_t nbin0_16_24  = 1 ; //bins in pt from 16 to 24 GeV
	const Int_t nbin3_0_6  = 6 ; //bins in ptPi from 0 to 6 GeV
	const Int_t nbin3_6_8  = 1 ; //bins in ptPi from 6 to 8 GeV
	const Int_t nbin3_8_16  = 2 ; //bins in ptPi from 8 to 16 GeV
	const Int_t nbin3_16_24  = 1 ; //bins in ptPi from 16 to 24 GeV
	const Int_t nbin4_0_6  = 6 ; //bins in ptK from 0 to 6 GeV
	const Int_t nbin4_6_8  = 1 ; //bins in ptK from 6 to 8 GeV
	const Int_t nbin4_8_16  = 2 ; //bins in ptK from 8 to 16 GeV
	const Int_t nbin4_16_24  = 1 ; //bins in ptK from 16 to 24 GeV
 	iBin[0]=nbin0_0_6+nbin0_6_8+nbin0_8_16+nbin0_16_24;
      	iBin[3]=nbin3_0_6+nbin3_6_8+nbin3_8_16+nbin3_16_24;
 	iBin[4]=nbin4_0_6+nbin4_6_8+nbin4_8_16+nbin4_16_24;
	Double_t *binLim0=new Double_t[iBin[0]+1];
	Double_t *binLim3=new Double_t[iBin[3]+1];
	Double_t *binLim4=new Double_t[iBin[4]+1];

	// values for bin lower bounds
	// pt
	for(Int_t i=0; i<=nbin0_0_6; i++) binLim0[i]=(Double_t)ptmin_0_6 + (ptmax_0_6-ptmin_0_6)/nbin0_0_6*(Double_t)i ; 
	if (binLim0[nbin0_0_6] != ptmin_6_8)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 1st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_6_8; i++) binLim0[i+nbin0_0_6]=(Double_t)ptmin_6_8 + (ptmax_6_8-ptmin_6_8)/nbin0_6_8*(Double_t)i ; 
	if (binLim0[nbin0_0_6+nbin0_6_8] != ptmin_8_16)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_8_16; i++) binLim0[i+nbin0_0_6+nbin0_6_8]=(Double_t)ptmin_8_16 + (ptmax_8_16-ptmin_8_16)/nbin0_8_16*(Double_t)i ; 
	if (binLim0[nbin0_0_6+nbin0_6_8+nbin0_8_16] != ptmin_16_24)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_16_24; i++) binLim0[i+nbin0_0_6+nbin0_6_8+nbin0_8_16]=(Double_t)ptmin_16_24 + (ptmax_16_24-ptmin_16_24)/nbin0_16_24*(Double_t)i ; 

	// ptPi
	for(Int_t i=0; i<=nbin3_0_6; i++) binLim3[i]=(Double_t)ptmin_0_6 + (ptmax_0_6-ptmin_0_6)/nbin3_0_6*(Double_t)i ; 
	if (binLim3[nbin3_0_6] != ptmin_6_8)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 1st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin3_6_8; i++) binLim3[i+nbin3_0_6]=(Double_t)ptmin_6_8 + (ptmax_6_8-ptmin_6_8)/nbin3_6_8*(Double_t)i ; 
	if (binLim3[nbin3_0_6+nbin3_6_8] != ptmin_8_16)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin3_8_16; i++) binLim3[i+nbin3_0_6+nbin0_6_8]=(Double_t)ptmin_8_16 + (ptmax_8_16-ptmin_8_16)/nbin3_8_16*(Double_t)i ; 
	if (binLim3[nbin3_0_6+nbin3_6_8+nbin3_8_16] != ptmin_16_24)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin3_16_24; i++) binLim3[i+nbin3_0_6+nbin3_6_8+nbin3_8_16]=(Double_t)ptmin_16_24 + (ptmax_16_24-ptmin_16_24)/nbin3_16_24*(Double_t)i ; 

	// ptKa
	for(Int_t i=0; i<=nbin4_0_6; i++) binLim4[i]=(Double_t)ptmin_0_6 + (ptmax_0_6-ptmin_0_6)/nbin4_0_6*(Double_t)i ; 
	if (binLim4[nbin4_0_6] != ptmin_6_8)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 1st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin4_6_8; i++) binLim4[i+nbin4_0_6]=(Double_t)ptmin_6_8 + (ptmax_6_8-ptmin_6_8)/nbin4_6_8*(Double_t)i ; 
	if (binLim4[nbin4_0_6+nbin4_6_8] != ptmin_8_16)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin4_8_16; i++) binLim4[i+nbin4_0_6+nbin0_6_8]=(Double_t)ptmin_8_16 + (ptmax_8_16-ptmin_8_16)/nbin4_8_16*(Double_t)i ; 
	if (binLim4[nbin4_0_6+nbin4_6_8+nbin4_8_16] != ptmin_16_24)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin4_16_24; i++) binLim4[i+nbin4_0_6+nbin4_6_8+nbin4_8_16]=(Double_t)ptmin_16_24 + (ptmax_16_24-ptmin_16_24)/nbin4_16_24*(Double_t)i ; 
	
	//OPTION 2: ...or from the cuts file

	//const Int_t nbin0 = cutsD0toKpi->GetNPtBins(); // bins in pT
	//iBin[0]=nbin0;
 	//iBin[3]=nbin0;
 	//iBin[4]=nbin0;
	// values for bin lower bounds
	//Float_t* floatbinLim0 = cutsD0toKpi->GetPtBinLimits();
	//for (Int_t ibin0 = 0 ; ibin0<iBin[0]+1; ibin0++){
	//	binLim0[ibin0] = (Double_t)floatbinLim0[ibin0];
	//	binLim3[ibin0] = (Double_t)floatbinLim0[ibin0];
	//	binLim4[ibin0] = (Double_t)floatbinLim0[ibin0];
	//}
	//for(Int_t i=0; i<=nbin0; i++) printf("binLim0[%d]=%f\n",i,binLim0[i]);  

	//printf("pT: nbin (from cuts file) = %d\n",nbin0);

	// defining now the binning for the other variables:

	const Int_t nbin1  = 42 ; //bins in y
	const Int_t nbin2  = 42 ; //bins in cosThetaStar 
	const Int_t nbin5  = 24 ; //bins in cT
	const Int_t nbin6  = 24 ; //bins in dca
	const Int_t nbin7  = 100 ; //bins in d0pi
	const Int_t nbin8  = 100 ; //bins in d0K
	const Int_t nbin9  = 80 ; //bins in d0xd0
	const Int_t nbin10  = 1050 ; //bins in cosPointingAngle
	const Int_t nbin11  = 20 ; //bins in Phi
	const Int_t nbin12  = 60 ; //bins in z vertex

	iBin[1]=nbin1;
	iBin[2]=nbin2;
	iBin[5]=nbin5;
	iBin[6]=nbin6;
	iBin[7]=nbin7;
	iBin[8]=nbin8;
	iBin[9]=nbin9;
	iBin[10]=nbin10;
	iBin[11]=nbin11;
	iBin[12]=nbin12;
	
	//arrays for lower bounds :
	Double_t *binLim1=new Double_t[iBin[1]+1];
	Double_t *binLim2=new Double_t[iBin[2]+1];
	Double_t *binLim5=new Double_t[iBin[5]+1];
	Double_t *binLim6=new Double_t[iBin[6]+1];
	Double_t *binLim7=new Double_t[iBin[7]+1];
	Double_t *binLim8=new Double_t[iBin[8]+1];
	Double_t *binLim9=new Double_t[iBin[9]+1];
	Double_t *binLim10=new Double_t[iBin[10]+1];
	Double_t *binLim11=new Double_t[iBin[11]+1];
	Double_t *binLim12=new Double_t[iBin[12]+1];

	// y
	for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)ymin  + (ymax-ymin)  /nbin1*(Double_t)i ;

	// cosThetaStar
	for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)cosmin  + (cosmax-cosmin)  /nbin2*(Double_t)i ;
	
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

	// Phi
	for(Int_t i=0; i<=nbin11; i++) binLim11[i]=(Double_t)phimin  + (phimax-phimin)  /nbin11*(Double_t)i ;

	// z Primary Vertex
	for(Int_t i=0; i<=nbin12; i++) {
		binLim12[i]=(Double_t)zmin  + (zmax-zmin)  /nbin12*(Double_t)i ;
	}

	//one "container" for MC
	TString nameContainer="";
	if(!isKeepDfromB) {
		nameContainer="CFHFccontainer0_CommonFramework";
	}
	else  if(isKeepDfromBOnly){
		nameContainer="CFHFccontainer0DfromB_CommonFramework";
	}
	else  {
		nameContainer="CFHFccontainer0allD_CommonFramework";	  
	}

	AliCFContainer* container = new AliCFContainer(nameContainer,"container for tracks",nstep,nvar,iBin);
	//setting the bin limits
	printf("pt\n");
	container -> SetBinLimits(ipt,binLim0);
	printf("y\n");
	container -> SetBinLimits(iy,binLim1);
	printf("cts\n");
	container -> SetBinLimits(icosThetaStar,binLim2);
	printf("ptPi\n");
	container -> SetBinLimits(ipTpi,binLim3);
	printf("ptK\n");
	container -> SetBinLimits(ipTk,binLim4);
	printf("cT\n");
	container -> SetBinLimits(icT,binLim5);
	printf("dca\n");
	container -> SetBinLimits(idca,binLim6);
	printf("d0Pi\n");
	container -> SetBinLimits(id0pi,binLim7);
	printf("d0K\n");
	container -> SetBinLimits(id0K,binLim8);
	printf("d0xd0\n");
	container -> SetBinLimits(id0xd0,binLim9);
	printf("pointing\n");
	container -> SetBinLimits(ipointing,binLim10);
	printf("phi\n");
	container -> SetBinLimits(iphi,binLim11);
	printf("z\n");
	container -> SetBinLimits(iz,binLim12);
	
	container -> SetStepTitle(0, "MCLimAcc");
	container -> SetStepTitle(1, "MC");
        container -> SetStepTitle(2, "MCAcc");
        container -> SetStepTitle(3, "RecoVertex");
        container -> SetStepTitle(4, "RecoRefit");
        container -> SetStepTitle(5, "Reco");
        container -> SetStepTitle(6, "RecoAcc");
	container -> SetStepTitle(7, "RecoITSCluster");
	container -> SetStepTitle(8, "RecoCuts");
	container -> SetStepTitle(8, "RecoPID");

        container -> SetVarTitle(ipt,"pt");
	container -> SetVarTitle(iy,"y");
        container -> SetVarTitle(icosThetaStar, "cosThetaStar");
        container -> SetVarTitle(ipTpi, "ptpi");
	container -> SetVarTitle(ipTk, "ptK");
        container -> SetVarTitle(icT, "ct");
        container -> SetVarTitle(idca, "dca");
        container -> SetVarTitle(id0pi, "d0pi");
        container -> SetVarTitle(id0K, "d0K");
	container -> SetVarTitle(id0xd0, "d0xd0");
	container -> SetVarTitle(ipointing, "piointing");
	container -> SetVarTitle(iphi, "phi");
	container -> SetVarTitle(iz, "z");


	//CREATE THE  CUTS -----------------------------------------------
	
	// Gen-Level kinematic cuts
	AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
	
	//Particle-Level cuts:  
	AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
	Bool_t useAbsolute = kTRUE;
	if (isSign != 2){
		useAbsolute = kFALSE;
	}
	mcGenCuts->SetRequirePdgCode(pdgCode, useAbsolute);  // kTRUE set in order to include D0_bar
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
	man->SetParticleContainer(container);
	man->SetParticleCutsList(0 , mcList); // MC, Limited Acceptance
	man->SetParticleCutsList(1 , mcList); // MC
	man->SetParticleCutsList(2 , accList); // Acceptance 
	man->SetParticleCutsList(3 , emptyList); // Vertex 
	man->SetParticleCutsList(4 , emptyList); // Refit 
	man->SetParticleCutsList(5 , emptyList); // AOD
	man->SetParticleCutsList(6 , emptyList); // AOD in Acceptance
	man->SetParticleCutsList(7 , emptyList); // AOD with required n. of ITS clusters
	man->SetParticleCutsList(8 , emptyList); // AOD Reco (PPR cuts implemented in Task)
	man->SetParticleCutsList(9 , emptyList); // AOD Reco PID
	
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
	AliCFTaskVertexingHF *task = new AliCFTaskVertexingHF("AliCFTaskVertexingHF",cutsD0toKpi);
	task->SetFillFromGenerated(kFALSE);
	task->SetCFManager(man); //here is set the CF manager
	task->SetDecayChannel(21);
	task->SetUseWeight(kFALSE);
	task->SetSign(isSign);
	task->SetDebugLevel(10);
	if (isKeepDfromB && !isKeepDfromBOnly) task->SetDselection(2);
	if (isKeepDfromB && isKeepDfromBOnly) task->SetDselection(1);		

	Printf("***************** CONTAINER SETTINGS *****************");
	Printf("decay channel = %d",(Int_t)task->GetDecayChannel());
	Printf("FillFromGenerated = %d",(Int_t)task->GetFillFromGenerated());
	Printf("Dselection = %d",(Int_t)task->GetDselection());
	Printf("UseWeight = %d",(Int_t)task->GetUseWeight());
	Printf("Sign = %d",(Int_t)task->GetSign());
	Printf("***************END CONTAINER SETTINGS *****************\n");

        //-----------------------------------------------------------//
        //   create correlation matrix for unfolding - only eta-pt   //
        //-----------------------------------------------------------//

        Bool_t AcceptanceUnf = kTRUE; // unfold at acceptance level, otherwise PPR

        Int_t thnDim[4];
        
        //first half  : reconstructed 
        //second half : MC

        thnDim[0] = iBin[0];
        thnDim[2] = iBin[0];
        thnDim[1] = iBin[1];
        thnDim[3] = iBin[1];

	TString nameCorr="";
	if(!isKeepDfromB) {
		nameCorr="CFHFcorr0_CommonFramework";
	}
	else  if(isKeepDfromBOnly){
		nameCorr= "CFHFcorr0KeepDfromBOnly_CommonFramework";
	}
	else  {
		nameCorr="CFHFcorr0allD_CommonFramework";

	}

        THnSparseD* correlation = new THnSparseD(nameCorr,"THnSparse with correlations",4,thnDim);
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
	TString output1name="", output2name="", output3name="",output4name="";
	output2name=nameContainer;
	output3name=nameCorr;
	if(!isKeepDfromB) {
		outputfile += ":PWG3_D2H_CFtaskD0toKpi_CommonFramework";
		output1name="CFHFchist0_CommonFramework";
	}
	else  if(isKeepDfromBOnly){
		outputfile += ":PWG3_D2H_CFtaskD0toKpiKeepDfromBOnly_CommonFramework";
		output1name="CFHFchist0DfromB_CommonFramework";
	}
	else{
		outputfile += ":PWG3_D2H_CFtaskD0toKpiKeepDfromB_CommonFramework";
		output1name="CFHFchist0allD_CommonFramework";
	}
	output4name= "Cuts_CommonFramework";

	//now comes user's output objects :
	// output TH1I for event counting
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(output1name, TH1I::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	// output Correction Framework Container (for acceptance & efficiency calculations)
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(output2name, AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	// Unfolding - correlation matrix
        AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(output3name, THnSparseD::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	// cuts
	AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(output4name, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());

	mgr->AddTask(task);
	
	mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task,1,coutput1);
	mgr->ConnectOutput(task,2,coutput2);
        mgr->ConnectOutput(task,3,coutput3);
	mgr->ConnectOutput(task,4,coutput4);
	return task;
}

