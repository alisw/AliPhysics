//DEFINITION OF A FEW CONSTANTS
const Double_t ymin  = -2.1 ;
const Double_t ymax  =  2.1 ;
// const Double_t ptmin_0_4 =  0.0 ;
// const Double_t ptmax_0_4 =  4.0 ;
// const Double_t ptmin_4_8 =  4.0 ;
// const Double_t ptmax_4_8 =  8.0 ;
// const Double_t ptmin_8_10 =  8.0 ;
// const Double_t ptmax_8_10 =  10.0 ;
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
//const Double_t phimax = 2Pi;  // defined in the macro!!!!!!!!!!!!!!  
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
const Int_t    minITSClusters = 5;

const Float_t centmin = 0.;
const Float_t centmax = 100.;
const Float_t fakemin = -0.5;
const Float_t fakemax = 2.5.;


//----------------------------------------------------

AliCFTaskVertexingHF *AddTaskCFVertexingHF3Prong(const char* cutFile = "./DplustoKpipiCuts.root",Bool_t isKeepDfromB=kFALSE, Bool_t isKeepDfromBOnly=kFALSE, Int_t pdgCode = 411, Char_t isSign = 2)
{
	printf("Addig CF task using cuts from file %s\n",cutFile);

	// isSign = 0 --> D0 only
	// isSign = 1 --> D0bar only
	// isSign = 2 --> D0 + D0bar
	
	TString expected;
	if (isSign == 0 && pdgCode < 0){
		AliError(Form("Error setting PDG code (%d) and sign (0 --> particle (%d) only): they are not compatible, returning",pdgCode));
		return 0x0;
	}
	else if (isSign == 1 && pdgCode > 0){
		AliError(Form("Error setting PDG code (%d) and sign (1 --> antiparticle (%d) only): they are not compatible, returning",pdgCode));
		return 0x0;
	}
	else if (isSign > 2 || isSign < 0){
		AliError(Form("Sign not valid (%d, possible values are 0, 1, 2), returning"));
		return 0x0;
	}

	TFile* fileCuts = new TFile(cutFile);
	AliRDHFCutsDplustoKpipi *cutsDplustoKpipi = (AliRDHFCutsDplustoKpipi*)fileCuts->Get("DplustoKpipiCuts");
	
	// check that the fKeepD0fromB flag is set to true when the fKeepD0fromBOnly flag is true
	//  for now the binning is the same than for all D's
	if(isKeepDfromBOnly) isKeepDfromB = true;
	

	/*
	  Double_t ptmin_0_4;
	  Double_t ptmax_0_4;
	  Double_t ptmin_4_8;
	  Double_t ptmax_4_8;
	  Double_t ptmin_8_10;
	  Double_t ptmax_8_10;
	  
	  if(!isKeepDfromB){
	  ptmin_0_4 =  0.0 ;
	  ptmax_0_4 =  4.0 ;
	  ptmin_4_8 =  4.0 ;
	  ptmax_4_8 =  8.0 ;
	  ptmin_8_10 =  8.0 ;
	  ptmax_8_10 =  10.0 ;
	  } else{
	  ptmin_0_4 =  0.0 ;
	  ptmax_0_4 =  3.0 ;
	  ptmin_4_8 =  3.0 ;
	  ptmax_4_8 =  5.0 ;
	  ptmin_8_10 =  5.0 ;
	  ptmax_8_10 =  10.0 ;
	  }
	*/

	//CONTAINER DEFINITION
	Info("AliCFTaskVertexingHF","SETUP CONTAINER");
	//the sensitive variables, their indices
	UInt_t ipt = 0;
	UInt_t iy  = 1;
	UInt_t iphi  = 2;
	UInt_t icT  = 3;
	UInt_t ipointing  = 4;
	UInt_t iptpi  = 5;
	UInt_t iptK  = 6;
	UInt_t iptpi2  = 7;
	UInt_t id0pi  = 8;
	UInt_t id0K  = 9;
	UInt_t id0pi2  = 10;
	UInt_t iz  = 11;
	UInt_t icent = 12;
	UInt_t ifake = 13;

	const Double_t phimax = 2*TMath::Pi();

	//Setting up the container grid... 
	UInt_t nstep = 10; //number of selection steps: MC with limited acceptance, MC, Acceptance, Vertex, Refit, Reco (no cuts), RecoAcceptance, RecoITSClusters (RecoAcceptance included), RecoPPR (RecoAcceptance+RecoITSCluster included), RecoPID 
	const Int_t nvar   = 13 ; //number of variables on the grid:pt, y, cosThetaStar, pTpi, pTk, cT, dca, d0pi, d0K, d0xd0, cosPointingAngle, phi 
// 	const Int_t nbin0_0_4  = 8 ; //bins in pt from 0 to 4 GeV
// 	const Int_t nbin0_4_8  = 4 ; //bins in pt from 4 to 8 GeV
// 	const Int_t nbin0_8_10  = 1 ; //bins in pt from 8 to 10 GeV

/*
	Int_t nbin0_0_4;
	Int_t nbin0_4_8;
	Int_t nbin0_8_10;
	if (!isKeepDfromB){
	  nbin0_0_4  = 8 ; //bins in pt from 0 to 4 GeV
	  nbin0_4_8  = 4 ; //bins in pt from 4 to 8 GeV
	  nbin0_8_10  = 1 ; //bins in pt from 8 to 10 GeV
	}else{
	  nbin0_0_4  = 3 ; //bins in pt from 0 to 3 GeV
	  nbin0_4_8  = 1 ; //bins in pt from 3 to 5 GeV
	  nbin0_8_10  = 1 ; //bins in pt from 5 to 10 GeV
	}
*/
	const Int_t nbin0 = cutsDplustoKpipi->GetNPtBins(); // bins in pT
	printf("pT: nbin (from cuts file) = %d\n",nbin0);
	const Int_t nbin1  = 42 ; //bins in y
	const Int_t nbin2  = 20 ; //bins in phi
	const Int_t nbin3  = 24 ; //bins in cT 
	const Int_t nbin4  = 1050 ; //bins in cosPointingAngle	
	const Int_t nbin5_0_4  = 8 ; //bins in ptPi from 0 to 4 GeV
	const Int_t nbin5_4_8  = 4 ; //bins in ptPi from 4 to 8 GeV
	const Int_t nbin5_8_10  = 1 ; //bins in ptPi from 8 to 10 GeV
	const Int_t nbin6_0_4  = 8 ; //bins in ptKa from 0 to 4 GeV
	const Int_t nbin6_4_8  = 4 ; //bins in ptKa from 4 to 8 GeV
	const Int_t nbin6_8_10  = 1 ; //bins in ptKa from 8 to 10 GeV
	const Int_t nbin7_0_4  = 8 ; //bins in ptpi2 from 0 to 4 GeV
	const Int_t nbin7_4_8  = 4 ; //bins in ptpi2 from 4 to 8 GeV
	const Int_t nbin7_8_10  = 1 ; //bins in ptpi2 from 8 to 10 GeV
	const Int_t nbin8  = 100 ; //bins in d0pi
	const Int_t nbin9  = 100 ; //bins in d0K
	const Int_t nbin10  = 100 ; //bins in d0pi2
	const Int_t nbin11  = 60 ; //bins in z vertex
	const Int_t nbin12 = 10; //bins in centrality
	const Int_t nbin13 = 3;  //bins in fake
	
	//arrays for the number of bins in each dimension
	Int_t iBin[nvar];
 	//iBin[0]=nbin0_0_4+nbin0_4_8+nbin0_8_10;
	iBin[0]=nbin0;
	iBin[1]=nbin1;
	iBin[2]=nbin2;
	// 	iBin[3]=nbin3_0_4+nbin3_4_8+nbin3_8_10;
 	//iBin[4]=nbin4_0_4+nbin4_4_8+nbin4_8_10;
 	iBin[3]=nbin3;
 	iBin[4]=nbin4;
	iBin[5]=nbin0;
	iBin[6]=nbin0;
	iBin[7]=nbin0;
	iBin[8]=nbin8;
	iBin[9]=nbin9;
	iBin[10]=nbin10;
	iBin[11]=nbin11;
	iBin[12]=nbin12;
	
	
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
	Double_t *binLim11=new Double_t[iBin[11]+1];
	Double_t *binLim12=new Double_t[iBin[12]+1];

	// checking limits
	/*
	if (ptmax_0_4 != ptmin_4_8) {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","max lim 1st range != min lim 2nd range, please check!");
	}
	if (ptmax_4_8 != ptmin_8_10) {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","max lim 2nd range != min lim 3rd range, please check!");
	}
	*/
	// values for bin lower bounds
	// pt
	Float_t* floatbinLim0 = cutsDplustoKpipi->GetPtBinLimits();
	for (Int_t ibin0 = 0 ; ibin0<iBin[0]+1; ibin0++){
		binLim0[ibin0] = (Double_t)floatbinLim0[ibin0];
		binLim5[ibin0] = (Double_t)floatbinLim0[ibin0];
		binLim6[ibin0] = (Double_t)floatbinLim0[ibin0];
		binLim7[ibin0] = (Double_t)floatbinLim0[ibin0];
	}
	for(Int_t i=0; i<=nbin0; i++) printf("binLim0[%d]=%f\n",i,binLim0[i]);  

	/*
	for(Int_t i=0; i<=nbin0_0_4; i++) binLim0[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbin0_0_4*(Double_t)i ; 
	if (binLim0[nbin0_0_4] != ptmin_4_8)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 1st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_4_8; i++) binLim0[i+nbin0_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbin0_4_8*(Double_t)i ; 
	if (binLim0[nbin0_0_4+nbin0_4_8] != ptmin_8_10)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbin0_8_10; i++) binLim0[i+nbin0_0_4+nbin0_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbin0_8_10*(Double_t)i ; 
	*/

	// y
	for(Int_t i=0; i<=nbin1; i++) binLim1[i]=(Double_t)ymin  + (ymax-ymin)  /nbin1*(Double_t)i ;

	// cosThetaStar
	//	for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)cosmin  + (cosmax-cosmin)  /nbin2*(Double_t)i ;
	// Phi
	for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)phimin  + (phimax-phimin)  /nbin2*(Double_t)i ;

	// cT
	for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)cTmin  + (cTmax-cTmin)  /nbin3*(Double_t)i ;

	// cosPointingAngle
	for(Int_t i=0; i<=nbin4; i++) binLim4[i]=(Double_t)cosmin  + (cosmax-cosmin)  /nbin4*(Double_t)i ;

	/*
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
	*/


	// dca
	//for(Int_t i=0; i<=nbin6; i++) binLim6[i]=(Double_t)dcamin  + (dcamax-dcamin)  /nbin6*(Double_t)i ;

	// d0pi
	for(Int_t i=0; i<=nbin8; i++) binLim8[i]=(Double_t)d0min  + (d0max-d0min)  /nbin8*(Double_t)i ;

	// d0K
	for(Int_t i=0; i<=nbin9; i++) binLim9[i]=(Double_t)d0min  + (d0max-d0min)  /nbin9*(Double_t)i ;

	// d0pi2
	for(Int_t i=0; i<=nbin10; i++) binLim10[i]=(Double_t)d0min  + (d0max-d0min)  /nbin10*(Double_t)i ;





	// z Primary Vertex
	for(Int_t i=0; i<=nbin11; i++) {
		binLim11[i]=(Double_t)zmin  + (zmax-zmin)  /nbin11*(Double_t)i ;
		//		Info("AliCFHeavyFlavourTaskMultiVarMultiStep",Form("i-th bin, lower limit = %f", binLim12[i]));
	}

	// centrality
	for(Int_t i=0; i<=nbin12; i++) {
	  binLim12[i]=(Double_t)centmin  + (centmax-centmin)/nbin12 * (Double_t)i;
	}

	// fake
	for(Int_t i=0; i<=nbin13; i++) {
	  binLim13[i]=(Double_t)fakemin  + (fakemax-fakemin)/nbin13 * (Double_t)i;
	}

	

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
	TString nameContainer="";
	if(!isKeepDfromB) {
	  nameContainer="CFHFccontainer0_3Prong_CommonFramework";
	}
	else  if(isKeepDfromBOnly){
	  nameContainer="CFHFccontainer0DfromB_3Prong_CommonFramework";
	}
	else  {
	  nameContainer="CFHFccontainer0allD_3Prong_CommonFramework";          
	}
	
	AliCFContainer* container = new AliCFContainer(nameContainer,"container for tracks",nstep,nvar,iBin);
	//setting the bin limits
	printf("pt\n");
	container -> SetBinLimits(ipt,binLim0);
	printf("y\n");
	container -> SetBinLimits(iy,binLim1);
	printf("Phi\n");
	container -> SetBinLimits(iphi,binLim2);
	printf("cT\n");
	container -> SetBinLimits(icT,binLim3);
	printf("pointing angle\n");
	container -> SetBinLimits(ipointing,binLim4);
	printf("ptpi\n");
	container -> SetBinLimits(iptpi,binLim5);
	printf("ptK\n");
	container -> SetBinLimits(iptK,binLim6);
	printf("ptpi2\n");
	container -> SetBinLimits(iptpi2,binLim7);
	printf("d0pi\n");
	container -> SetBinLimits(id0pi,binLim8);
	printf("d0K\n");
	container -> SetBinLimits(id0K,binLim9);
	printf("d0pi2\n");
	container -> SetBinLimits(id0pi2,binLim10);
	printf("z \n");
	container -> SetBinLimits(iz,binLim11);
	printf("cent\n");
	container -> SetBinLimits(icent,binLim12);
	printf("fake\n");
	container -> SetBinLimits(ifake,binLim13);
	
	
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
        container -> SetVarTitle(iphi, "phi");
        container -> SetVarTitle(icT, "ct");
	container -> SetVarTitle(ipointing, "pionting");        
	container -> SetVarTitle(iptpi, "ptpi");
	container -> SetVarTitle(iptK, "ptK");
	container -> SetVarTitle(iptpi2, "ptpi2");
	container -> SetVarTitle(id0pi, "d0pi");
        container -> SetVarTitle(id0K, "d0K");
	container -> SetVarTitle(id0pi2, "d0pi2");
	container -> SetVarTitle(iz, "z");
	container -> SetVarTitle(icent, "centrality");
	container -> SetVarTitle(ifake, "fake");


	//CREATE THE  CUTS -----------------------------------------------
	
	// Gen-Level kinematic cuts
	AliCFTrackKineCuts *mcKineCuts = new AliCFTrackKineCuts("mcKineCuts","MC-level kinematic cuts");
	
	//Particle-Level cuts:  
	AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts("mcGenCuts","MC particle generation cuts");
	Bool_t useAbsolute = kTRUE;
	if (isSign != 2){
		useAbsolute = kFALSE;
	}
	mcGenCuts->SetRequirePdgCode(pdgCode, useAbsolute);  // kTRUE set in order to include antiparticle
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
	AliCFTaskVertexingHF *task = new AliCFTaskVertexingHF("AliCFTaskVertexingHF",cutsDplustoKpipi);
	task->SetFillFromGenerated(kFALSE);
	task->SetDecayChannel(31);
	task->SetUseWeight(kFALSE);
	task->SetCFManager(man); //here is set the CF manager
	task->SetSign(isSign);
	task->SetCentralitySelection(kTRUE);

	if (isKeepDfromB && !isKeepDfromBOnly) task->SetDselection(2);
	if (isKeepDfromB && isKeepDfromBOnly) task->SetDselection(1);		

	Printf("***************** CONTAINER SETTINGS *****************");
	Printf("decay channel = %d",(Int_t)task->GetDecayChannel());
	Printf("FillFromGenerated = %d",(Int_t)task->GetFillFromGenerated());
	Printf("Dselection = %d",(Int_t)task->GetDselection());
	Printf("UseWeight = %d",(Int_t)task->GetUseWeight());
	Printf("Sign = %d",(Int_t)task->GetSign());
	Printf("Centrality selection = %d",(Int_t)task->GetCentralitySelection());
	Printf("Fake selection = %d",(Int_t)task->GetFakeSelection());
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
		nameCorr="CFHFcorr0_3Prong_CommonFramework";
	}
	else  if(isKeepDfromBOnly){
		nameCorr= "CFHFcorr0KeepDfromBOnly_3Prong_CommonFramework";
	}
	else  {
		nameCorr="CFHFcorr0allD_3Prong_CommonFramework";		
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
	TString output1name="", output2name="", output3name="", output4name="";;
	output2name=nameContainer;
	output3name=nameCorr;
	if(!isKeepDfromB) {
		outputfile += ":PWG3_D2H_CFtaskDplustoKpipi_CommonFramework";
		output1name="CFHFchist0_3Prong_CommonFramework";
	}
	else  if(isKeepDfromBOnly){
		outputfile += ":PWG3_D2H_CFtaskDplustoKpipiKeepDfromBOnly_CommonFramework";
		output1name="CFHFchist0DfromB_3Prong_CommonFramework";
	}
	else{
		outputfile += ":PWG3_D2H_CFtaskDplustoKpipiKeepDfromB_CommonFramework";
		output1name="CFHFchist0allD_3Prong_CommonFramework";
	}

	output4name= "Cuts_3Prong_CommonFramework";

	//now comes user's output objects :
	// output TH1I for event counting
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(output1name, TH1I::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	// output Correction Framework Container (for acceptance & efficiency calculations)
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(output2name, AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	// Unfolding - correlation matrix
        AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(output3name, THnSparseD::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(output4name, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());

	mgr->AddTask(task);
	
	mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task,1,coutput1);
	mgr->ConnectOutput(task,2,coutput2);
        mgr->ConnectOutput(task,3,coutput3);
	mgr->ConnectOutput(task,4,coutput4);

	return task;
}

