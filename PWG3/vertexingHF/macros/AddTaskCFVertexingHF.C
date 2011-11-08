//DEFINITION OF A FEW CONSTANTS
const Double_t ymin  = -1.2 ;
const Double_t ymax  =  1.2 ;
const Double_t cosminTS = -1.05;
const Double_t cosmaxTS =  1.05;
const Double_t cosmin = 0.7;
const Double_t cosmax =  1.02;
const Double_t cTmin = 0;  // micron
const Double_t cTmax = 300;  // micron
const Double_t dcamin = 0;  // micron
const Double_t dcamax = 600;  // micron
const Double_t d0xd0min = -80000;  // micron
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
const Int_t    minITSClusters = 5;

const Float_t centmin_0_10 = 0.;
const Float_t centmax_0_10 = 10.;
const Float_t centmin_10_100 = 10.;
const Float_t centmax_10_100 = 100.;
const Float_t centmax = 100.;
const Float_t fakemin = -0.5;
const Float_t fakemax = 2.5.;
const Float_t cosminXY = 0.95;
const Float_t cosmaxXY = 1.0;
const Float_t normDecLXYmin = 0;
const Float_t normDecLXYmax = 20;
const Float_t multmin_0_20 = 0;
const Float_t multmax_0_20 = 20;
const Float_t multmin_20_50 = 20;
const Float_t multmax_20_50 = 50;
const Float_t multmin_50_102 = 50;
const Float_t multmax_50_102 = 102;

//----------------------------------------------------

AliCFTaskVertexingHF *AddTaskCFVertexingHF(const char* cutFile = "./D0toKpiCuts.root",Int_t configuration = AliCFTaskVertexingHF::kSnail, Bool_t isKeepDfromB=kFALSE, Bool_t isKeepDfromBOnly=kFALSE, Int_t pdgCode = 421, Char_t isSign = 2)
//AliCFContainer *AddTaskCFVertexingHF(const char* cutFile = "./D0toKpiCuts.root", Int_t configuration = AliCFTaskVertexingHF::kSnail, Bool_t isKeepDfromB=kFALSE, Bool_t isKeepDfromBOnly=kFALSE, Int_t pdgCode = 421, Char_t isSign = 2)
{
	printf("Adding CF task using cuts from file %s\n",cutFile);
	if (configuration == AliCFTaskVertexingHF::kSnail){
		printf("The configuration is set to be SLOW --> all the variables will be used to fill the CF\n");
	}
	else if (configuration == AliCFTaskVertexingHF::kCheetah){
		printf("The configuration is set to be FAST --> using only pt, y, ct, phi, zvtx, centrality, fake, multiplicity to fill the CF\n");
	}
	else{
		printf("The configuration is not defined! returning\n");
		return;
	}
	       
	gSystem->Sleep(2000);

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
	AliRDHFCutsD0toKpi *cutsD0toKpi = (AliRDHFCutsD0toKpi*)fileCuts->Get("D0toKpiCutsStandard");
	
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
	const Double_t phimax = 2*TMath::Pi();
	UInt_t nstep = 10; //number of selection steps: MC with limited acceptance, MC, Acceptance, Vertex, Refit, Reco (no cuts), RecoAcceptance, RecoITSClusters (RecoAcceptance included), RecoPPR (RecoAcceptance+RecoITSCluster included), RecoPID 

	//const UInt_t ipT, iy, icosThetaStar, ipTpi, ipTk, icT, idca, id0xd0, ipointing, iphi, izvtx, icent, ifake, ipointingXY, iNormDecayLXY, imult;
	const Int_t nbiny  = 24 ; //bins in y
	const Int_t nbincosThetaStar  = 42 ; //bins in cosThetaStar 
	const Int_t nbincT  = 15 ; //bins in cT
	const Int_t nbindca  = 20 ; //bins in dca
	const Int_t nbind0xd0  = 90 ; //bins in d0xd0
	const Int_t nbinpointing  = 50 ; //bins in cosPointingAngle
	const Int_t nbinphi  = 18 ; //bins in Phi
	const Int_t nbinzvtx  = 30 ; //bins in z vertex
	const Int_t nbincent = 11;  //bins in centrality
	const Int_t nbincent_0_10 = 2;  //bins in centrality between 0 and 10
	const Int_t nbincent_10_100 = 9;  //bins in centrality between 10 and 100
	const Int_t nbinfake = 3;  //bins in fake
	const Int_t nbinpointingXY = 50;  //bins in cosPointingAngleXY
	const Int_t nbinnormDecayLXY = 20;  //bins in NormDecayLengthXY
	const Int_t nbinmult = 48;  //bins in multiplicity (total number)
	const Int_t nbinmult_0_20 = 20; //bins in multiplicity between 0 and 20
	const Int_t nbinmult_20_50 = 15; //bins in multiplicity between 20 and 50
	const Int_t nbinmult_50_102 = 13; //bins in multiplicity between 50 and 102

	//the sensitive variables, their indices

	const Int_t nvarTot   = 16 ; //number of variables on the grid:pt, y, cosThetaStar, pTpi, pTk, cT, dca, d0pi, d0K, d0xd0, cosPointingAngle, phi, z, centrality, fake, cosPointingAngleXY, normDecayLengthXY, multiplicity

	// variables' indices
	const UInt_t ipT = 0;
	const UInt_t iy  = 1;
	const UInt_t icosThetaStar  = 2;
	const UInt_t ipTpi  = 3;
	const UInt_t ipTk  = 4;
	const UInt_t icT  = 5;
	const UInt_t idca  = 6;
	const UInt_t id0xd0  = 7;
	const UInt_t ipointing  = 8;
	const UInt_t iphi  = 9;
	const UInt_t izvtx  = 10;
	const UInt_t icent = 11;
	const UInt_t ifake = 12;
	const UInt_t ipointingXY = 13;
	const UInt_t inormDecayLXY = 14;
	const UInt_t imult = 15;
	
	//Setting the bins: pt, ptPi, and ptK are considered seprately because for them you can either define the binning by hand, or using the cuts file
	
	//arrays for the number of bins in each dimension
	Int_t iBin[nvarTot];
	
	//OPTION 1: defining the pt, ptPi, ptK bins by hand...		
	/*
	  const Int_t nbinpt_0_6  = 6 ; //bins in pt from 0 to 6 GeV
	  const Int_t nbinpt_6_8  = 1 ; //bins in pt from 6 to 8 GeV
	  const Int_t nbinpt_8_16  = 2 ; //bins in pt from 8 to 16 GeV
	  const Int_t nbinpt_16_24  = 1 ; //bins in pt from 16 to 24 GeV
	  const Int_t nbinpTpi_0_6  = 6 ; //bins in ptPi from 0 to 6 GeV
	  const Int_t nbinpTpi_6_8  = 1 ; //bins in ptPi from 6 to 8 GeV
	  const Int_t nbinpTpi_8_16  = 2 ; //bins in ptPi from 8 to 16 GeV
	  const Int_t nbinpTpi_16_24  = 1 ; //bins in ptPi from 16 to 24 GeV
	  const Int_t nbinpTk_0_6  = 6 ; //bins in ptK from 0 to 6 GeV
	  const Int_t nbinpTk_6_8  = 1 ; //bins in ptK from 6 to 8 GeV
	  const Int_t nbinpTk_8_16  = 2 ; //bins in ptK from 8 to 16 GeV
	  const Int_t nbinpTk_16_24  = 1 ; //bins in ptK from 16 to 24 GeV
	  iBin[ipT]=nbinpt_0_6+nbinpt_6_8+nbinpt_8_16+nbinpt_16_24;
	  iBin[ipTpi]=nbinpTpi_0_6+nbinpTpi_6_8+nbinpTpi_8_16+nbinpTpi_16_24;
	  iBin[ipTk]=nbinpTk_0_6+nbinpTk_6_8+nbinpTk_8_16+nbinpTk_16_24;
	  Double_t *binLimpT=new Double_t[iBin[0]+1];
	  Double_t *binLimpTpi=new Double_t[iBin[3]+1];
	  Double_t *binLimpTk=new Double_t[iBin[4]+1];
	  
	  // values for bin lower bounds
	  // pt
	  for(Int_t i=0; i<=nbinpt_0_6; i++) binLimpT[i]=(Double_t)ptmin_0_6 + (ptmax_0_6-ptmin_0_6)/nbinpt_0_6*(Double_t)i ; 
	  if (binLimpT[nbinpt_0_6] != ptmin_6_8)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 1st range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpt_6_8; i++) binLimpT[i+nbinpt_0_6]=(Double_t)ptmin_6_8 + (ptmax_6_8-ptmin_6_8)/nbinpt_6_8*(Double_t)i ; 
	  if (binLimpT[nbinpt_0_6+nbinpt_6_8] != ptmin_8_16)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpt_8_16; i++) binLimpT[i+nbinpt_0_6+nbinpt_6_8]=(Double_t)ptmin_8_16 + (ptmax_8_16-ptmin_8_16)/nbinpt_8_16*(Double_t)i ; 
	  if (binLimpT[nbinpt_0_6+nbinpt_6_8+nbinpt_8_16] != ptmin_16_24)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpt_16_24; i++) binLimpT[i+nbinpt_0_6+nbinpt_6_8+nbinpt_8_16]=(Double_t)ptmin_16_24 + (ptmax_16_24-ptmin_16_24)/nbinpt_16_24*(Double_t)i ; 
	  
	  // ptPi
	  for(Int_t i=0; i<=nbinpTpi_0_6; i++) binLimpTpi[i]=(Double_t)ptmin_0_6 + (ptmax_0_6-ptmin_0_6)/nbinpTpi_0_6*(Double_t)i ; 
	  if (binLimpTpi[nbinpTpi_0_6] != ptmin_6_8)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 1st range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpTpi_6_8; i++) binLimpTpi[i+nbinpTpi_0_6]=(Double_t)ptmin_6_8 + (ptmax_6_8-ptmin_6_8)/nbinpTpi_6_8*(Double_t)i ; 
	  if (binLimpTpi[nbinpTpi_0_6+nbinpTpi_6_8] != ptmin_8_16)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpTpi_8_16; i++) binLimpTpi[i+nbinpTpi_0_6+nbinpt_6_8]=(Double_t)ptmin_8_16 + (ptmax_8_16-ptmin_8_16)/nbinpTpi_8_16*(Double_t)i ; 
	  if (binLimpTpi[nbinpTpi_0_6+nbinpTpi_6_8+nbinpTpi_8_16] != ptmin_16_24)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpTpi_16_24; i++) binLimpTpi[i+nbinpTpi_0_6+nbinpTpi_6_8+nbinpTpi_8_16]=(Double_t)ptmin_16_24 + (ptmax_16_24-ptmin_16_24)/nbinpTpi_16_24*(Double_t)i ; 
	  
	  // ptKa
	  for(Int_t i=0; i<=nbinpTk_0_6; i++) binLimpTk[i]=(Double_t)ptmin_0_6 + (ptmax_0_6-ptmin_0_6)/nbinpTk_0_6*(Double_t)i ; 
	  if (binLimpTk[nbinpTk_0_6] != ptmin_6_8)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 1st range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpTk_6_8; i++) binLimpTk[i+nbinpTk_0_6]=(Double_t)ptmin_6_8 + (ptmax_6_8-ptmin_6_8)/nbinpTk_6_8*(Double_t)i ; 
	  if (binLimpTk[nbinpTk_0_6+nbinpTk_6_8] != ptmin_8_16)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpTk_8_16; i++) binLimpTk[i+nbinpTk_0_6+nbinpt_6_8]=(Double_t)ptmin_8_16 + (ptmax_8_16-ptmin_8_16)/nbinpTk_8_16*(Double_t)i ; 
	  if (binLimpTk[nbinpTk_0_6+nbinpTk_6_8+nbinpTk_8_16] != ptmin_16_24)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpTk_16_24; i++) binLimpTk[i+nbinpTk_0_6+nbinpTk_6_8+nbinpTk_8_16]=(Double_t)ptmin_16_24 + (ptmax_16_24-ptmin_16_24)/nbinpTk_16_24*(Double_t)i ; 
	*/
	
	//OPTION 2: ...or from the cuts file
	
	const Int_t nbinpt = cutsD0toKpi->GetNPtBins(); // bins in pT
	iBin[ipT]=nbinpt;
	iBin[ipTpi]=nbinpt;
	iBin[ipTk]=nbinpt;
	Double_t *binLimpT=new Double_t[iBin[ipT]+1];
	Double_t *binLimpTpi=new Double_t[iBin[ipTpi]+1];
	Double_t *binLimpTk=new Double_t[iBin[ipTk]+1];
	// values for bin lower bounds
	Float_t* floatbinLimpT = cutsD0toKpi->GetPtBinLimits();
	for (Int_t ibin0 = 0 ; ibin0<iBin[ipT]+1; ibin0++){
		binLimpT[ibin0] = (Double_t)floatbinLimpT[ibin0];
		binLimpTpi[ibin0] = (Double_t)floatbinLimpT[ibin0];
		binLimpTk[ibin0] = (Double_t)floatbinLimpT[ibin0];
	}
	for(Int_t i=0; i<=nbinpt; i++) printf("binLimpT[%d]=%f\n",i,binLimpT[i]);  
	
	printf("pT: nbin (from cuts file) = %d\n",nbinpt);
	
	// defining now the binning for the other variables:
	
	iBin[iy]=nbiny;
	iBin[icosThetaStar]=nbincosThetaStar;
	iBin[icT]=nbincT;
	iBin[idca]=nbindca;
	iBin[id0xd0]=nbind0xd0;
	iBin[ipointing]=nbinpointing;
	iBin[iphi]=nbinphi;
	iBin[izvtx]=nbinzvtx;
	iBin[icent]=nbincent;
	iBin[ifake]=nbinfake;
	iBin[ipointingXY]=nbinpointingXY;
	iBin[inormDecayLXY]=nbinnormDecayLXY;
	iBin[imult]=nbinmult;
	
	//arrays for lower bounds :
	Double_t *binLimy=new Double_t[iBin[iy]+1];
	Double_t *binLimcosThetaStar=new Double_t[iBin[icosThetaStar]+1];
	Double_t *binLimcT=new Double_t[iBin[icT]+1];
	Double_t *binLimdca=new Double_t[iBin[idca]+1];
	Double_t *binLimd0xd0=new Double_t[iBin[id0xd0]+1];
	Double_t *binLimpointing=new Double_t[iBin[ipointing]+1];
	Double_t *binLimphi=new Double_t[iBin[iphi]+1];
	Double_t *binLimzvtx=new Double_t[iBin[izvtx]+1];
	Double_t *binLimcent=new Double_t[iBin[icent]+1];
	Double_t *binLimfake=new Double_t[iBin[ifake]+1];
	Double_t *binLimpointingXY=new Double_t[iBin[ipointingXY]+1];
	Double_t *binLimnormDecayLXY=new Double_t[iBin[inormDecayLXY]+1];
	Double_t *binLimmult=new Double_t[iBin[imult]+1];

	// y
	for(Int_t i=0; i<=nbiny; i++) binLimy[i]=(Double_t)ymin  + (ymax-ymin)  /nbiny*(Double_t)i ;

	// cosThetaStar
	for(Int_t i=0; i<=nbincosThetaStar; i++) binLimcosThetaStar[i]=(Double_t)cosminTS  + (cosmaxTS-cosminTS)  /nbincosThetaStar*(Double_t)i ;
	
	// cT
	for(Int_t i=0; i<=nbincT; i++) binLimcT[i]=(Double_t)cTmin  + (cTmax-cTmin)  /nbincT*(Double_t)i ;

	// dca
	for(Int_t i=0; i<=nbindca; i++) binLimdca[i]=(Double_t)dcamin  + (dcamax-dcamin)  /nbindca*(Double_t)i ;

	// d0xd0
	for(Int_t i=0; i<=nbind0xd0; i++) binLimd0xd0[i]=(Double_t)d0xd0min  + (d0xd0max-d0xd0min)  /nbind0xd0*(Double_t)i ;

	// cosPointingAngle
	for(Int_t i=0; i<=nbinpointing; i++) binLimpointing[i]=(Double_t)cosmin  + (cosmax-cosmin)  /nbinpointing*(Double_t)i ;

	// Phi
	for(Int_t i=0; i<=nbinphi; i++) binLimphi[i]=(Double_t)phimin  + (phimax-phimin)  /nbinphi*(Double_t)i ;

	// z Primary Vertex
	for(Int_t i=0; i<=nbinzvtx; i++) {
		binLimzvtx[i]=(Double_t)zmin  + (zmax-zmin)  /nbinzvtx*(Double_t)i ;
	}

	// centrality
	for(Int_t i=0; i<=nbincent_0_10; i++) binLimcent[i]=(Double_t)centmin_0_10 + (centmax_0_10-centmin_0_10)/nbincent_0_10*(Double_t)i ; 
	if (binLimcent[nbincent_0_10] != centmin_10_100)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for cent - 1st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbincent_10_100; i++) binLimcent[i+nbincent_0_10]=(Double_t)centmin_10_100 + (centmax_10_100-centmin_10_100)/nbincent_10_100*(Double_t)i ; 

	// fake
	for(Int_t i=0; i<=nbinfake; i++) {
	  binLimfake[i]=(Double_t)fakemin  + (fakemax-fakemin)/nbinfake * (Double_t)i;
	}

	// cosPointingAngleXY
	for(Int_t i=0; i<=nbinpointingXY; i++) binLimpointingXY[i]=(Double_t)cosminXY  + (cosmaxXY-cosminXY)  /nbinpointingXY*(Double_t)i ;

	// normDecayLXY
	for(Int_t i=0; i<=nbinnormDecayLXY; i++) binLimnormDecayLXY[i]=(Double_t)normDecLXYmin  + (normDecLXYmax-normDecLXYmin)  /nbinnormDecayLXY*(Double_t)i ;

	// multiplicity
	for(Int_t i=0; i<=nbinmult_0_20; i++) binLimmult[i]=(Double_t)multmin_0_20 + (multmax_0_20-multmin_0_20)/nbinmult_0_20*(Double_t)i ; 
	if (binLimmult[nbinmult_0_20] != multmin_20_50)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for mult - 1st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbinmult_20_50; i++) binLimmult[i+nbinmult_0_20]=(Double_t)multmin_20_50 + (multmax_20_50-multmin_20_50)/nbinmult_20_50*(Double_t)i ; 
	if (binLimmult[nbinmult_0_20+nbinmult_20_50] != multmin_50_102)  {
		Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for mult - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbinmult_50_102; i++) binLimmult[i+nbinmult_0_20+nbinmult_20_50]=(Double_t)multmin_50_102 + (multmax_50_102-multmin_50_102)/nbinmult_50_102*(Double_t)i ; 

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

	//Setting up the container grid... 

	AliCFContainer* container;

	if (configuration == AliCFTaskVertexingHF::kSnail){
		container = new AliCFContainer(nameContainer,"container for tracks",nstep,nvarTot,iBin);
		//setting the bin limits
		printf("pt\n");
		container -> SetBinLimits(ipT,binLimpT);
		printf("y\n");
		container -> SetBinLimits(iy,binLimy);
		printf("cts\n");
		container -> SetBinLimits(icosThetaStar,binLimcosThetaStar);
		printf("ptPi\n");
		container -> SetBinLimits(ipTpi,binLimpTpi);
		printf("ptK\n");
		container -> SetBinLimits(ipTk,binLimpTk);
		printf("cT\n");
		container -> SetBinLimits(icT,binLimcT);
		printf("dca\n");
		container -> SetBinLimits(idca,binLimdca);
		printf("d0xd0\n");
		container -> SetBinLimits(id0xd0,binLimd0xd0);
		printf("pointing\n");
		container -> SetBinLimits(ipointing,binLimpointing);
		printf("phi\n");
		container -> SetBinLimits(iphi,binLimphi);
		printf("z\n");
		container -> SetBinLimits(izvtx,binLimzvtx);
		printf("cent\n");
		container -> SetBinLimits(icent,binLimcent);
		printf("fake\n");
		container -> SetBinLimits(ifake,binLimfake);
		printf("pointingXY\n");
		container -> SetBinLimits(ipointingXY,binLimpointingXY);
		printf("normDecayLXY\n");
		container -> SetBinLimits(inormDecayLXY,binLimnormDecayLXY);
		printf("multiplicity\n");
		container -> SetBinLimits(imult,binLimmult);

		container -> SetVarTitle(ipT,"pt");
		container -> SetVarTitle(iy,"y");
		container -> SetVarTitle(icosThetaStar, "cosThetaStar");
		container -> SetVarTitle(ipTpi, "ptpi");
		container -> SetVarTitle(ipTk, "ptK");
		container -> SetVarTitle(icT, "ct");
		container -> SetVarTitle(idca, "dca");
		container -> SetVarTitle(id0xd0, "d0xd0");
		container -> SetVarTitle(ipointing, "pointing");
		container -> SetVarTitle(iphi, "phi");
		container -> SetVarTitle(izvtx, "zvtx");
		container -> SetVarTitle(icent, "centrality");
		container -> SetVarTitle(ifake, "fake");
		container -> SetVarTitle(ipointingXY, "piointingXY");
		container -> SetVarTitle(inormDecayLXY, "normDecayLXY");
		container -> SetVarTitle(imult, "multiplicity");
	}
	else if (configuration == AliCFTaskVertexingHF::kCheetah){
		//arrays for the number of bins in each dimension
		const Int_t nvar = 8;

		const UInt_t ipTFast = 0;
		const UInt_t iyFast = 1;
		const UInt_t icTFast = 2;
		const UInt_t iphiFast = 3;
		const UInt_t izvtxFast = 4;
		const UInt_t icentFast = 5;
		const UInt_t ifakeFast = 6;
		const UInt_t imultFast = 7;

		Int_t iBinFast[nvar];
		iBinFast[ipTFast] = iBin[ipT];
		iBinFast[iyFast] = iBin[iy];
		iBinFast[icTFast] = iBin[icT];
		iBinFast[iphiFast] = iBin[iphi];
		iBinFast[izvtxFast] = iBin[izvtx];
		iBinFast[icentFast] = iBin[icent];
		iBinFast[ifakeFast] = iBin[ifake];
		iBinFast[imultFast] = iBin[imult];

		container = new AliCFContainer(nameContainer,"container for tracks",nstep,nvar,iBinFast);
		printf("pt\n");
		container -> SetBinLimits(ipTFast,binLimpT);
		printf("y\n");
		container -> SetBinLimits(iyFast,binLimy);
		printf("ct\n");
		container -> SetBinLimits(icTFast,binLimcT);
		printf("phi\n");
		container -> SetBinLimits(iphiFast,binLimphi);
		printf("zvtx\n");
		container -> SetBinLimits(izvtxFast,binLimzvtx);
		printf("centrality\n");
		container -> SetBinLimits(icentFast,binLimcent);
		printf("fake\n");
		container -> SetBinLimits(ifakeFast,binLimfake);
		printf("multiplicity\n");
		container -> SetBinLimits(imultFast,binLimmult);

		container -> SetVarTitle(ipTFast,"pt");
		container -> SetVarTitle(iyFast,"y");
		container -> SetVarTitle(icTFast, "ct");
		container -> SetVarTitle(iphiFast, "phi");
		container -> SetVarTitle(izvtxFast, "zvtx");
		container -> SetVarTitle(icentFast, "centrality");
		container -> SetVarTitle(ifakeFast, "fake");
		container -> SetVarTitle(imultFast, "multiplicity");
	}

	container -> SetStepTitle(0, "MCLimAcc");
	container -> SetStepTitle(1, "MC");
        container -> SetStepTitle(2, "MCAcc");
        container -> SetStepTitle(3, "RecoVertex");
        container -> SetStepTitle(4, "RecoRefit");
        container -> SetStepTitle(5, "Reco");
        container -> SetStepTitle(6, "RecoAcc");
	container -> SetStepTitle(7, "RecoITSCluster");
	container -> SetStepTitle(8, "RecoCuts");
	container -> SetStepTitle(9, "RecoPID");

	//return container;

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
	task->SetConfiguration(configuration);
	task->SetFillFromGenerated(kFALSE);
	task->SetCFManager(man); //here is set the CF manager
	task->SetDecayChannel(2);
	task->SetUseWeight(kFALSE);
	task->SetSign(isSign);
	task->SetCentralitySelection(kFALSE);
	task->SetFakeSelection(0);
	task->SetRejectCandidateIfNotFromQuark(kTRUE); // put to false if you want to keep HIJING D0!!
	task->SetUseMCVertex(kFALSE); // put to true if you want to do studies on pp

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
	Printf("RejectCandidateIfNotFromQuark selection = %d",(Int_t)task->GetRejectCandidateIfNotFromQuark());
	Printf("UseMCVertex selection = %d",(Int_t)task->GetUseMCVertex());
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

        binEdges[0]= binLimpT;
        binEdges[1]= binLimy;

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

