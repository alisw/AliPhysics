


//----------------------------------------------------

AliCFTaskVertexingHF *AddTaskCFVertexingHF3Prong(TString suffixName="", const char* cutFile = "./DplustoKpipiCuts.root", Int_t configuration = AliCFTaskVertexingHF::kCheetah, Bool_t isKeepDfromB=kFALSE, Bool_t isKeepDfromBOnly=kFALSE, Int_t pdgCode = 411, Char_t isSign = 2, Bool_t useWeight=kFALSE, TString multFile="", Bool_t useNchWeight=kFALSE, Bool_t useNtrkWeight=kFALSE, TString estimatorFilename="", Int_t multiplicityEstimator = AliCFTaskVertexingHF::kNtrk10, Bool_t isPPbData = kFALSE, Double_t refMult=9.26, Bool_t isFineNtrkBin=kFALSE, TString multweighthistoname = "hNtrUnCorrEvWithCandWeight")
//AliCFContainer *AddTaskCFVertexingHF3Prong(const char* cutFile = "./DplustoKpipiCuts.root", Int_t configuration = AliCFTaskVertexingHF::kSnail, Bool_t isKeepDfromB=kFALSE, Bool_t isKeepDfromBOnly=kFALSE, Int_t pdgCode = 411, Char_t isSign = 2)
{

//DEFINITION OF A FEW CONSTANTS
const Double_t ymin  = -1.2 ;
const Double_t ymax  =  1.2 ;
const Double_t cosmin = -0.7;
const Double_t cosmax =  1.05;
const Double_t cTmin = 0;  // micron
const Double_t cTmax = 500;  // micron
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
const Double_t zvtxmin = -15;
const Double_t zvtxmax = 15;
const Int_t    minITSClusters = 5;

const Float_t centmin_0_10 = 0.;
const Float_t centmax_0_10 = 10.;
const Float_t centmin_10_60 = 10.;
const Float_t centmax_10_60 = 60.;
const Float_t centmin_60_100 = 60.;
const Float_t centmax_60_100 = 100.;
const Float_t centmax = 100.;
const Float_t fakemin = -0.5;
const Float_t fakemax = 2.5;
const Float_t cosminXY = 0.95;
const Float_t cosmaxXY = 1.0;
const Float_t normDecLXYmin = 0;
const Float_t normDecLXYmax = 20;
const Float_t multmin_0_20 = 0;
const Float_t multmax_0_20 = 20;
const Float_t multmin_20_50 = 20;
const Float_t multmax_20_50 = 50;
const Float_t multmin_50_80 = 50;
const Float_t multmax_50_80 = 80;
const Float_t multmin_80_100 = 80;
const Float_t multmax_80_100 = 100;
const Float_t multmin_100_400 = 100;
const Float_t multmax_100_400 = 400;

	printf("Addig CF task using cuts from file %s\n",cutFile);
	if (configuration == AliCFTaskVertexingHF::kSnail){
		printf("The configuration is set to be SLOW --> all the variables will be used to fill the CF\n");
	}
	else if (configuration == AliCFTaskVertexingHF::kCheetah){
		printf("The configuration is set to be FAST --> using only pt, y, ct, phi, zvtx, centrality, fake, multiplicity to fill the CF\n");
	}
  else if (configuration == AliCFTaskVertexingHF::kFalcon){
    printf("The configuration is set to be FAST --> using only pt, y, centrality, multiplicity to fill the CF\n");
  }
  else if (configuration == AliCFTaskVertexingHF::kESE){
    printf("The configuration is set to be for ESE analysis --> using pt, y, centrality, multiplicity, local multiplicity and q2 to fill the CF\n");
  }
  else if (configuration == AliCFTaskVertexingHF::kRT) {
    printf("The configuration is set to be for RT analysis --> using pt, y, multiplicity, RT, delta-phi leading to fill the CF\n");
  }
	else{
		printf("The configuration is not defined! returning\n");
		return NULL;
	}

	gSystem->Sleep(2000);

	// isSign = 0 --> D+ only
	// isSign = 1 --> D- only
	// isSign = 2 --> both

	TString expected;
	if (isSign == 0 && pdgCode < 0){
		Printf("ERROR: Error setting PDG code (%d) and sign (0 --> particle (%d) only): they are not compatible, returning",pdgCode,isSign);
		return 0x0;
	}
	else if (isSign == 1 && pdgCode > 0){
		Printf("ERROR: Error setting PDG code (%d) and sign (1 --> antiparticle (%d) only): they are not compatible, returning",pdgCode,isSign);
		return 0x0;
	}
	else if (isSign > 2 || isSign < 0){
		Printf("ERROR: Sign not valid (%d, possible values are 0, 1, 2), returning",isSign);
		return 0x0;
	}

	TFile* fileCuts = TFile::Open(cutFile);
	if(!fileCuts || (fileCuts && !fileCuts->IsOpen())){
	  Printf("FATAL:  Cut file not found");
	  return 0x0;
	}
	AliRDHFCutsDplustoKpipi *cutsDplustoKpipi = (AliRDHFCutsDplustoKpipi*)fileCuts->Get("AnalysisCuts");
	TH1F *hMult=0x0;
	if(multFile.EqualTo("") ) {
	  printf("Will not be corrected with weights \n");
	}else{
	  TFile *fileMult = TFile::Open(multFile.Data());
	  if(isPPbData && multweighthistoname=="") {
			hMult = (TH1F*)fileMult->Get("hNtrUnCorrEvWithCandWeight");
		}
	  else if(!isPPbData && multweighthistoname==""){
	    hMult=(TH1F*)fileMult->Get("hGenPrimaryParticlesInelGt0");
	  }
		else {
			hMult=(TH1F*)fileMult->Get(multweighthistoname.Data());
		}
	}

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

	const Double_t phimax = 2*TMath::Pi();

	//Setting up the container grid...
	UInt_t nstep = 10; //number of selection steps: MC with limited acceptance, MC, Acceptance, Vertex, Refit, Reco (no cuts), RecoAcceptance, RecoITSClusters (RecoAcceptance included), RecoPPR (RecoAcceptance+RecoITSCluster included), RecoPID
// 	const Int_t nbinpt_0_4  = 8 ; //bins in pt from 0 to 4 GeV
// 	const Int_t nbinpt_4_8  = 4 ; //bins in pt from 4 to 8 GeV
// 	const Int_t nbinpt_8_10  = 1 ; //bins in pt from 8 to 10 GeV

/*
	Int_t nbinpt_0_4;
	Int_t nbinpt_4_8;
	Int_t nbinpt_8_10;
	if (!isKeepDfromB){
	  nbinpt_0_4  = 8 ; //bins in pt from 0 to 4 GeV
	  nbinpt_4_8  = 4 ; //bins in pt from 4 to 8 GeV
	  nbinpt_8_10  = 1 ; //bins in pt from 8 to 10 GeV
	}else{
	  nbinpt_0_4  = 3 ; //bins in pt from 0 to 3 GeV
	  nbinpt_4_8  = 1 ; //bins in pt from 3 to 5 GeV
	  nbinpt_8_10  = 1 ; //bins in pt from 5 to 10 GeV
	}
*/
	const Int_t nbinpt = cutsDplustoKpipi->GetNPtBins(); // bins in pT
	printf("pT: nbin (from cuts file) = %d\n",nbinpt);
	const Int_t nbiny  = 24 ; //bins in y
	const Int_t nbinphi  = 18 ; //bins in phi
	const Int_t nbincT  = 2 ; //bins in cT
	const Int_t nbinpointing  = 35 ; //bins in cosPointingAngle
	const Int_t nbinpTpi_0_4  = 8 ; //bins in ptPi from 0 to 4 GeV
	const Int_t nbinpTpi_4_8  = 4 ; //bins in ptPi from 4 to 8 GeV
	const Int_t nbinpTpi_8_10  = 1 ; //bins in ptPi from 8 to 10 GeV
	const Int_t nbinpTk_0_4  = 8 ; //bins in ptKa from 0 to 4 GeV
	const Int_t nbinpTk_4_8  = 4 ; //bins in ptKa from 4 to 8 GeV
	const Int_t nbinpTk_8_10  = 1 ; //bins in ptKa from 8 to 10 GeV
	const Int_t nbinpTpi2_0_4  = 8 ; //bins in ptpi2 from 0 to 4 GeV
	const Int_t nbinpTpi2_4_8  = 4 ; //bins in ptpi2 from 4 to 8 GeV
	const Int_t nbinpTpi2_8_10  = 1 ; //bins in ptpi2 from 8 to 10 GeV
	const Int_t nbinzvtx  = 30 ; //bins in z vertex
	const Int_t nbincent = 18; //bins in centrality
	const Int_t nbincent_0_10 = 4;  //bins in centrality between 0 and 10
	const Int_t nbincent_10_60 = 10;  //bins in centrality between 10 and 60
	const Int_t nbincent_60_100 = 4;  //bins in centrality between 60 and 100
	const Int_t nbinfake = 3;  //bins in fake
	const Int_t nbinpointingXY = 50;  //bins in cosPointingAngleXY
	const Int_t nbinnormDecayLXY = 20;  //bins in NormDecayLengthXY
	Int_t nbinmult = 49;  //bins in multiplicity (total number)
	const Int_t nbinmult_0_20 = 20; //bins in multiplicity between 0 and 20
	const Int_t nbinmult_20_50 = 15; //bins in multiplicity between 20 and 50
	const Int_t nbinmult_50_80 = 10; //bins in multiplicity between 50 and 80
	const Int_t nbinmult_80_100 = 4; //bins in multiplicity between 80 and 100
	const Int_t nbinmult_100_400 = 3; //bins in multiplicity between 100 and 400
	if(isPPbData) nbinmult += nbinmult_100_400;

	// Fine Ntrk bining setting
	Double_t *binLimmultFine;
	Int_t nbinmultTmp=nbinmult;
	if(isFineNtrkBin){
	  Int_t nbinLimmultFine=100;
	  if(isPPbData) nbinLimmultFine = 200;
	  const UInt_t nbinMultFine = nbinLimmultFine;
	  binLimmultFine = new Double_t[nbinMultFine+1];
	  for (Int_t ibin0 = 0 ; ibin0<nbinMultFine+1; ibin0++){
	    binLimmultFine[ibin0] = ibin0;
	  }
	  nbinmultTmp=nbinLimmultFine;
	}
	const Int_t nbinmultTot=nbinmultTmp;


	//the sensitive variables, their indices
	const UInt_t ipT = 0;
	const UInt_t iy  = 1;
	const UInt_t iphi  = 2;
	const UInt_t icT  = 3;
	const UInt_t ipointing  = 4;
	const UInt_t ipTpi  = 5;
	const UInt_t ipTk  = 6;
	const UInt_t ipTpi2  = 7;
	const UInt_t izvtx  = 8;
	const UInt_t icent = 9;
	const UInt_t ifake = 10;
	const UInt_t ipointingXY = 11;
	const UInt_t inormDecayLXY = 12;
	const UInt_t imult = 13;

	const Int_t nvarTot   = 14 ; //number of variables on the grid:pt, y, cosThetaStar, pTpi, pTk, cT, dca, d0pi, d0K, d0xd0, cosPointingAngle, phi, zvtx, centrality, fake, cosPointingAngleXY, normDecayLengthXY, multiplicity

	//arrays for the number of bins in each dimension
	Int_t iBin[nvarTot];
 	//iBin[ipT]=nbinpt_0_4+nbinpt_4_8+nbinpt_8_10;
	iBin[ipT]=nbinpt;
	iBin[iy]=nbiny;
	iBin[iphi]=nbinphi;
	// 	iBin[icT]=nbincT_0_4+nbincT_4_8+nbincT_8_10;
 	//iBin[4]=nbinpointing_0_4+nbinpointing_4_8+nbinpointing_8_10;
 	iBin[icT]=nbincT;
 	iBin[ipointing]=nbinpointing;
	iBin[ipTpi]=nbinpt;
	iBin[ipTk]=nbinpt;
	iBin[ipTpi2]=nbinpt;
	iBin[izvtx]=nbinzvtx;
	iBin[icent]=nbincent;
	iBin[ifake]=nbinfake;
	iBin[ipointingXY]=nbinpointingXY;
	iBin[inormDecayLXY]=nbinnormDecayLXY;
	iBin[imult]=nbinmultTot;

	//arrays for lower bounds :
	Double_t *binLimpT=new Double_t[iBin[ipT]+1];
	Double_t *binLimy=new Double_t[iBin[iy]+1];
	Double_t *binLimphi=new Double_t[iBin[iphi]+1];
	Double_t *binLimcT=new Double_t[iBin[icT]+1];
	Double_t *binLimpointing=new Double_t[iBin[ipointing]+1];
	Double_t *binLimpTpi=new Double_t[iBin[ipTpi]+1];
	Double_t *binLimpTk=new Double_t[iBin[ipTk]+1];
	Double_t *binLimpTpi2=new Double_t[iBin[ipTpi2]+1];
	Double_t *binLimzvtx=new Double_t[iBin[izvtx]+1];
	Double_t *binLimcent=new Double_t[iBin[icent]+1];
	Double_t *binLimfake=new Double_t[iBin[ifake]+1];
	Double_t *binLimpointingXY=new Double_t[iBin[ipointingXY]+1];
	Double_t *binLimnormDecayLXY=new Double_t[iBin[inormDecayLXY]+1];
	Double_t *binLimmult=new Double_t[iBin[imult]+1];

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
	Float_t* floatbinLimpT = cutsDplustoKpipi->GetPtBinLimits();
	for (Int_t ibinpT = 0 ; ibinpT<iBin[ipT]+1; ibinpT++){
		binLimpT[ibinpT] = (Double_t)floatbinLimpT[ibinpT];
		binLimpTpi[ibinpT] = (Double_t)floatbinLimpT[ibinpT];
		binLimpTk[ibinpT] = (Double_t)floatbinLimpT[ibinpT];
		binLimpTpi2[ibinpT] = (Double_t)floatbinLimpT[ibinpT];
	}
	for(Int_t i=0; i<=nbinpt; i++) printf("binLimpT[%d]=%f\n",i,binLimpT[i]);

	/*
	  for(Int_t i=0; i<=nbinpt_0_4; i++) binLimpT[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbinpt_0_4*(Double_t)i ;
	  if (binLimpT[nbinpt_0_4] != ptmin_4_8)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 1st range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpt_4_8; i++) binLimpT[i+nbinpt_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbinpt_4_8*(Double_t)i ;
	  if (binLimpT[nbinpt_0_4+nbinpt_4_8] != ptmin_8_10)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for pt - 2nd range - differs from expected!\n");
	  }
	  for(Int_t i=0; i<=nbinpt_8_10; i++) binLimpT[i+nbinpt_0_4+nbinpt_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbinpt_8_10*(Double_t)i ;
	*/

	// y
	for(Int_t i=0; i<=nbiny; i++) binLimy[i]=(Double_t)ymin  + (ymax-ymin)  /nbiny*(Double_t)i ;

	// Phi
	for(Int_t i=0; i<=nbinphi; i++) binLimphi[i]=(Double_t)phimin  + (phimax-phimin)  /nbinphi*(Double_t)i ;

	// cT
	for(Int_t i=0; i<=nbincT; i++) binLimcT[i]=(Double_t)cTmin  + (cTmax-cTmin)  /nbincT*(Double_t)i ;

	// cosPointingAngle
	for(Int_t i=0; i<=nbinpointing; i++) binLimpointing[i]=(Double_t)cosmin  + (cosmax-cosmin)  /nbinpointing*(Double_t)i ;

	/*
	// ptPi
	for(Int_t i=0; i<=nbincT_0_4; i++) binLimcT[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbincT_0_4*(Double_t)i ;
	if (binLimcT[nbincT_0_4] != ptmin_4_8)  {
	Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for ptPi - 1st range - differs from expected!");
	}
	for(Int_t i=0; i<=nbincT_4_8; i++) binLimcT[i+nbincT_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbincT_4_8*(Double_t)i ;
	if (binLimcT[nbincT_0_4+nbincT_4_8] != ptmin_8_10)  {
	Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for ptPi - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbincT_8_10; i++) binLimcT[i+nbincT_0_4+nbincT_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbincT_8_10*(Double_t)i ;

	// ptKa
	for(Int_t i=0; i<=nbinpointing_0_4; i++) binLimpointing[i]=(Double_t)ptmin_0_4 + (ptmax_0_4-ptmin_0_4)/nbinpointing_0_4*(Double_t)i ;
	if (binLimpointing[nbinpointing_0_4] != ptmin_4_8)  {
	Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for ptKa - 1st range - differs from expected!");
	}
	for(Int_t i=0; i<=nbinpointing_4_8; i++) binLimpointing[i+nbinpointing_0_4]=(Double_t)ptmin_4_8 + (ptmax_4_8-ptmin_4_8)/nbinpointing_4_8*(Double_t)i ;
	if (binLimpointing[nbinpointing_0_4+nbinpointing_4_8] != ptmin_8_10)  {
	Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for ptKa - 2nd range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbinpointing_8_10; i++) binLimpointing[i+nbinpointing_0_4+nbinpointing_4_8]=(Double_t)ptmin_8_10 + (ptmax_8_10-ptmin_8_10)/nbinpointing_8_10*(Double_t)i ;
	*/

	// z Primary Vertex
	for(Int_t i=0; i<=nbinzvtx; i++) {
		binLimzvtx[i]=(Double_t)zvtxmin  + (zvtxmax-zvtxmin)  /nbinzvtx*(Double_t)i ;
	}

	// centrality

	for(Int_t i=0; i<=nbincent_0_10; i++) binLimcent[i]=(Double_t)centmin_0_10 + (centmax_0_10-centmin_0_10)/nbincent_0_10*(Double_t)i ;
	if (binLimcent[nbincent_0_10] != centmin_10_60)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for cent - 1st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbincent_10_60; i++) binLimcent[i+nbincent_0_10]=(Double_t)centmin_10_60 + (centmax_10_60-centmin_10_60)/nbincent_10_60*(Double_t)i ;
	if (binLimcent[nbincent_0_10+nbincent_10_60] != centmin_60_100)  {
	  Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for cent - 2st range - differs from expected!\n");
	}
	for(Int_t i=0; i<=nbincent_60_100; i++) binLimcent[i+nbincent_0_10+nbincent_10_60]=(Double_t)centmin_60_100 + (centmax_60_100-centmin_60_100)/nbincent_60_100*(Double_t)i ;

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
      if (binLimmult[nbinmult_0_20+nbinmult_20_50] != multmin_50_80)  {
               Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for mult - 2nd range - differs from expected!\n");
      }
      for(Int_t i=0; i<=nbinmult_50_80; i++) binLimmult[i+nbinmult_0_20+nbinmult_20_50]=(Double_t)multmin_50_80 + (multmax_50_80-multmin_50_80)/nbinmult_50_80*(Double_t)i ;
      if (binLimmult[nbinmult_0_20+nbinmult_20_50+nbinmult_50_80] != multmin_80_100)  {
               Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for mult - 3rd range - differs from expected!\n");
      }
      for(Int_t i=0; i<=nbinmult_80_100; i++) binLimmult[i+nbinmult_0_20+nbinmult_20_50+nbinmult_50_80]=(Double_t)multmin_80_100 + (multmax_80_100-multmin_80_100)/nbinmult_80_100*(Double_t)i ;
      if (binLimmult[nbinmult_0_20+nbinmult_20_50+nbinmult_50_80+nbinmult_80_100] != multmin_100_400) {
               Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for mult - 4th range - differs from expected!\n");
      }
      if(isPPbData){
	for (Int_t i = 0; i<=nbinmult_100_400; i++) binLimmult[i+nbinmult_0_20+nbinmult_20_50+nbinmult_50_80+nbinmult_80_100]= (Double_t)multmin_100_400 + (multmax_100_400-multmin_100_400)/nbinmult_100_400*(Double_t)i ;
      }

	//one "container" for MC
	TString nameContainer="";
	if(!isKeepDfromB) {
		nameContainer="CFHFcontainer_DplustoKpipi_Prompt";

	}
	else  if(isKeepDfromBOnly){
		nameContainer="CFHFcontainer_DplustoKpipi_FromB";

	}
	else  {
	  nameContainer="CFHFcontainer_DplustoKpipi_All";
	}
	nameContainer += suffixName.Data();


	AliCFContainer* container;
	if (configuration == AliCFTaskVertexingHF::kSnail){
		container = new AliCFContainer(nameContainer,"container for tracks",nstep,nvarTot,iBin);
		//setting the bin limits
		printf("pt\n");
		container -> SetBinLimits(ipT,binLimpT);
		printf("y\n");
		container -> SetBinLimits(iy,binLimy);
		printf("Phi\n");
		container -> SetBinLimits(iphi,binLimphi);
		printf("cT\n");
		container -> SetBinLimits(icT,binLimcT);
		printf("pointing angle\n");
		container -> SetBinLimits(ipointing,binLimpointing);
		printf("ptpi\n");
		container -> SetBinLimits(ipTpi,binLimpTpi);
		printf("ptK\n");
		container -> SetBinLimits(ipTk,binLimpTk);
		printf("ptpi2\n");
		container -> SetBinLimits(ipTpi2,binLimpTpi2);
		printf("zvtx \n");
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
		if(isFineNtrkBin) container -> SetBinLimits(imult,binLimmultFine);
		else container -> SetBinLimits(imult,binLimmult);

		container -> SetVarTitle(ipT,"pt");
		container -> SetVarTitle(iy,"y");
		container -> SetVarTitle(iphi, "phi");
		container -> SetVarTitle(icT, "ct");
		container -> SetVarTitle(ipointing, "pointing");
		container -> SetVarTitle(ipTpi, "ptpi");
		container -> SetVarTitle(ipTk, "ptK");
		container -> SetVarTitle(ipTpi2, "ptpi2");
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
		if(isFineNtrkBin) container -> SetBinLimits(imultFast,binLimmultFine);
		else container -> SetBinLimits(imultFast,binLimmult);

		container -> SetVarTitle(ipTFast,"pt");
		container -> SetVarTitle(iyFast,"y");
		container -> SetVarTitle(icTFast, "ct");
		container -> SetVarTitle(iphiFast, "phi");
		container -> SetVarTitle(izvtxFast, "zvtx");
		container -> SetVarTitle(icentFast, "centrality");
		container -> SetVarTitle(ifakeFast, "fake");
		container -> SetVarTitle(imultFast, "multiplicity");
	}
	else if (configuration == AliCFTaskVertexingHF::kFalcon){
		//arrays for the number of bins in each dimension
		const Int_t nvar = 4;

		const UInt_t ipTSuperFast = 0;
		const UInt_t iySuperFast = 1;
		const UInt_t icentSuperFast = 2;
		const UInt_t imultSuperFast = 3;

		Int_t iBinSuperFast[nvar];
		iBinSuperFast[ipTSuperFast] = iBin[ipT];
		iBinSuperFast[iySuperFast] = iBin[iy];
		iBinSuperFast[icentSuperFast] = iBin[icent];
		iBinSuperFast[imultSuperFast] = iBin[imult];

		container = new AliCFContainer(nameContainer,"container for tracks",nstep,nvar,iBinSuperFast);
		printf("pt\n");
		container -> SetBinLimits(ipTSuperFast,binLimpT);
		printf("y\n");
		container -> SetBinLimits(iySuperFast,binLimy);
		printf("centrality\n");
		container -> SetBinLimits(icentSuperFast,binLimcent);
		printf("multiplicity\n");
		if(isFineNtrkBin) container -> SetBinLimits(imultSuperFast,binLimmultFine);
		else container -> SetBinLimits(imultSuperFast,binLimmult);

		container -> SetVarTitle(ipTSuperFast,"pt");
		container -> SetVarTitle(iySuperFast,"y");
		container -> SetVarTitle(icentSuperFast, "centrality");
		container -> SetVarTitle(imultSuperFast, "multiplicity");
	}
  else if (configuration == AliCFTaskVertexingHF::kESE){
    //arrays for the number of bins in each dimension
    const Int_t nvar = 6;
    
    const UInt_t ipTESE = 0;
    const UInt_t iyESE = 1;
    const UInt_t icentESE = 2;
    const UInt_t imultESE = 3;
    const UInt_t ilocalmultESE = 4;
    const UInt_t iq2ESE = 5;

    const Int_t iBinESE[nvar] = {iBin[ipT],iBin[iy],100,50,50,100};

    Double_t binLimcentESE[iBinESE[icentESE]+1];
    for(Int_t iCent=0; iCent<iBinESE[icentESE]+1; iCent++) {
      binLimcentESE[iCent] = iCent;
    }
    Double_t binLimmultESE[iBinESE[imultESE]+1];
    for(Int_t iMult=0; iMult<iBinESE[imultESE]+1; iMult++) {
      binLimmultESE[iMult] = -0.5+iMult*5000./iBinESE[imultESE];
    }
    Double_t binLimlocalmultESE[iBinESE[ilocalmultESE]+1];
    for(Int_t iLocalMult=0; iLocalMult<iBinESE[ilocalmultESE]+1; iLocalMult++) {
      binLimlocalmultESE[iLocalMult] = -0.5+iLocalMult*200./iBinESE[ilocalmultESE];
    }
    Double_t binLimq2ESE[iBinESE[iq2ESE]+1];
    for(Int_t iq2=0; iq2<iBinESE[iq2ESE]+1; iq2++) {
      binLimq2ESE[iq2] = iq2*5./iBinESE[iq2ESE];
    }

    container = new AliCFContainer(nameContainer,"container for tracks",nstep,nvar,iBinESE);
    printf("pt\n");
    container -> SetBinLimits(ipTESE,binLimpT);
    printf("y\n");
    container -> SetBinLimits(iyESE,binLimy);
    printf("centrality\n");
    container -> SetBinLimits(icentESE,binLimcentESE);
    printf("multiplicity\n");
    container -> SetBinLimits(imultESE,binLimmultESE);
    printf("local multiplicity\n");
    container -> SetBinLimits(ilocalmultESE,binLimlocalmultESE);
    printf("q2\n");
    container -> SetBinLimits(iq2ESE,binLimq2ESE);

    container -> SetVarTitle(ipTESE,"pt");
    container -> SetVarTitle(iyESE,"y");
    container -> SetVarTitle(icentESE, "centrality");
    container -> SetVarTitle(imultESE, "multiplicity");
    container -> SetVarTitle(ilocalmultESE, "local multiplicity");
    container -> SetVarTitle(iq2ESE, "q2");
  }
 else if (configuration == AliCFTaskVertexingHF::kRT) {
    //arrays for number of bins in each dimension
    const Int_t nvar = 5;

    const UInt_t ipTRT = 0;
    const UInt_t iyRT = 1;
    const UInt_t imultRT = 2;
    const UInt_t iRT = 3;
    const UInt_t idelphiRT = 4;

    const Int_t iBinRT[nvar] = {iBin[ipT], iBin[iy], iBin[imult], 100, 100};

    Double_t binLimRT[iBinRT[iRT]+1];
    for (Int_t jRT = 0; jRT < iBinRT[iRT]+1; jRT++) {
       binLimRT[jRT] = jRT / 10.;
    }
    
    Double_t binLimDeltaPhi[iBinRT[idelphiRT]+1];
    for (Int_t jDelPhi =0; jDelPhi < iBinRT[idelphiRT]+1; jDelPhi++) {
       binLimDeltaPhi[jDelPhi] = -TMath::PiOver2() + (jDelPhi * TMath::TwoPi()/iBinRT[idelphiRT]);
    }
    
    container = new AliCFContainer(nameContainer,"container for tracks",nstep,nvar,iBinRT);
    
    container -> SetBinLimits(ipTRT,binLimpT);
    printf("pt\n");
    container -> SetBinLimits(iyRT,binLimRT);
    printf("y\n");
    container -> SetBinLimits(imultRT,binLimmult);
    printf("multiplicity\n");
    container -> SetBinLimits(iRT,binLimRT);
    printf("RT\n");
    container -> SetBinLimits(idelphiRT,binLimDeltaPhi);
    printf("delta phi leading\n");
    
    container -> SetVarTitle(ipTRT,"pt");
    container -> SetVarTitle(iyRT,"y");
    container -> SetVarTitle(imultRT,"multiplicity");
    container -> SetVarTitle(iRT,"rt");
    container -> SetVarTitle(idelphiRT,"deltaphileading");
    
  }

	//return container;

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
	task->SetConfiguration(configuration);
	task->SetFillFromGenerated(kFALSE);
	task->SetDecayChannel(31);
	task->SetUseWeight(useWeight);
	task->SetCFManager(man); //here is set the CF manager
	task->SetSign(isSign);
	task->SetCentralitySelection(kFALSE);
	task->SetFakeSelection(0);
	task->SetRejectCandidateIfNotFromQuark(kTRUE); // put to false if you want to keep HIJING D0!!
	task->SetUseMCVertex(kFALSE); // put to true if you want to do studies on pp
	task->SetUseNchWeight(useNchWeight); //correction with mult weight
	if(useNchWeight || useNtrkWeight){
	  if(hMult) task->SetMCNchHisto(hMult);
	  else{
	    Printf("FATAL: Histogram for multiplicity weights not found");
	    return 0x0;
	  }
	  task->SetUseNchWeight(kTRUE);
	  if(useNtrkWeight) task->SetUseNchTrackletsWeight();
	}
	if(isPPbData) {
	  task->SetIsPPbData(kTRUE);
	}
	if (isKeepDfromB && !isKeepDfromBOnly) task->SetDselection(2);
	if (isKeepDfromB && isKeepDfromBOnly) task->SetDselection(1);

	TF1* funcWeight = 0x0;
	if (task->GetUseWeight()) {
		funcWeight = (TF1*)fileCuts->Get("funcWeight");
		if (funcWeight == 0x0){
			Printf("FONLL Weights will be used");
		}
		else {
			task->SetWeightFunction(funcWeight);
			Printf("User-defined Weights will be used. The function being:");
			task->GetWeightFunction()->Print();
		}
	}

   task->SetMultiplicityEstimator(multiplicityEstimator);
   if(estimatorFilename.EqualTo("") ) {
     printf("Estimator file not provided, multiplicity corrected histograms will not be filled\n");
     task->SetUseZvtxCorrectedNtrkEstimator(kFALSE);
   } else{
     TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
     if(!fileEstimator)  {
       Printf("FATAL: File with multiplicity estimator not found");
       return NULL;
     }
     task->SetUseZvtxCorrectedNtrkEstimator(kTRUE);
     task->SetReferenceMultiplcity(refMult);

      if (isPPbData) {  //load multiplicity estimators for pPb
         const Char_t* periodNames[2] = {"LHC13b", "LHC13c"};
         TProfile *multEstimatorAvg[2];
         for (Int_t ip=0; ip < 2; ip++) {
            multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("SPDmult10_%s",periodNames[ip]))->Clone(Form("SPDmult10_%s_clone",periodNames[ip])));
            if (!multEstimatorAvg[ip]) {
               Printf("FATAL: Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
               return NULL;
            }
         }
         task->SetMultiplVsZProfileLHC13b(multEstimatorAvg[0]);
         task->SetMultiplVsZProfileLHC13c(multEstimatorAvg[1]);

      } else {    //load multiplicity estimators for pp
         const Char_t* periodNames[4] = {"LHC10b", "LHC10c", "LHC10d", "LHC10e"};
         TProfile* multEstimatorAvg[4];

         for(Int_t ip=0; ip<4; ip++) {
            multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("SPDmult10_%s",periodNames[ip]))->Clone(Form("SPDmult10_%s_clone",periodNames[ip])));
            if (!multEstimatorAvg[ip]) {
               Printf("FATAL: Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
               return NULL;
            }
         }
         task->SetMultiplVsZProfileLHC10b(multEstimatorAvg[0]);
         task->SetMultiplVsZProfileLHC10c(multEstimatorAvg[1]);
         task->SetMultiplVsZProfileLHC10d(multEstimatorAvg[2]);
         task->SetMultiplVsZProfileLHC10e(multEstimatorAvg[3]);
     }
   }

	Printf("***************** CONTAINER SETTINGS *****************");
	Printf("decay channel = %d",(Int_t)task->GetDecayChannel());
	Printf("FillFromGenerated = %d",(Int_t)task->GetFillFromGenerated());
	Printf("Dselection = %d",(Int_t)task->GetDselection());
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

        thnDim[0] = iBin[ipT];
        thnDim[2] = iBin[ipT];
        thnDim[1] = iBin[iy];
        thnDim[3] = iBin[iy];

	TString nameCorr="";
	if(!isKeepDfromB) {
		nameCorr="CFHFcorr_DplustoKpipi_Prompt";
	}
	else  if(isKeepDfromBOnly){
		nameCorr="CFHFcorr_DplustoKpipi_FromB";
	}
	else  {
	        nameCorr="CFHFcorr_DplustoKpipi_All";
	}
	nameCorr += suffixName.Data();


        THnSparseD* correlation = new THnSparseD(nameCorr,"THnSparse with correlations",4,thnDim);
        Double_t** binEdges = new Double_t*[2];

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
	TString output1name="", output2name="", output3name="", output4name="", output5name="";
	output2name=nameContainer;
	output3name=nameCorr;
	output5name= "coutProfDp";
	if(!isKeepDfromB) {
		outputfile += ":PWG3_D2H_CFtaskDplustoKpipi_Prompt";
		output1name="CFHFhist_DplustoKpipi_Prompt";
		output3name+="_Prompt";
		output4name= "Cuts_DplustoKpipi_Prompt";
		output5name+="_Prompt";
	}
	else  if(isKeepDfromBOnly){
   	        outputfile += ":PWG3_D2H_CFtaskDplustoKpipi_FromB";
	        output1name="CFHFhist_DplustoKpipi_FromB";
		output3name+="_FromB";
		output4name= "Cuts_DplustoKpipi_FromB";
		output5name+="_FromB";
	}
	else{
        	outputfile += ":PWG3_D2H_CFtaskDplustoKpipi_All";
		output1name="CFHFhist_DplustoKpipi_All";
		output3name+="_All";
		output4name= "Cuts_DplustoKpipi_All";
		output5name+="_All";
	}
	outputfile += suffixName.Data();
	output1name += suffixName.Data();
	output4name += suffixName.Data();
	output5name += suffixName.Data();

	//now comes user's output objects :
	// output TH1I for event counting
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(output1name, TH1I::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	// output Correction Framework Container (for acceptance & efficiency calculations)
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(output2name, AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	// Unfolding - correlation matrix
        AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(output3name, THnSparseD::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
	AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(output4name, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
	// estimators list
	AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(output5name, TList::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());

	mgr->AddTask(task);

	mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task,1,coutput1);
	mgr->ConnectOutput(task,2,coutput2);
        mgr->ConnectOutput(task,3,coutput3);
	mgr->ConnectOutput(task,4,coutput4);
	mgr->ConnectOutput(task,5,coutput5);

	return task;
}
