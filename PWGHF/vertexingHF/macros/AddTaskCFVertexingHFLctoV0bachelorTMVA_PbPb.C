
//----------------------------------------------------

AliCFTaskVertexingHF *AddTaskCFVertexingHFLctoV0bachelorTMVA_PbPb(Bool_t isCent = kFALSE, const char* cutFile = "CutsLc2pK0STMVA_20140311.root", 
								  TString cutObjectName = "LctoV0AnalysisCuts", 
								  TString suffix = "TMVA", 
								  Int_t configuration = AliCFTaskVertexingHF::kCheetah,
								  Bool_t rejectIfNotFromQuark = kTRUE,
								  Bool_t isKeepDfromB = kTRUE, 
								  Bool_t isKeepDfromBOnly = kFALSE, 
								  Float_t cutOnMomConservation = 1000.,
								  Int_t pdgCode = 4122, 
								  Char_t isSign = 2, 
								  Bool_t useWeight = kFALSE, 
								  Bool_t useFlatPtWeight = kFALSE, 
								  Bool_t useZWeight = kFALSE,
								  const char* weightHistoName = "")
{
const Double_t ymin  = -1.2 ;
const Double_t ymax  =  1.2 ;
const Double_t cosminTS = -1.05;
const Double_t cosmaxTS =  1.05;
const Double_t cosmin = 0.7;
const Double_t cosmax =  1.02;
const Double_t cTmin = 0;  // micron
const Double_t cTmax = 300000;  // micron
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

 printf("Adding CF task using cuts from file %s\n",cutFile);
  if (configuration == AliCFTaskVertexingHF::kSnail){
    printf("The configuration is set to be SLOW --> all the variables will be used to fill the CF\n");
  }
  else if (configuration == AliCFTaskVertexingHF::kCheetah){
    printf("The configuration is set to be FAST --> using only pt, y, ct, phi, zvtx, centrality, fake, multiplicity to fill the CF\n");
  }
  else if (configuration == AliCFTaskVertexingHF::kFalcon){
    printf("The configuration is set to be FAST --> using only pt, y, centrality, multiplicity to fill the CF\n");
  }
  else{
    printf("The configuration is not defined! returning\n");
    return NULL;
  }
	       
  gSystem->Sleep(2000);

  // isSign = 0 --> D* only
  // isSign = 1 --> D*bar only
  // isSign = 2 --> D* + D*bar

  TString expected;
  if (isSign == 0 && pdgCode < 0){
    Printf("ERROR:Error setting PDG code (%d) and sign (0 --> D* only): they are not compatible, returning",pdgCode);
    return 0x0;
  }
  else if (isSign == 1 && pdgCode > 0){
    Printf("ERROR:Error setting PDG code (%d) and sign (1 --> D*bar only): they are not compatible, returning",pdgCode);
    return 0x0;
  }
  else if (isSign > 2 || isSign < 0){
    Printf("ERROR:Sign not valid (%d, possible values are 0, 1, 2), returning",isSign);
    return 0x0;
  }

  TFile* fileCuts = TFile::Open(cutFile);
  if(!fileCuts || (fileCuts && !fileCuts->IsOpen())){ 
    Printf("ERROR: Wrong cut file");
    return 0x0;
  }

  AliRDHFCutsLctoV0 *cutsLcTMVA = (AliRDHFCutsLctoV0*)fileCuts->Get(cutObjectName.Data());
	
  // check that the fKeepD0fromB flag is set to true when the fKeepD0fromBOnly flag is true
  //  for now the binning is the same than for all D's
  if (isKeepDfromBOnly) isKeepDfromB = true;
	
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
  const Int_t nbincent = 28;  //bins in centrality
  const Int_t nbincent_0_10 = 4;  //bins in centrality between 0 and 10
  const Int_t nbincent_10_60 = 20;  //bins in centrality between 10 and 60
  const Int_t nbincent_60_100 = 4;  //bins in centrality between 60 and 100
  const Int_t nbinfake = 3;  //bins in fake
  const Int_t nbinpointingXY = 50;  //bins in cosPointingAngleXY
  const Int_t nbinnormDecayLXY = 20;  //bins in NormDecayLengthXY
  const Int_t nbinmult = 21;  //bins in NormDecayLengthXY

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
		
  const Int_t nbinpt = cutsLcTMVA->GetNPtBins(); // bins in pT
  iBin[ipT] = nbinpt;
  iBin[ipTpi] = nbinpt;
  iBin[ipTk] = nbinpt;
  Double_t *binLimpT = new Double_t[iBin[ipT]+1];
  Double_t *binLimpTpi = new Double_t[iBin[ipTpi]+1];
  Double_t *binLimpTk = new Double_t[iBin[ipTk]+1];
  // values for bin lower bounds
  Float_t* floatbinLimpT = cutsLcTMVA->GetPtBinLimits();
  for (Int_t ibin0 = 0 ; ibin0<iBin[ipT]+1; ibin0++){
    binLimpT[ibin0] = (Double_t)floatbinLimpT[ibin0];
    binLimpTpi[ibin0] = (Double_t)floatbinLimpT[ibin0];
    binLimpTk[ibin0] = (Double_t)floatbinLimpT[ibin0];
  }
  for(Int_t i=0; i<=nbinpt; i++) printf("binLimpT[%d]=%f\n",i,binLimpT[i]);  
	
  printf("pT: nbin (from cuts file) = %d\n",nbinpt);

  // defining now the binning for the other variables:
	
  AliLog::SetClassDebugLevel("AliCFManager",AliLog::kInfo);

  iBin[iy] = nbiny;
  iBin[icosThetaStar] = nbincosThetaStar;
  iBin[icT] = nbincT;
  iBin[idca] = nbindca;
  iBin[id0xd0] = nbind0xd0;
  iBin[ipointing] = nbinpointing;
  iBin[iphi] = nbinphi;
  iBin[izvtx] = nbinzvtx;
  iBin[icent] = nbincent;
  iBin[ifake] = nbinfake;
  iBin[ipointingXY] = nbinpointingXY;
  iBin[inormDecayLXY] = nbinnormDecayLXY;
  iBin[imult] = nbinmult;
	
  //arrays for lower bounds :
  Double_t *binLimy = new Double_t[iBin[iy]+1];
  Double_t *binLimcosThetaStar = new Double_t[iBin[icosThetaStar]+1];
  Double_t *binLimcT = new Double_t[iBin[icT]+1];
  Double_t *binLimdca = new Double_t[iBin[idca]+1];
  Double_t *binLimd0xd0 = new Double_t[iBin[id0xd0]+1];
  Double_t *binLimpointing = new Double_t[iBin[ipointing]+1];
  Double_t *binLimphi = new Double_t[iBin[iphi]+1];
  Double_t *binLimzvtx = new Double_t[iBin[izvtx]+1];
  Double_t *binLimcent = new Double_t[iBin[icent]+1];
  Double_t *binLimfake = new Double_t[iBin[ifake]+1];
  Double_t *binLimpointingXY = new Double_t[iBin[ipointingXY]+1];
  Double_t *binLimnormDecayLXY = new Double_t[iBin[inormDecayLXY]+1];
  Double_t *binLimmult = new Double_t[iBin[imult]+1];

  // muliplicity
  for (Int_t i = 0; i <= 8; i++) binLimmult[i] = 300 + i*100;
  for (Int_t i = 9; i <= 18; i++) binLimmult[i] = binLimmult[8] + (i-8)*100;
  for (Int_t i = 19; i <= 21; i++) binLimmult[i] = binLimmult[18] + (i-18)*500;
  if (isCent) {
    for (Int_t i = 0; i <= 18; i++) binLimmult[i] = 2000 + i*200;
    for (Int_t i = 19; i <= 21; i++) binLimmult[i] = binLimmult[18] + (i-18)*500;
  }
  
  for (Int_t i = 0; i < iBin[imult]+1; i++){ Printf("binLimmult[%d] = %f", i, binLimmult[i]);}
  // y
  for(Int_t i=0; i<=nbiny; i++) binLimy[i] = (Double_t)ymin  + (ymax-ymin)  /nbiny*(Double_t)i ;

  // cosThetaStar
  for(Int_t i=0; i<=nbincosThetaStar; i++) binLimcosThetaStar[i] = (Double_t)cosminTS  + (cosmaxTS-cosminTS)  /nbincosThetaStar*(Double_t)i ;
        
  // cT
  for(Int_t i=0; i<=nbincT; i++) binLimcT[i] = (Double_t)cTmin  + (cTmax-cTmin)  /nbincT*(Double_t)i ;

  // dca
  for(Int_t i=0; i<=nbindca; i++) binLimdca[i] = (Double_t)dcamin  + (dcamax-dcamin)  /nbindca*(Double_t)i ;

  // d0xd0
  for(Int_t i=0; i<=nbind0xd0; i++) binLimd0xd0[i] = (Double_t)d0xd0min  + (d0xd0max-d0xd0min)  /nbind0xd0*(Double_t)i ;

  // cosPointingAngle
  for(Int_t i=0; i<=nbinpointing; i++) binLimpointing[i] = (Double_t)cosmin  + (cosmax-cosmin)  /nbinpointing*(Double_t)i ;

  // Phi
  for(Int_t i=0; i<=nbinphi; i++) binLimphi[i] = (Double_t)phimin  + (phimax-phimin)  /nbinphi*(Double_t)i ;

  // z Primary Vertex
  for(Int_t i=0; i<=nbinzvtx; i++) {
    binLimzvtx[i] = (Double_t)zmin  + (zmax-zmin)  /nbinzvtx*(Double_t)i ;
  }

  // centrality
  for(Int_t i=0; i<=nbincent_0_10; i++) binLimcent[i] = (Double_t)centmin_0_10 + (centmax_0_10-centmin_0_10)/nbincent_0_10*(Double_t)i ; 
  if (binLimcent[nbincent_0_10] != centmin_10_60)  {
    Printf("Error: Calculated bin lim for cent - 1st range - differs from expected!\n");
  }
  for(Int_t i=0; i<=nbincent_10_60; i++) binLimcent[i+nbincent_0_10] = (Double_t)centmin_10_60 + (centmax_10_60-centmin_10_60)/nbincent_10_60*(Double_t)i ;
  if (binLimcent[nbincent_0_10+nbincent_10_60] != centmin_60_100)  {
    Printf("Error: Calculated bin lim for cent - 2st range - differs from expected!\n");
  }
  for(Int_t i=0; i<=nbincent_60_100; i++) binLimcent[i+nbincent_0_10+nbincent_10_60] = (Double_t)centmin_60_100 + (centmax_60_100-centmin_60_100)/nbincent_60_100*(Double_t)i ;

  // fake
  for(Int_t i=0; i<=nbinfake; i++) {
    binLimfake[i] = (Double_t)fakemin  + (fakemax-fakemin)/nbinfake * (Double_t)i;
  }

  // cosPointingAngleXY
  for(Int_t i=0; i<=nbinpointingXY; i++) binLimpointingXY[i] = (Double_t)cosminXY  + (cosmaxXY-cosminXY)  /nbinpointingXY*(Double_t)i ;

  // normDecayLXY
  for(Int_t i=0; i<=nbinnormDecayLXY; i++) binLimnormDecayLXY[i] = (Double_t)normDecLXYmin  + (normDecLXYmax-normDecLXYmin)  /nbinnormDecayLXY*(Double_t)i ;

  //one "container" for MC
  TString nameContainer="";
  if(!isKeepDfromB) {
    nameContainer = "CFHFccontainer0";
  }
  else  if(isKeepDfromBOnly){
    nameContainer = "CFHFccontainer0DfromB";
  }
  else  {
    nameContainer = "CFHFccontainer0allD";	  
  }
  nameContainer += suffix;
  //Setting up the container grid... 

  AliCFContainer* container;

  if (configuration == AliCFTaskVertexingHF::kSnail){
    container = new AliCFContainer(nameContainer,"container for tracks",nstep,nvarTot,iBin);
    //setting the bin limits
    printf("pt\n");
    container -> SetBinLimits(ipT, binLimpT);
    printf("y\n");
    container -> SetBinLimits(iy, binLimy);
    printf("cts\n");
    container -> SetBinLimits(icosThetaStar, binLimcosThetaStar);
    printf("ptPi\n");
    container -> SetBinLimits(ipTpi, binLimpTpi);
    printf("ptK\n");
    container -> SetBinLimits(ipTk, binLimpTk);
    printf("cT\n");
    container -> SetBinLimits(icT, binLimcT);
    printf("dca\n");
    container -> SetBinLimits(idca, binLimdca);
    printf("d0xd0\n");
    container -> SetBinLimits(id0xd0, binLimd0xd0);
    printf("pointing\n");
    container -> SetBinLimits(ipointing, binLimpointing);
    printf("phi\n");
    container -> SetBinLimits(iphi, binLimphi);
    printf("z\n");
    container -> SetBinLimits(izvtx, binLimzvtx);
    printf("cent\n");
    container -> SetBinLimits(icent, binLimcent);
    printf("fake\n");
    container -> SetBinLimits(ifake, binLimfake);
    printf("pointingXY\n");
    container -> SetBinLimits(ipointingXY, binLimpointingXY);
    printf("normDecayLXY\n");
    container -> SetBinLimits(inormDecayLXY, binLimnormDecayLXY);
    printf("multiplicity\n");
		
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

    container = new AliCFContainer(nameContainer,"container for tracks", nstep, nvar, iBinFast);
    printf("pt\n");
    container -> SetBinLimits(ipTFast, binLimpT);
    printf("y\n");
    container -> SetBinLimits(iyFast, binLimy);
    printf("ct\n");
    container -> SetBinLimits(icTFast, binLimcT);
    printf("phi\n");
    container -> SetBinLimits(iphiFast, binLimphi);
    printf("zvtx\n");
    container -> SetBinLimits(izvtxFast, binLimzvtx);
    printf("centrality\n");
    container -> SetBinLimits(icentFast, binLimcent);
    printf("fake\n");
    container -> SetBinLimits(ifakeFast, binLimfake);
    printf("multiplicity\n");

    container -> SetVarTitle(ipTFast, "pt");
    container -> SetVarTitle(iyFast, "y");
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

    container = new AliCFContainer(nameContainer,"container for tracks", nstep, nvar, iBinSuperFast);
    printf("pt\n");
    container -> SetBinLimits(ipTSuperFast, binLimpT);
    printf("y\n");
    container -> SetBinLimits(iySuperFast, binLimy);
    printf("centrality\n");
    container -> SetBinLimits(icentSuperFast, binLimcent);
    printf("multiplicity\n");

    container -> SetVarTitle(ipTSuperFast, "pt");
    container -> SetVarTitle(iySuperFast, "y");
    container -> SetVarTitle(icentSuperFast, "centrality");
    container -> SetVarTitle(imultSuperFast, "multiplicity");
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
  kineAccCuts->SetPtRange(ptmin, ptmax);
  kineAccCuts->SetEtaRange(etamin, etamax);

  // Rec-Level kinematic cuts
  AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts", "rec-level kine cuts");
	
  AliCFTrackQualityCuts *recQualityCuts = new AliCFTrackQualityCuts("recQualityCuts", "rec-level quality cuts");
	
  AliCFTrackIsPrimaryCuts *recIsPrimaryCuts = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts", "rec-level isPrimary cuts");
	
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
  AliCFTaskVertexingHF *task = new AliCFTaskVertexingHF("AliCFTaskVertexingHF", cutsLcTMVA);
  task->SetConfiguration(configuration);
  task->SetFillFromGenerated(kFALSE);
  task->SetCFManager(man); //here is set the CF manager
  task->SetDecayChannel(22);
  task->SetUseCascadeTaskForLctoV0bachelor(kTRUE);
  task->SetUseAdditionalCuts(kFALSE);
  task->SetUseCutsForTMVA(kFALSE);
  task->SetUseFlatPtWeight(useFlatPtWeight); 
  task->SetUseWeight(useWeight);
  task->SetUseZWeight(useZWeight);
  task->SetCutOnMomConservation(cutOnMomConservation);
  Printf("CutOn mom = %f", cutOnMomConservation);
  task->SetSign(isSign);
  task->SetCentralitySelection(kFALSE);
  task->SetFakeSelection(0);
  task->SetRejectCandidateIfNotFromQuark(rejectIfNotFromQuark); // put to false if you want to keep HIJING D0!!
  task->SetUseMCVertex(kFALSE); // put to true if you want to do studies on pp
  //task->SetPtWeightsFromDataPbPb276overLHC12a17a();

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

  if (useWeight){
    TH1F *fHistoPtWeight_010_1 = (TH1F*)fileCuts->Get("fHistoPtWeight_010_1");
    TH1F *fHistoPtWeight_010_2 = (TH1F*)fileCuts->Get("fHistoPtWeight_010_2");
    TH1F *fHistoPtWeight_3050_1 = (TH1F*)fileCuts->Get("fHistoPtWeight_3050_1");
    TH1F *fHistoPtWeight_3050_2 = (TH1F*)fileCuts->Get("fHistoPtWeight_3050_2");
    if(!fHistoPtWeight_010_1) {
      Printf("FATAL: Histogram for pt weights not found");
      return 0x0;
    }
    if(!fHistoPtWeight_010_2) {
      Printf("FATAL: Histogram for pt weights not found");
      return 0x0;
    }
    if(!fHistoPtWeight_3050_1) {
      Printf("FATAL: Histogram for pt weights not found");
      return 0x0;
    }
    if(!fHistoPtWeight_3050_2) {
      Printf("FATAL: Histogram for pt weights not found");
      return 0x0;
    }
    TH1F* weightHisto = (TH1F*)fileCuts->Get(weightHistoName);
    task->SetWeightHistogram(weightHisto);
  }

 

  Printf("***************** TASK SETTINGS *****************");       
  Printf("UseCascadeTask = %d", (Int_t)task->GetUseCascadeTaskForLctoV0bachelor());
  Printf("decay channel = %d", (Int_t)task->GetDecayChannel());
  Printf("configuration = %d", task->GetConfiguration());
  Printf("cut on momentum conservation = %f", task->GetCutOnMomConservation());
  Printf("Running for central (for binning in multiplicity) = %d", (Int_t)isCent);
  Printf("FillFromGenerated = %d", (Int_t)task->GetFillFromGenerated());
  Printf("RejectCandidateIfNotFromQuark selection = %d", (Int_t)task->GetRejectCandidateIfNotFromQuark());
  Printf("Dselection = %d", (Int_t)task->GetDselection());
  Printf("UseWeight = %d", (Int_t)task->GetUseWeight());
  if (task->GetUseWeight()) {
    if(funcWeight) Printf("User-defined Weight function");
    else Printf("FONLL will be used for the weights");
  }
  Printf("Use Nch weight = %d", (Int_t)task->GetUseNchWeight());
  Printf("Sign = %d", (Int_t)task->GetSign());
  Printf("Centrality selection = %d", (Int_t)task->GetCentralitySelection());
  Printf("Fake selection = %d", (Int_t)task->GetFakeSelection());
  Printf("UseMCVertex selection = %d", (Int_t)task->GetUseMCVertex());
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
    nameCorr="CFHFcorr0";
  }
  else  if(isKeepDfromBOnly){
    nameCorr= "CFHFcorr0KeepDfromBOnly";
  }
  else  {
    nameCorr="CFHFcorr0allD";

  }
  nameCorr += suffix;

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
  TString output1name="", output2name="", output3name="",output4name="", output5name="";
  output2name=nameContainer;
  output3name=nameCorr;
  output4name= "Cuts";
  output5name= "coutProfDst";
  if(!isKeepDfromB) {
    outputfile += ":PWG3_D2H_CFtaskLctoK0SpCascade";
    output1name="CFHFchist0";
    output3name+="_cOnly";
    output4name+="_cOnly";
    output5name+="_cOnly";
  }
  else  if(isKeepDfromBOnly){
    outputfile += ":PWG3_D2H_CFtaskLctoK0SpCascadeKeepDfromBOnly";
    output1name="CFHFchist0DfromB";
    output3name+="_bOnly";
    output4name+="_bOnly";
    output5name+="_bOnly";
  }
  else{
    outputfile += ":PWG3_D2H_CFtaskLctoK0SpCascadeKeepDfromB";
    output1name="CFHFchist0allD";
    output3name+="_all";
    output4name+="_all";
    output5name+="_all";
  }

  outputfile += suffix;
  output1name += suffix;
  output4name += suffix;
  output5name += suffix;

  //now comes user's output objects :
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(output1name, TH1I::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  // output Correction Framework Container (for acceptance & efficiency calculations)
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(output2name, AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  // Unfolding - correlation matrix
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(output3name, THnSparseD::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  // cuts
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
