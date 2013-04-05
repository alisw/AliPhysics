//DEFINITION OF A FEW CONSTANTS
const Double_t ptmin = 0.0;
const Double_t ptmax = 12.0;
const Double_t ymin  = -1.2 ;
const Double_t ymax  =  1.2 ;
const Double_t cosPAV0min = 0.99;
const Double_t cosPAV0max = 1.02;
const Int_t onFlymin = 0;
const Int_t onFlymax = 2;
const Float_t centmin =   0.;
//const Float_t centmin_0_10   =   0.;
//const Float_t centmax_0_10   =  10.;
//const Float_t centmin_10_60  =  10.;
//const Float_t centmax_10_60  =  60.;
//const Float_t centmin_60_100 =  60.;
//const Float_t centmax_60_100 = 100.;
const Float_t centmax = 100.;
const Int_t fakemin = 0;
const Int_t fakemax = 3;
const Float_t multmin =   0;
//const Float_t multmin_0_20 =     0;
//const Float_t multmax_0_20 =    20;
//const Float_t multmin_20_50 =   20;
//const Float_t multmax_20_50 =   50;
//const Float_t multmin_50_102 =  50;
//const Float_t multmax_50_102 = 102;
const Float_t multmax = 102;

const Double_t ptBachmin  = 0.0;
const Double_t ptBachmax  = 30.0;
const Double_t ptV0posmin = 0.0;
const Double_t ptV0posmax = 30.0;
const Double_t ptV0negmin = 0.0;
const Double_t ptV0negmax = 30.0;
const Double_t dcaV0min   = 0.0; // nSigma
const Double_t dcaV0max   = 1.5; // nSigma
const Double_t cTV0min    = 0.0; // micron
const Double_t cTV0max    = 300; // micron
const Double_t cTmin      = 0.0; // micron
const Double_t cTmax      = 300; // micron
const Float_t cosPAmin    =-1.02;
const Float_t cosPAmax    = 1.02;

const Double_t etamin = -0.9;
const Double_t etamax =  0.9;
//const Double_t zmin = -15.;
//const Double_t zmax =  15.;


//----------------------------------------------------

AliCFTaskVertexingHF *AddTaskCFVertexingHFLctoV0bachelor(const char* cutFile = "./LctoV0bachelorCuts.root",
							 Int_t configuration = AliCFTaskVertexingHF::kCheetah, Bool_t isKeepDfromB = kTRUE,
							 Bool_t isKeepDfromBOnly = kFALSE, Int_t pdgCode = 4122, Char_t isSign = 2,
							 Char_t lcToV0bachelorDecayMode = 0,
							 TString usercomment = "username")
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

  // isSign = 0 --> Lc+ only
  // isSign = 1 --> Lc- only
  // isSign = 2 --> Lc+ and Lc-

  TString expected;
  if (isSign == 0 && pdgCode < 0){
    AliError(Form("Error setting PDG code (%d) and sign (0 --> Lc+ only): they are not compatible, returning",pdgCode));
    return 0x0;
  }
  else if (isSign == 1 && pdgCode > 0){
    AliError(Form("Error setting PDG code (%d) and sign (1 --> Lc- only): they are not compatible, returning",pdgCode));
    return 0x0;
  }
  else if (isSign > 2 || isSign < 0){
    AliError(Form("Sign not valid (%d, possible values are 0, 1, 2), returning",isSign));
    return 0x0;
  }

  TFile* fileCuts = TFile::Open(cutFile);
  AliRDHFCuts *cutsLctoV0 = (AliRDHFCutsLctoV0*)fileCuts->Get("LctoV0AnalysisCuts");

  // check that the fKeepD0fromB flag is set to true when the fKeepD0fromBOnly flag is true
  //  for now the binning is the same than for all D's
  if (isKeepDfromBOnly) isKeepDfromB = true;

  Double_t massV0min = 0.47;
  Double_t massV0max = 1.14;
  if (lcToV0bachelorDecayMode==0) {
    massV0min = 0.47 ;
    massV0max = 0.53 ;
  } else if (lcToV0bachelorDecayMode==1) {
    massV0min = 1.09;
    massV0max = 1.14;
  }

  const Double_t phimin = 0.0;
  const Double_t phimax = 2.*TMath::Pi();

  const Int_t nbinpt          =  8; //bins in pt from 0 to 12 GeV
  const Int_t nbiny           = 24; //bins in y
  const Int_t nbininvMassV0   = 60; //bins in invMassV0
  const Int_t nbinpointingV0  = 12; //bins in cosPointingAngleV0
  const Int_t nbinonFly       =  2; //bins in onFlyStatus x V0

  const Int_t nbincent        = 18; //bins in centrality (total number)
  //const Int_t nbincent_0_10   =  4; //bins in centrality between 0 and 10
  //const Int_t nbincent_10_60  = 10; //bins in centrality between 10 and 60
  //const Int_t nbincent_60_100 =  4; //bins in centrality between 60 and 100

  const Int_t nbinfake        =  3; //bins in fake

  const Int_t nbinmult        = 48; //bins in multiplicity (total number)
  //const Int_t nbinmult_0_20   = 20; //bins in multiplicity between 0 and 20
  //const Int_t nbinmult_20_50  = 15; //bins in multiplicity between 20 and 50
  //const Int_t nbinmult_50_102 = 13; //bins in multiplicity between 50 and 102


  const Int_t nbinptBach      = 300; //bins in pt from 0 to 30 GeV
  const Int_t nbinptV0pos     = 300; //bins in pt from 0 to 30 GeV
  const Int_t nbinptV0neg     = 300; //bins in pt from 0 to 30 GeV
  const Int_t nbinphi         = 18; //bins in Phi
  const Int_t nbindcaV0       = 15; //bins in dcaV0
  const Int_t nbincTV0        = 15; //bins in cTV0
  const Int_t nbincT          = 15; //bins in cT
  const Int_t nbinpointing    = 51; //bins in cosPointingAngle

  //the sensitive variables, their indices

  // variables' indices
  const UInt_t ipT        = 0;
  const UInt_t iy         = 1;
  const UInt_t iphi       = 2;
  const UInt_t icosPAxV0  = 3;
  const UInt_t ionFly     = 4;
  const UInt_t imult      = 5;
  const UInt_t ifake      = 6;
  const UInt_t icent      = 7;

  const UInt_t ipTbach   =  8;
  const UInt_t ipTposV0  =  9;
  const UInt_t ipTnegV0  = 10;
  const UInt_t iinvMassV0= 11;
  const UInt_t idcaV0    = 12;
  const UInt_t icTv0     = 13;
  const UInt_t icT       = 14;
  const UInt_t icosPA    = 15;

  //Setting the bins: pt, ptPi, and ptK are considered seprately because for them you can either define the binning by hand, or using the cuts file

  //arrays for the number of bins in each dimension

  //if ( configuration ==AliCFTaskVertexingHF::kSnail)
  const Int_t nvarTot   = 16 ; //number of variables on the grid:pt, y, cosThetaStar, pTpi, pTk, cT, dca, d0pi, d0K, d0xd0, cosPointingAngle, phi, z, centrality, fake, cosPointingAngleXY, normDecayLengthXY, multiplicity
  //if ( configuration ==AliCFTaskVertexingHF::kCheetah)
  //const Int_t nvarTot   =  8 ; //number of variables on the grid:pt, y, cosThetaStar, pTpi, pTk, cT, dca, d0pi, d0K, d0xd0, cosPointingAngle, phi, z, centrality, fake, cosPointingAngleXY, normDecayLengthXY, multiplicity

  Int_t iBin[nvarTot];

  //OPTION 1: defining the pt, ptPi, ptK bins by hand...		
  iBin[ipT]=nbinpt;
  iBin[iy]=nbiny;
  iBin[iphi]=nbinphi;
  iBin[icosPAxV0]=nbinpointingV0;
  iBin[ionFly]=nbinonFly;
  iBin[imult]=nbinmult;
  iBin[ifake]=nbinfake;
  iBin[icent]=nbincent;

  iBin[ipTbach]=nbinptBach;
  iBin[ipTposV0]=nbinptV0pos;
  iBin[ipTnegV0]=nbinptV0neg;
  iBin[iinvMassV0]=nbininvMassV0;
  iBin[idcaV0]=nbindcaV0;
  iBin[icTv0]=nbincTV0;
  iBin[icT]=nbincT;
  iBin[icosPA]=nbinpointing;

  // values for bin lower bounds

  // pt
  Double_t *binLimpT=new Double_t[iBin[0]+1];
  //for(Int_t ii=0; ii<=iBin[0]; ii++) binLimpT[ii]=(Double_t)ptmin + (ptmax-ptmin)/iBin[0]*(Double_t)ii ; 
  for(Int_t ii=0; ii<=iBin[0]-2; ii++) binLimpT[ii]=(Double_t)ptmin + (Double_t)ii;
  binLimpT[iBin[0]-1]=8.;
  binLimpT[iBin[0]]=12.;

  // y
  Double_t *binLimy=new Double_t[iBin[1]+1];
  for(Int_t i=0; i<=iBin[1]; i++) binLimy[i]=(Double_t)ymin + (ymax-ymin)/iBin[1]*(Double_t)i ;

  // phi
  Double_t *binLimphi=new Double_t[iBin[2]+1];
  for(Int_t i=0; i<=iBin[2]; i++) binLimphi[i]=(Double_t)phimin  + (phimax-phimin)/iBin[2]*(Double_t)i ;

  // cosPointingAngleV0
  Double_t *binLimcosPAV0=new Double_t[iBin[3]+1];
  for(Int_t i=0; i<=iBin[3]; i++) binLimcosPAV0[i]=(Double_t)cosPAV0min + (cosPAV0max-cosPAV0min)/iBin[3]*(Double_t)i ;

  // onTheFlyV0
  Double_t *binLimonFlyV0=new Double_t[iBin[4]+1];
  for(Int_t i=0; i<=iBin[4]; i++) binLimonFlyV0[i]=(Double_t)onFlymin + (onFlymax-onFlymin)/iBin[4]*(Double_t)i ;

  // centrality
  Double_t *binLimcent=new Double_t[iBin[5]+1];
  for(Int_t i=0; i<=iBin[5]; i++) binLimcent[i]=(Double_t)centmin + (centmax-centmin)/iBin[5]*(Double_t)i ; 
  /*
  Double_t *binLimcent_0_10=new Double_t[iBin[icent]+1];
  for(Int_t i=0; i<=nbincent_0_10; i++) binLimcent[i]=(Double_t)centmin_0_10 + (centmax_0_10-centmin_0_10)/nbincent_0_10*(Double_t)i ; 
  if (binLimcent[nbincent_0_10] != centmin_10_60)  {
    Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for cent - 1st range - differs from expected!\n");
  }

  Double_t *binLimcent_10_60=new Double_t[iBin[icent]+1];
  for(Int_t i=0; i<=nbincent_10_60; i++) binLimcent[i+nbincent_0_10]=(Double_t)centmin_10_60 + (centmax_10_60-centmin_10_60)/nbincent_10_60*(Double_t)i ;
  if (binLimcent[nbincent_0_10+nbincent_10_60] != centmin_60_100)  {
    Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for cent - 2st range - differs from expected!\n");
  }

  Double_t *binLimcent_60_100=new Double_t[iBin[icent]+1];
  for(Int_t i=0; i<=nbincent_60_100; i++) binLimcent[i+nbincent_10_60]=(Double_t)centmin_60_100 + (centmax_60_100-centmin_60_100)/nbincent_60_100*(Double_t)i ;
  */
        
  // fake
  Double_t *binLimfake=new Double_t[iBin[6]+1];
  for(Int_t i=0; i<=iBin[6]; i++) binLimfake[i]=(Double_t)fakemin  + (fakemax-fakemin)/iBin[6] * (Double_t)i;

  // multiplicity
  Double_t *binLimmult=new Double_t[iBin[7]+1];
  for(Int_t i=0; i<=iBin[7]; i++) binLimmult[i]=(Double_t)multmin + (multmax-multmin)/iBin[7]*(Double_t)i ; 
  /*
  for(Int_t i=0; i<=nbinmult_0_20; i++) binLimmult[i]=(Double_t)multmin_0_20 + (multmax_0_20-multmin_0_20)/nbinmult_0_20*(Double_t)i ; 
  if (binLimmult[nbinmult_0_20] != multmin_20_50)  {
    Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for mult - 1st range - differs from expected!\n");
  }
  for(Int_t i=0; i<=nbinmult_20_50; i++) binLimmult[i+nbinmult_0_20]=(Double_t)multmin_20_50 + (multmax_20_50-multmin_20_50)/nbinmult_20_50*(Double_t)i ; 
  if (binLimmult[nbinmult_0_20+nbinmult_20_50] != multmin_50_102)  {
    Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for mult - 2nd range - differs from expected!\n");
  }
  for(Int_t i=0; i<=nbinmult_50_102; i++) binLimmult[i+nbinmult_0_20+nbinmult_20_50]=(Double_t)multmin_50_102 + (multmax_50_102-multmin_50_102)/nbinmult_50_102*(Double_t)i ; 
  */


  // ptBach
  Double_t *binLimpTbach=new Double_t[iBin[8]+1];
  for(Int_t i=0; i<=iBin[8]; i++) binLimpTbach[i]=(Double_t)ptBachmin + (ptBachmax-ptBachmin)/iBin[8]*(Double_t)i ; 

  // ptV0pos
  Double_t *binLimpTV0pos=new Double_t[iBin[9]+1];
  for(Int_t i=0; i<=iBin[9]; i++) binLimpTV0pos[i]=(Double_t)ptV0posmin + (ptV0posmax-ptV0posmin)/iBin[9]*(Double_t)i ; 

  // ptV0neg
  Double_t *binLimpTV0neg=new Double_t[iBin[10]+1];
  for(Int_t i=0; i<=iBin[10]; i++) binLimpTV0neg[i]=(Double_t)ptV0negmin + (ptV0negmax-ptV0negmin)/iBin[10]*(Double_t)i ; 

  // invMassV0
  Double_t *binLimInvMassV0=new Double_t[iBin[11]+1];
  for(Int_t i=0; i<=iBin[11]; i++) binLimInvMassV0[i]=(Double_t)massV0min + (massV0max-massV0min)/iBin[11]*(Double_t)i ;

  // dcaV0
  Double_t *binLimdcaV0=new Double_t[iBin[12]+1];
  for(Int_t i=0; i<=iBin[12]; i++) binLimdcaV0[i]=(Double_t)dcaV0min + (dcaV0max-dcaV0min)/iBin[12]*(Double_t)i ; 

  // cTV0
  Double_t *binLimcTV0=new Double_t[iBin[13]+1];
  for(Int_t i=0; i<=iBin[13]; i++) binLimcTV0[i]=(Double_t)cTV0min + (cTV0max-cTV0min)/iBin[13]*(Double_t)i ; 

  // cT
  Double_t *binLimcT=new Double_t[iBin[14]+1];
  for(Int_t i=0; i<=iBin[14]; i++) binLimcT[i]=(Double_t)cTmin + (cTmax-cTmin)/iBin[14]*(Double_t)i ;

  // cosPointingAngle
  Double_t *binLimcosPA=new Double_t[iBin[15]+1];
  for(Int_t i=0; i<=iBin[15]; i++) binLimcosPA[i]=(Double_t)cosPAmin + (cosPAmax-cosPAmin)/iBin[15]*(Double_t)i ;



  // z Primary Vertex
  //for(Int_t i=0; i<=nbinzvtx; i++) binLimzvtx[i]=(Double_t)zmin  + (zmax-zmin)  /nbinzvtx*(Double_t)i ;


  //OPTION 2: ...or from the cuts file
  /*
  const Int_t nbinpt = cutsLctoV0->GetNPtBins(); // bins in pT
  iBin[ipT]=nbinpt;
  iBin[ipTpi]=nbinpt;
  iBin[ipTk]=nbinpt;
  Double_t *binLimpT=new Double_t[iBin[ipT]+1];
  Double_t *binLimpTpi=new Double_t[iBin[ipTpi]+1];
  Double_t *binLimpTk=new Double_t[iBin[ipTk]+1];
  // values for bin lower bounds
  Float_t* floatbinLimpT = cutsLctoV0->GetPtBinLimits();
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
  if (binLimcent[nbincent_0_10] != centmin_10_60)  {
    Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for cent - 1st range - differs from expected!\n");
  }
  for(Int_t i=0; i<=nbincent_10_60; i++) binLimcent[i+nbincent_0_10]=(Double_t)centmin_10_60 + (centmax_10_60-centmin_10_60)/nbincent_10_60*(Double_t)i ;
  if (binLimcent[nbincent_0_10+nbincent_10_60] != centmin_60_100)  {
    Error("AliCFHeavyFlavourTaskMultiVarMultiStep","Calculated bin lim for cent - 2st range - differs from expected!\n");
  }
  for(Int_t i=0; i<=nbincent_60_100; i++) binLimcent[i+nbincent_10_60]=(Double_t)centmin_60_100 + (centmax_60_100-centmin_60_100)/nbincent_60_100*(Double_t)i ;
        
  // fake
  for(Int_t i=0; i<=nbinfake; i++) {
    binLimfake[i]=(Double_t)fakemin  + (fakemax-fakemin)/nbinfake * (Double_t)i;
  }

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
  */



  //one "container" for MC
  TString nameContainer="";
  if (!isKeepDfromB) {
    nameContainer="CFHFccontainer0_CommonFramework_"+usercomment;
  }
  else if (isKeepDfromBOnly) {
    nameContainer="CFHFccontainer0LcfromB_CommonFramework_"+usercomment;
  }
  else  {
    nameContainer="CFHFccontainer0allLc_CommonFramework_"+usercomment;	  
  }

  //Setting up the container grid... 

  //CONTAINER DEFINITION
  Info("AliCFTaskVertexingHF","SETUP CONTAINER");
  UInt_t nstep = 10; //number of selection steps: MC with limited acceptance, MC, Acceptance, Vertex, Refit, Reco (no cuts), RecoAcceptance, RecoITSClusters (RecoAcceptance included), RecoPPR (RecoAcceptance+RecoITSCluster included), RecoPID 

  AliCFContainer* container;
  if (configuration == AliCFTaskVertexingHF::kSnail) {
    container = new AliCFContainer(nameContainer,"container for tracks",nstep,nvarTot,iBin);
  }
  else if (configuration == AliCFTaskVertexingHF::kCheetah) {
    container = new AliCFContainer(nameContainer,"container for tracks",nstep,8,iBin);
  }

  //setting the bin limits
  container -> SetBinLimits(0,binLimpT);
  container -> SetBinLimits(1,binLimy);
  container -> SetBinLimits(2,binLimphi);
  container -> SetBinLimits(3,binLimcosPAV0);
  container -> SetBinLimits(4,binLimonFlyV0);
  container -> SetBinLimits(5,binLimcent);
  container -> SetBinLimits(6,binLimfake);
  container -> SetBinLimits(7,binLimmult);

  container -> SetVarTitle(0,"pt");
  container -> SetVarTitle(1,"y");
  container -> SetVarTitle(2,"phi");
  container -> SetVarTitle(3,"cosPA -V0-");
  container -> SetVarTitle(4,"onFlyV0");
  container -> SetVarTitle(5,"centrality");
  container -> SetVarTitle(6,"fake");
  container -> SetVarTitle(7,"multiplicity");

  if (configuration == AliCFTaskVertexingHF::kSnail) {
    container -> SetBinLimits(8,binLimpTbach);
    container -> SetBinLimits(9,binLimpTV0pos);
    container -> SetBinLimits(10,binLimpTV0neg);
    container -> SetBinLimits(11,binLimInvMassV0);
    container -> SetBinLimits(12,binLimdcaV0);
    container -> SetBinLimits(13,binLimcTV0);
    container -> SetBinLimits(14,binLimcT);
    container -> SetBinLimits(15,binLimcosPA);

    container -> SetVarTitle(8,"ptBachelor");
    container -> SetVarTitle(9,"ptV0pos");
    container -> SetVarTitle(10,"ptV0neg");
    container -> SetVarTitle(11,"mV0");
    container -> SetVarTitle(12,"DCA -V0-");
    container -> SetVarTitle(13,"c#tau -V0-");
    container -> SetVarTitle(14,"c#tau");
    container -> SetVarTitle(15,"cosPA");
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
  if (isSign != 2) {
    useAbsolute = kFALSE;
  }
  mcGenCuts->SetRequirePdgCode(pdgCode, useAbsolute);  // kTRUE set in order to include Lc-
  mcGenCuts->SetAODMC(1); //special flag for reading MC in AOD tree (important)

  // Acceptance cuts:
  AliCFAcceptanceCuts* accCuts = new AliCFAcceptanceCuts("accCuts", "Acceptance cuts");
  AliCFTrackKineCuts * kineAccCuts = new AliCFTrackKineCuts("kineAccCuts","Kine-Acceptance cuts");
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
  AliCFTaskVertexingHF *task = new AliCFTaskVertexingHF("AliCFTaskVertexingHF",cutsLctoV0);
  task->SetConfiguration(configuration);
  task->SetFillFromGenerated(kFALSE);
  task->SetCFManager(man); //here is set the CF manager
  task->SetDecayChannel(22);//kLctoV0bachelor
  switch (lcToV0bachelorDecayMode) {
  case 0:
    task->SetCountLctoK0Sp();
    break;
  case 1:
    task->SetCountLctoLambdapi();
    break;
  }
  task->SetUseWeight(kFALSE);
  task->SetSign(isSign);
  task->SetCentralitySelection(kFALSE);
  task->SetFakeSelection(0);
  task->SetRejectCandidateIfNotFromQuark(kTRUE); // put to false if you want to keep HIJING D0!!
  task->SetUseMCVertex(kFALSE); // put to true if you want to do studies on pp

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
      task->GetWeightFunction(funcWeight)->Print();
    }
  }

  Printf("***************** CONTAINER SETTINGS *****************");	
  Printf("decay channel = %d",(Int_t)task->GetDecayChannel());
  Printf("FillFromGenerated = %d",(Int_t)task->GetFillFromGenerated());
  Printf("Dselection = %d",(Int_t)task->GetDselection());
  Printf("UseWeight = %d",(Int_t)task->GetUseWeight());
  if (task->GetUseWeight()) {
    Printf("User-defined Weight function:");
    task->GetWeightFunction(funcWeight)->Print();
  }
  else {
    Printf("FONLL will be used for the weights");
  }
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
  if (!isKeepDfromB) {
    nameCorr="CFHFcorr0_CommonFramework_"+usercomment;
  }
  else if (isKeepDfromBOnly) {
    nameCorr= "CFHFcorr0KeepDfromBOnly_CommonFramework_"+usercomment;
  }
  else {
    nameCorr="CFHFcorr0allLc_CommonFramework_"+usercomment;
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
  if (!isKeepDfromB) {
    outputfile += ":PWG3_D2H_CFtaskLctoK0Sp_CommonFramework_"+usercomment;
    output1name="CFHFchist0_CommonFramework_"+usercomment;
    output4name= "Cuts_CommonFramework_"+usercomment;
  }
  else  if (isKeepDfromBOnly) {
    outputfile += ":PWG3_D2H_CFtaskLctoK0SpKeepDfromBOnly_CommonFramework_"+usercomment;
    output1name="CFHFchist0DfromB_CommonFramework_"+usercomment;
    output4name= "Cuts_CommonFramework_DfromB_"+usercomment;
  }
  else {
    outputfile += ":PWG3_D2H_CFtaskLctoK0SpKeepDfromB_CommonFramework_"+usercomment;
    output1name="CFHFchist0allLc_CommonFramework_"+usercomment;
    output4name= "Cuts_CommonFramework_allLc_"+usercomment;
  }

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
