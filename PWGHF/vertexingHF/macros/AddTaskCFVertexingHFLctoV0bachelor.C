//DEFINITION OF A FEW CONSTANTS
const Double_t ptmin =    0.0;
const Double_t ptmax = 9999.0;
const Double_t etamin = -0.9;
const Double_t etamax =  0.9;
const Double_t ymin  = -1.2 ;
const Double_t ymax  =  1.2 ;

const Double_t zmin   = -15.;
const Double_t zmax   =  15.;
const Float_t centmin = 0.;
const Float_t centmax = 100.;
const Float_t fakemin = -0.5;
const Float_t fakemax =  2.5;
const Float_t multmin =   0.;
const Float_t multmax = 102.;

const Double_t cosPAV0min = 0.985;
const Double_t cosPAV0max = 1.005;
const Float_t onFlymin = -0.5;
const Float_t onFlymax =  1.5;
const Double_t pBachmin   =  0.0;
const Double_t pBachmax   = 30.0;
const Double_t ptV0min    =  0.0;
const Double_t ptV0max    = 30.0;
const Double_t yV0min     =-1.2;
const Double_t yV0max     = 1.2;
const Double_t dcaV0min   = 0.0; // nSigma
const Double_t dcaV0max   = 1.5; // nSigma
const Double_t cTV0min    = 0.0; // micron
const Double_t cTV0max    = 300; // micron
const Double_t cTmin      = 0.0; // micron
const Double_t cTmax      = 300; // micron
const Float_t cosPAmin    =-1.05;
const Float_t cosPAmax    = 1.05;


//----------------------------------------------------

AliCFTaskVertexingHF *AddTaskCFVertexingHFLctoV0bachelor(const char* cutFile = "./LctoV0bachelorCuts.root",
							 Bool_t rejectIfNotFromQuark=kTRUE,
							 //Bool_t isKeepDfromB = kTRUE, Bool_t isKeepDfromBOnly = kFALSE, // all in
							 Bool_t isKeepDfromB = kFALSE, Bool_t isKeepDfromBOnly = kFALSE, // prompt
							 //Bool_t isKeepDfromB = kTRUE, Bool_t isKeepDfromBOnly = kTRUE, // no-prompt
							 Int_t configuration = AliCFTaskVertexingHF::kCheetah,
							 Int_t pdgCode = 4122, Char_t isSign = 2, Char_t lcToV0bachelorDecayMode = 0,
							 TString usercomment = "username",
							 Bool_t useWeight=kFALSE)
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

  const Double_t phimin = 0.;
  const Double_t phimax = 2.*TMath::Pi();
  const Double_t phiV0min = 0.;
  const Double_t phiV0max = 2.*TMath::Pi();

  const Int_t nbinZvtx       =  30; //bins in centrality (total number)
  const Int_t nbincent       =18+1; //bins in centrality (total number)
  const Int_t nbinfake       =   3; //bins in fake
  const Int_t nbinmult       =  48; //bins in multiplicity (total number)

  const Int_t nbinpt         =  11; //bins in pt from 0,1,2,3,4,5,6,8,12,17,25,35 GeV
  const Int_t nbiny          =  24; //bins in y Lc
  const Int_t nbinphi        =  18; //bins in phi Lc
  const Int_t nbinonFly      =   2; //bins in onFlyStatus x V0

  const Int_t nbinpBach      = 300; //bins in pt from 0 to 30 GeV
  const Int_t nbinptV0       = 300; //bins in pt from 0 to 30 GeV
  const Int_t nbinyV0        =  24; //bins in y V0
  const Int_t nbinphiV0      =  18; //bins in phi V0
  const Int_t nbindcaV0      = 150; //bins in dcaV0
  const Int_t nbininvMassV0  =  60; //bins in invMassV0
  const Int_t nbinpointingV0 =  42; //bins in cosPointingAngleV0
  const Int_t nbinpointing   =  42; //bins in cosPointingAngle

  //const Int_t nbincTV0      =  15; //bins in cTV0
  //const Int_t nbincT        =  15; //bins in cT

  //the sensitive variables, their indices

  // variables' indices
  const UInt_t ipT        = 0;
  const UInt_t iy         = 1;
  const UInt_t iphi       = 2;
  const UInt_t ionFly     = 3;
  const UInt_t iZvtx      = 4;
  const UInt_t icent      = 5;
  const UInt_t ifake      = 6;
  const UInt_t imult      = 7;

  const UInt_t ipbach    =  8;
  const UInt_t ipTV0     =  9;
  const UInt_t iyV0      = 10;
  const UInt_t iphiV0    = 11;
  const UInt_t iinvMassV0= 12;
  const UInt_t idcaV0    = 13;
  const UInt_t icosPAxV0 = 14;
  const UInt_t icosPA    = 15;
  //const UInt_t icTv0     = 16;
  //const UInt_t icT       = 17;

  //Setting the bins: pt, ptPi, and ptK are considered seprately because for them you can either define the binning by hand, or using the cuts file

  //arrays for the number of bins in each dimension

  //if ( configuration ==AliCFTaskVertexingHF::kSnail)
  const Int_t nvarTot   = 16 ; //number of variables on the grid
  //if ( configuration ==AliCFTaskVertexingHF::kCheetah)
  //const Int_t nvarTot =  8 ; //number of variables on the grid

  Int_t iBin[nvarTot];

  //OPTION 1: defining the pt, ptPi, ptK bins by hand...		
  iBin[ipT]=nbinpt;
  iBin[iy]=nbiny;
  iBin[iphi]=nbinphi;
  iBin[ionFly]=nbinonFly;
  iBin[iZvtx]=nbinZvtx;
  iBin[icent]=nbincent;
  iBin[ifake]=nbinfake;
  iBin[imult]=nbinmult;

  iBin[ipbach]=nbinpBach;
  iBin[ipTV0]=nbinptV0;
  iBin[iyV0]=nbinyV0;
  iBin[iphiV0]=nbinphiV0;
  iBin[iinvMassV0]=nbininvMassV0;
  iBin[idcaV0]=nbindcaV0;
  iBin[icosPAxV0]=nbinpointingV0;
  iBin[icosPA]=nbinpointing;
  //iBin[icTv0]=nbincTV0;
  //iBin[icT]=nbincT;

  // values for bin lower bounds

  // pt
  Double_t *binLimpT=new Double_t[iBin[ipT]+1];
  binLimpT[ 0]= 0.;
  binLimpT[ 1]= 1.;
  binLimpT[ 2]= 2.;
  binLimpT[ 3]= 3.;
  binLimpT[ 4]= 4.;
  binLimpT[ 5]= 5.;
  binLimpT[ 6]= 6.;
  binLimpT[ 7]= 8.;
  binLimpT[ 8]=12.;
  binLimpT[ 9]=17.;
  binLimpT[10]=25.;
  binLimpT[11]=35.;

  // y
  Double_t *binLimy=new Double_t[iBin[iy]+1];
  for(Int_t i=0; i<=iBin[iy]; i++) binLimy[i]=(Double_t)ymin + (ymax-ymin)/iBin[iy]*(Double_t)i ;

  // phi
  Double_t *binLimphi=new Double_t[iBin[iphi]+1];
  for(Int_t i=0; i<=iBin[iphi]; i++) binLimphi[i]=(Double_t)phimin + (phimax-phimin)/iBin[iphi]*(Double_t)i ;

  // onTheFlyV0
  Double_t *binLimonFlyV0=new Double_t[iBin[ionFly]+1];
  for(Int_t i=0; i<=iBin[ionFly]; i++) binLimonFlyV0[i]=(Double_t)onFlymin + (onFlymax-onFlymin)/iBin[ionFly]*(Double_t)i ;

  // z Primary Vertex
  Double_t *binLimzvtx=new Double_t[iBin[iZvtx]+1];
  for(Int_t i=0; i<=nbinZvtx; i++) binLimzvtx[i]=(Double_t)zmin  + (zmax-zmin)/iBin[iZvtx]*(Double_t)i ;

  // centrality
  Double_t *binLimcent=new Double_t[iBin[icent]+1];
  binLimcent[0]=centmin;
  for(Int_t i=1; i<=iBin[icent]; i++) binLimcent[i]=centmin + (centmax-centmin)/iBin[icent]*(Double_t)i;

  // fake
  Double_t *binLimfake=new Double_t[iBin[ifake]+1];
  for(Int_t i=0; i<=iBin[ifake]; i++) binLimfake[i]=(Double_t)fakemin  + (fakemax-fakemin)/iBin[ifake] * (Double_t)i;

  // multiplicity
  Double_t *binLimmult=new Double_t[iBin[imult]+1];
  for(Int_t i=0; i<=iBin[imult]; i++) binLimmult[i]=(Double_t)multmin + (multmax-multmin)/iBin[imult]*(Double_t)i ; 


  // pBach
  Double_t *binLimpbach=new Double_t[iBin[ipbach]+1];
  for(Int_t i=0; i<=iBin[ipbach]; i++) binLimpbach[i]=(Double_t)pBachmin + (pBachmax-pBachmin)/iBin[ipbach]*(Double_t)i ; 

  // ptV0
  Double_t *binLimpTV0=new Double_t[iBin[ipTV0]+1];
  for(Int_t i=0; i<=iBin[ipTV0]; i++) binLimpTV0[i]=(Double_t)ptV0min + (ptV0max-ptV0min)/iBin[ipTV0]*(Double_t)i ; 

  // yV0
  Double_t *binLimyV0=new Double_t[iBin[iyV0]+1];
  for(Int_t i=0; i<=iBin[iyV0]; i++) binLimyV0[i]=(Double_t)yV0min + (yV0max-yV0min)/iBin[iyV0]*(Double_t)i ; 

  // phiV0
  Double_t *binLimphiV0=new Double_t[iBin[iphiV0]+1];
  for(Int_t i=0; i<=iBin[iphiV0]; i++) binLimphiV0[i]=(Double_t)phiV0min + (phiV0max-phiV0min)/iBin[iphiV0]*(Double_t)i ; 

  // invMassV0
  Double_t *binLimInvMassV0=new Double_t[iBin[iinvMassV0]+1];
  for(Int_t i=0; i<=iBin[iinvMassV0]; i++) binLimInvMassV0[i]=(Double_t)massV0min + (massV0max-massV0min)/iBin[iinvMassV0]*(Double_t)i ;

  // dcaV0
  Double_t *binLimdcaV0=new Double_t[iBin[idcaV0]+1];
  for(Int_t i=0; i<=iBin[idcaV0]; i++) binLimdcaV0[i]=(Double_t)dcaV0min + (dcaV0max-dcaV0min)/iBin[idcaV0]*(Double_t)i ; 

  // cosPointingAngleV0
  Double_t *binLimcosPAV0=new Double_t[iBin[icosPAxV0]+1];
  for(Int_t i=0; i<=iBin[icosPAxV0]; i++) binLimcosPAV0[i]=(Double_t)cosPAV0min + (cosPAV0max-cosPAV0min)/iBin[icosPAxV0]*(Double_t)i ;

  // cosPointingAngle
  Double_t *binLimcosPA=new Double_t[iBin[icosPA]+1];
  for(Int_t i=0; i<=iBin[icosPA]; i++) binLimcosPA[i]=(Double_t)cosPAmin + (cosPAmax-cosPAmin)/iBin[icosPA]*(Double_t)i ;

  /*
  // cTV0
  Double_t *binLimcTV0=new Double_t[iBin[icTv0]+1];
  for(Int_t i=0; i<=iBin[icTv0]; i++) binLimcTV0[i]=(Double_t)cTV0min + (cTV0max-cTV0min)/iBin[icTv0]*(Double_t)i ; 

  // cT
  Double_t *binLimcT=new Double_t[iBin[icT]+1];
  for(Int_t i=0; i<=iBin[icT]; i++) binLimcT[i]=(Double_t)cTmin + (cTmax-cTmin)/iBin[icT]*(Double_t)i ;
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
  container -> SetBinLimits(ipT,binLimpT);
  container -> SetBinLimits(iy,binLimy);
  container -> SetBinLimits(iphi,binLimphi);
  container -> SetBinLimits(ionFly,binLimonFlyV0);
  container -> SetBinLimits(iZvtx,binLimzvtx);
  container -> SetBinLimits(icent,binLimcent);
  container -> SetBinLimits(ifake,binLimfake);
  container -> SetBinLimits(imult,binLimmult);

  container -> SetVarTitle(ipT,"p_{T}(#Lambda_{c}) [GeV/c]");
  container -> SetVarTitle(iy,"y(#Lambda_{c})");
  container -> SetVarTitle(iphi,"#phi(#Lambda_{c}) [rad]");
  container -> SetVarTitle(ionFly,"onTheFlyStatusV0");
  container -> SetVarTitle(iZvtx,"z_{vtx} [cm]");
  container -> SetVarTitle(icent,"centrality");
  container -> SetVarTitle(ifake,"fake");
  container -> SetVarTitle(imult,"multiplicity");

  if (configuration == AliCFTaskVertexingHF::kSnail) {
    container -> SetBinLimits(ipbach,binLimpbach);
    container -> SetBinLimits(ipTV0,binLimpTV0);
    container -> SetBinLimits(iyV0,binLimyV0);
    container -> SetBinLimits(iphiV0,binLimphiV0);
    container -> SetBinLimits(iinvMassV0,binLimInvMassV0);
    container -> SetBinLimits(idcaV0,binLimdcaV0);
    container -> SetBinLimits(icosPAxV0,binLimcosPAV0);
    container -> SetBinLimits(icosPA,binLimcosPA);
    //container -> SetBinLimits(,binLimcTV0);
    //container -> SetBinLimits(,binLimcT);

    container -> SetVarTitle(ipbach,"p(bachelor) [GeV/c]");
    container -> SetVarTitle(ipTV0,"p_{T}(V0) [GeV/c]");
    container -> SetVarTitle(iyV0,"y(V0)");
    container -> SetVarTitle(iphiV0,"#varphi(V0) [rad]");
    container -> SetVarTitle(iinvMassV0,"m_{inv}(#pi^{+},#pi^{-}) [GeV/c^{2}]");
    container -> SetVarTitle(idcaV0,"DCA(V0) [n#sigma]");
    container -> SetVarTitle(icosPAxV0,"cosine pointing angle(V0)");
    container -> SetVarTitle(icosPA,"cosine pointing angle (#Lambda_{c})");
    //container -> SetVarTitle(,"c#tau -V0-");
    //container -> SetVarTitle(,"c#tau");
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
  task->SetUseWeight(useWeight);
  task->SetSign(isSign);
  task->SetCentralitySelection(kFALSE);
  task->SetFakeSelection(0);
  task->SetRejectCandidateIfNotFromQuark(rejectIfNotFromQuark); // put to false if you want to keep HIJING D0!!
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

  thnDim[0] = iBin[ipT];
  thnDim[2] = iBin[ipT];
  thnDim[1] = iBin[iy];
  thnDim[3] = iBin[iy];

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
