
void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);

TString names=("nocut;TPCrefit;ESDcuts;ESDcuts+SPDany+TPCpid;ESDcuts+SPDfirst+TPCpid");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();

AliDielectron* ConfigJpsi2ee(Int_t cutDefinition, Bool_t isAOD=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  if (isAOD) return 0x0;
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die =
    new AliDielectron(Form("%s",name.Data()),
                      Form("Track cuts: %s",name.Data()));
  
  // cut setup
  SetupTrackCuts(die,cutDefinition);
  SetupPairCuts(die,cutDefinition);
  
  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  //  dielectron framework histograms will be filled
  //
  InitHistograms(die,cutDefinition);
  InitCF(die,cutDefinition);
  
  return die;
}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  //ESD quality cuts
  if (cutDefinition>0){
    die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));
  }

  if (cutDefinition>=3){
    //ESD pid cuts
    AliDielectronPID *pid=new AliDielectronPID("MC_Prod_Data","MC to reproduce data");
    //proton cut to reproduce data parametrisation
    Double_t resolution=0.058;
    Double_t nSigma=3.;
    TF1 *ff=new TF1(Form("fBethe%d",AliPID::kProton), Form("(%f*%f+(AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])-AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])))/%f", nSigma,resolution, AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kElectron), resolution), 0.05,200.);
    ff->SetParameters(0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,3);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,ff,10,0,3);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3,3,kTRUE);
    die->GetTrackFilter().AddCuts(pid);
  } 
}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //

  //reject conversions
  AliDielectronVarCuts *openingAngleCut=new AliDielectronVarCuts("OpeningAngle","Opening angle>0.35rad");
  openingAngleCut->AddCut(AliDielectronVarManager::kOpeningAngle,.035,4.);
  if(cutDefinition>1)
    die->GetPairFilter().AddCuts(openingAngleCut);  
}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition)
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
  if(cutDefinition==1) {
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
  }
  if(cutDefinition>1) {
    //    esdTrackCuts->SetEtaRange(-0.9,0.9);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetMaxDCAToVertexZ(3.0);
    esdTrackCuts->SetMaxDCAToVertexXY(1);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(80);
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition==3) {
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }
  if(cutDefinition==4) {
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  }

  return esdTrackCuts;
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram classes
  AliDielectronHistos *histos=
    new AliDielectronHistos(die->GetName(),
                            die->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  
  //Event class
  histos->AddClass("Event");
  
  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
  
  //Pair classes
  // to fill also mixed event histograms loop until 10
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
  }
  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",
                        1,0.,1.,AliDielectronVarManager::kNevents);
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Eta","Eta;Eta;#tracks",400,-2.0,-2.0,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","Phi;Phi;#tracks",650,0,6.5,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Theta","Theta;Theta;#tracks",350,0,3.5,AliDielectronVarManager::kTheta);
  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        400,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        400,0.2,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions;#tracks",
                        400,0.2,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","TPCnSigma_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        400,1e-2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        500,0.,4.,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-2.,2.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMass_OpeningAngle","Opening angle:Inv.Mass;Inv. Mass [GeV];angle",
                        100,0.,4.,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","ThetaHE","Theta (helicity system);#theta [rad.]",
                        500,0.,1.0,AliDielectronVarManager::kThetaHE);
  histos->UserHistogram("Pair","PhiHE","Phi (helicity system);#varphi [rad.]",
                        500,-6.5,6.5,AliDielectronVarManager::kPhiHE);
  histos->UserHistogram("Pair","ThetaCS","Theta (Collins-Soper system);#theta [rad.]",
                        500,0.0,1.0,AliDielectronVarManager::kThetaCS);
  histos->UserHistogram("Pair","PhiCS","Phi (Collins-Soper system);#varphi [rad.]",
                        500,-6.5,6.5,AliDielectronVarManager::kPhiCS);
  
  die->SetHistogramManager(histos);
}

//______________________________________________________________________________________
void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  //pair variables **********************************************************
  // j/psi pt ------------------
  TVectorD *binLimPt=new TVectorD(14);
  (*binLimPt)[0]=0.0; (*binLimPt)[1]=0.4; (*binLimPt)[2]=0.8; (*binLimPt)[3]=1.4;
  (*binLimPt)[4]=2.0; (*binLimPt)[5]=2.8; (*binLimPt)[6]=3.5; (*binLimPt)[7]=4.2;
  (*binLimPt)[8]=5.0; (*binLimPt)[9]=6.0; (*binLimPt)[10]=7.0; (*binLimPt)[11]=8.0;
  (*binLimPt)[12]=9.0; (*binLimPt)[13]=10.0;
  cf->AddVariable(AliDielectronVarManager::kPt,          binLimPt);
  // j/psi y -------------------
  TVectorD *binLimY=new TVectorD(13);
  (*binLimY)[0]=-1.0; (*binLimY)[1]=-0.88; (*binLimY)[2]=-0.8;
  (*binLimY)[3]=-0.6; (*binLimY)[4]=-0.4; (*binLimY)[5]=-0.2; (*binLimY)[6]=0.0;
  (*binLimY)[7]=0.2; (*binLimY)[8]=0.4; (*binLimY)[9]=0.6; (*binLimY)[10]=0.8;
  (*binLimY)[11]=0.88; (*binLimY)[12]=1.0;
  cf->AddVariable(AliDielectronVarManager::kY,            binLimY);
  // pair type -----------------
  cf->AddVariable(AliDielectronVarManager::kPairType,      3,  0.0, 3.0);  
  //leg variables **********************************************************
  // leg pseudo-rapidity --------------------------
  cf->AddVariable(AliDielectronVarManager::kEta,          binLimY, kTRUE);
  // leg TPC n-sigma electron ---------------------
  //  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle, 12, -3.0,   3.0, kTRUE);
  // leg pt ---------------------------------------
  TVectorD *binLimPtLeg=new TVectorD(7);
  (*binLimPtLeg)[0]=0.0; (*binLimPtLeg)[1]=0.8; (*binLimPtLeg)[2]=0.9;
  (*binLimPtLeg)[3]=1.0; (*binLimPtLeg)[4]=1.1; (*binLimPtLeg)[5]=1.2;
  (*binLimPtLeg)[6]=20.0;
  cf->AddVariable(AliDielectronVarManager::kPt,            binLimPtLeg, kTRUE);
  // leg Ncls TPC ---------------------------------
  TVectorD *binLimNclsTPC=new TVectorD(7);
  (*binLimNclsTPC)[0]=0.0;   (*binLimNclsTPC)[1]=80.0;  (*binLimNclsTPC)[2]=90.0;
  (*binLimNclsTPC)[3]=100.0; (*binLimNclsTPC)[4]=110.0; (*binLimNclsTPC)[5]=120.0;
  (*binLimNclsTPC)[6]=160.0;
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,       binLimNclsTPC, kTRUE);
  // -------------------------------------------------------------------------------

  if(cutDefinition==0) {
    cf->SetStepForMCtruth();
    cf->SetStepForAfterAllCuts(kFALSE);
    cf->SetStepsForSignal(kFALSE);
  }
  if(cutDefinition>0){
    //    cf->SetStepForNoCutsMCmotherPid();
    cf->SetStepForAfterAllCuts();
    //    cf->SetStepsForEachCut();
    cf->SetStepsForSignal();
    //    cf->SetStepsForBackground();
  }
  die->SetCFManagerPair(cf);
  
}
