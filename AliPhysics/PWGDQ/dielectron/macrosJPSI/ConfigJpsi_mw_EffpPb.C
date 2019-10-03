
void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);


TString names=("nocut;TPCrefit;ESDcuts+SPDany;ESDcuts+SPDfirst;ESDcuts+SPDany+TPCpid;ESDcuts+SPDfirst+TPCpid");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDieEff=arrNames->GetEntriesFast();

AliDielectron* ConfigJpsi_mw_EffpPb(Int_t cutDefinition)
{
  //
  // Setup the instance of AliDielectron
  //
  
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

  //tracking cuts
  AliDielectronTrackCuts *trackCuts=new AliDielectronTrackCuts("ITSandgeneral_trackCuts","ITSandgeneral_trackCuts");
  if(cutDefinition==1) {
    trackCuts->SetRequireTPCRefit(kTRUE);
  }
  if(cutDefinition>1) {
    trackCuts->SetRequireITSRefit(kTRUE);
    trackCuts->SetRequireTPCRefit(kTRUE);
  }
  if(cutDefinition==2 || cutDefinition==4) {
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }
  if(cutDefinition==3 || cutDefinition==5) {
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  }
  die->GetTrackFilter().AddCuts(trackCuts);
  
  if(cutDefinition>1){
    AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("trackkineandTPCQ","trackkine_and_TPC");
    
    varCuts->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
    varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    varCuts->AddCut(AliDielectronVarManager::kNclsTPC,80.,160.);
    varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    die->GetTrackFilter().AddCuts(varCuts);
  }
  if (cutDefinition>=4){
    //ESD pid cuts
    AliDielectronPID *pid=new AliDielectronPID("MC_Prod_Data","MC to reproduce data");
    //proton cut to reproduce data parametrisation
    Double_t resolution=0.055;
    Double_t nSigma=3.;
    TF1 *ffPro=new TF1(Form("fBethe%d",AliPID::kProton), Form("(%f*%f+(AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])-AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])))/%f", nSigma,resolution, AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kElectron), resolution), 0.05,200.);
    ffPro->SetParameters(0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00);
    
    TF1 *ffPio=new TF1(Form("fBethe%d",AliPID::kPion), Form("(%f*%f+(AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])-AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])))/%f", nSigma,resolution, AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kElectron), resolution), 0.05,200.);
    ffPio->SetParameters(0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00);
    
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,3);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,ffPro,10,0,3);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,ffPio,10);
//     pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3,3,0.,0.,kTRUE);
    
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
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",400,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","nClsoverfindablecluster","Number of found Clusters TPC over findably ;TPC number cluster over findable;#tracks",160,0.0,1.1,AliDielectronVarManager::kNFclsTPCrFrac);
  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        400,1e-2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  histos->UserHistogram("Track","TPCnSigma_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        400,1e-2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  

  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        125,0.,125*.04,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMass_OpeningAngle","Opening angle:Inv.Mass;Inv. Mass [GeV];angle",
                        100,0.,4.,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  
  die->SetHistogramManager(histos);
}

void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  //pair variables **********************************************************
  // j/psi pt ------------------
  cf->AddVariable(AliDielectronVarManager::kPt, "0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.5, 3.8, 4.2, 4.6, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0");
  // j/psi y -------------------
  cf->AddVariable(AliDielectronVarManager::kY, "-5,-1,-0.9,-0.8,-0.7,-0.5,-0.3,0.3,0.5,0.7,0.8,0.9,1.0,5");
  //j/psi mass
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  // pair type -----------------
  cf->AddVariable(AliDielectronVarManager::kPairType,      3,  0.0, 3.0);  //NOTE
  // cos theta* Collins-Soper
  // cf->AddVariable(AliDielectronVarManager::kThetaCS,      20, -1.0, 1.0);
  // cos theta* Helicity
  //cf->AddVariable(AliDielectronVarManager::kThetaHE,      20, -1.0, 1.0);
  //leg variables **********************************************************
  // leg pseudo-rapidity --------------------------
  cf->AddVariable(AliDielectronVarManager::kEta, "-1.0, -0.88, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.88, 1.0", kTRUE);
  // leg TPC n-sigma electron ---------------------
  //  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle, 12, -3.0,   3.0, kTRUE);
  // leg pt ---------------------------------------
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 100.0",kTRUE);
  // leg Ncls TPC ---------------------------------
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 70, 75, 80, 85, 90, 100, 120, 160",kTRUE);
  // leg ITS first cluster point
  cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,-0.5,5.5,kTRUE);
  //leg momentum
  cf->AddVariable(AliDielectronVarManager::kP,"0.0, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2,1.5,2.0,3.0,4.0,5.0, 100.0",kTRUE);
  // -------------------------------------------------------------------------------
  //event variables*********************************************************
  cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd10,101,-0.5,100.5);
  cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd10Corr,101,-0.5,100.5);

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
