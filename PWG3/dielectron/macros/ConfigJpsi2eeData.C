
void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);

TString names=("basicQ+SPDfirst+pt>.6+PID");

TObjArray *arrNames=names.Tokenize(";");

const Int_t nDie=arrNames->GetEntries();

AliDielectron* ConfigJpsi2ee(Int_t cutDefinition)
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
  // dielectron framework histograms will be filled
  //
  InitHistograms(die,cutDefinition);

  // the last definition uses no cuts and only the QA histograms should be filled!
//   if (cutDefinition<nDie-1)
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
  die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));

  if(cutDefinition==0) {
    //Pt cut ----------------------------------------------------------
    AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");
    pt->AddCut(AliDielectronVarManager::kPt,0.6,1e30);
    die->GetTrackFilter().AddCuts(pt);
    
    // PID cuts --------------------------------------------------------
    AliDielectronPID *pid = new AliDielectronPID("PID10","TPC nSigma |e|<3 + |Pi|>3 + |P|>3 + TOF nSigma |e|<3");
    pid->SetDefaults(10);
    die->GetTrackFilter().AddCuts(pid);
  }
}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //
  
  //Invariant mass selection
  AliDielectronVarCuts *invMassCut=new AliDielectronVarCuts("2<M<4+|Y|<.8","2<M<4 + |Y|<.8");
  invMassCut->AddCut(AliDielectronVarManager::kM,2.,4.);
  invMassCut->AddCut(AliDielectronVarManager::kY,-0.8,0.8);
  die->GetPairFilter().AddCuts(invMassCut);
}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition)
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
  
  // basic track quality cuts  (basicQ)
  esdTrackCuts->SetMaxDCAToVertexZ(3.0);
  esdTrackCuts->SetMaxDCAToVertexXY(1.0);
  
  esdTrackCuts->SetEtaRange( -0.8 , 0.8 );
  
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  
  esdTrackCuts->SetMinNClustersTPC(90);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  
  return esdTrackCuts;
}


//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  
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
  
  //legs from pair
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
  }
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",400,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                        200,-1,1,200,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        400,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        400,0.2,20.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
      
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        201,-.01,4.01,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  
  die->SetHistogramManager(histos);
}


void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  //pair variables
  TVectorD *binLimPt=new TVectorD(6);
  (*binLimPt)[0]=0.0; (*binLimPt)[1]=0.8; (*binLimPt)[2]=1.4; (*binLimPt)[3]=2.8; (*binLimPt)[4]=4.2; (*binLimPt)[5]=9.9;
  cf->AddVariable(AliDielectronVarManager::kPt,binLimPt);
  cf->AddVariable(AliDielectronVarManager::kY,40,-2,2);
  cf->AddVariable(AliDielectronVarManager::kM,150,0.,150*.027); //27Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,2,0.8,1.2,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,3,90,120,kTRUE);
  
  die->SetCFManagerPair(cf);
  
}

