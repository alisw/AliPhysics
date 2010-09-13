
void InitHistograms(AliDielectron *die);
void InitCF(AliDielectron* die);

void SetupTrackCuts(AliDielectron *die);
void SetupPairCuts(AliDielectron *die);

AliESDtrackCuts *SetupESDtrackCuts();

AliDielectron* ConfigJpsi2eeFilter()
{
  //
  // Setup the instance of AliDielectron
  //
  
  // create the actual framework object
  TString name="trackQ+Pt>0.5+60<dEdx<100";
  AliDielectron *die =
    new AliDielectron(Form("%s",name.Data()),
                      Form("Track cuts: %s",name.Data()));
  
  // cut setup
  SetupTrackCuts(die);
  SetupPairCuts(die);
  
  //
  // QA histogram setup
  //
  InitHistograms(die);
  
  return die;
}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die)
{
  //
  // Setup the track cuts
  //
  
  //ESD quality cuts
  die->GetTrackFilter().AddCuts(SetupESDtrackCuts());
  
  //Pt cut
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("Pt>.5+60<dEdx<100","Pt>.6 && 60<dEdx<100");
  pt->AddCut(AliDielectronVarManager::kPt,.6,1e30);
  pt->AddCut(AliDielectronVarManager::kTPCsignal,60.,100.);
  
  die->GetTrackFilter().AddCuts(pt);
}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die)
{
  //
  // Setup the pair cuts
  //
  
  
  //Invarian mass selection
  AliDielectronVarCuts *invMassCut=new AliDielectronVarCuts("InvMass","2<M<4");
  invMassCut->AddCut(AliDielectronVarManager::kM,2.,1e30);
//   invMassCut->AddCut(AliDielectronVarManager::kPairType,1.);
  die->GetPairFilter().AddCuts(invMassCut);

}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts()
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
  
  esdTrackCuts->SetMaxDCAToVertexZ(3.0);
  esdTrackCuts->SetMaxDCAToVertexXY(1.0); 
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  
  esdTrackCuts->SetMinNClustersTPC(90);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);

  return esdTrackCuts;
}


//______________________________________________________________________________________
void InitHistograms(AliDielectron *die)
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
