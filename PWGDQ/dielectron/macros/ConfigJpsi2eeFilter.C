void SetupTrackCutsDieleFilter(AliDielectron *diele, Bool_t isAOD);
void SetupPairCutsDieleFilter(AliDielectron *diele, Bool_t isAOD);
void SetupEventCutsDieleFilter(AliDielectron *diele, Int_t cutDefinition);

void InitHistogramsDieleFilter(AliDielectron *diele);

AliESDtrackCuts *SetupESDtrackCutsDieleFilter();


AliDielectron* ConfigJpsi2eeFilter(Bool_t isAOD=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //

  // create the actual framework object
  TString name="trackQ+Pt>0.6+60<dEdx<100";
  AliDielectron *diele = new AliDielectron(Form("%s",name.Data()),
                             Form("Track cuts: %s",name.Data()));
  
  // cut setup
  SetupEventCutsDieleFilter(diele);
  
  SetupTrackCutsDieleFilter(diele, isAOD);
  SetupPairCutsDieleFilter(diele, isAOD);
  
  //
  // QA histogram setup
  //
  InitHistogramsDieleFilter(diele, isAOD);
  
  return diele;
}

//______________________________________________________________________________________
void SetupEventCutsDieleFilter(AliDielectron *diele)
{
  //
  // Setup the event cuts
  //
  AliDielectronVarCuts *vtxZ = new AliDielectronVarCuts("vtxZ","Vertex z cut");
  vtxZ->AddCut(AliDielectronVarManager::kZvPrim,-15.,15.);
  diele->GetEventFilter().AddCuts(vtxZ);
}

//______________________________________________________________________________________
void SetupTrackCutsDieleFilter(AliDielectron *diele, Bool_t isAOD)
{
  //
  // Setup the track cuts
  //
  
  //ESD quality cuts DielectronTrackCuts
  if (!isAOD) {
    diele->GetTrackFilter().AddCuts(SetupESDtrackCutsDieleFilter());
  } else {
    AliDielectronTrackCuts *trackCuts=new AliDielectronTrackCuts("trackCuts","trackCuts");
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetRequireITSRefit(kTRUE);
    diele->GetTrackFilter().AddCuts(trackCuts);
  }
  
  //Pt cut
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("Pt>.5+60<dEdx<100","Pt>.6 && 60<dEdx<100");
  // pt > 0.7GeV
  pt->AddCut(AliDielectronVarManager::kPt,.7,1e30);
  pt->AddCut(AliDielectronVarManager::kTPCsignal,60.,100.);

  if (isAOD){
      // TPC #clusteres cut
    pt->AddCut(AliDielectronVarManager::kNclsTPC,90.,160.);
//     pt->AddCut(AliDielectronVarManager::kEta,-0.88,0.88);
    //TODO: DCA cuts to be investigated!!!
//     pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
//     pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }
  diele->GetTrackFilter().AddCuts(pt);

  // PID cuts ---------------------------------------------------
  AliDielectronPID *pid = new AliDielectronPID("PID","TPC nSigma |e|<3. + |Pi|>3 + |P|>3");
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3);
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.,3.,0,0,kTRUE);
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.,3.,0,0,kTRUE);
  
  diele->GetTrackFilter().AddCuts(pid);
}

//______________________________________________________________________________________
void SetupPairCutsDieleFilter(AliDielectron *diele, Bool_t isAOD)
{
  //
  // Setup the pair cuts
  //
  //Invarian mass selection
  AliDielectronVarCuts *invMassCut=new AliDielectronVarCuts("InvMass","2<M<4");
  // Minv > 1.8
  invMassCut->AddCut(AliDielectronVarManager::kM,1.8,1e30);
//invMassCut->AddCut(AliDielectronVarManager::kPairType,1.);
  // ptJpsi > 1GeV
  invMassCut->AddCut(AliDielectronVarManager::kPt,1.,1e30);
  diele->GetPairFilter().AddCuts(invMassCut);

}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCutsDieleFilter()
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
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
 
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);

  return esdTrackCuts;
}


//______________________________________________________________________________________
void InitHistogramsDieleFilter(AliDielectron *diele, Bool_t isAOD)
{
  //
  // Initialise the histograms
  //
  
//Setup histogram classes
  AliDielectronHistos *histos=
    new AliDielectronHistos(diele->GetName(),
                            diele->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");

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
  histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",300,-15.,15.,AliDielectronVarManager::kZvPrim);
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",400,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                        100,-1,1,144,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  
  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        200,0.2,20.,100,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        201,-.01,4.01,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  
  diele->SetHistogramManager(histos);
}
