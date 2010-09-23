

void SetupTrackCuts();
void SetupPairCuts();
void InitHistograms();
void InitCF();

AliESDtrackCuts *SetupESDtrackCuts();

TString names=("basicQ+SPDfirst+pt>.6+PID;basicQ+SPDany+pt>.6+PID");

TObjArray *arrNames=names.Tokenize(";");

const Int_t nDie=arrNames->GetEntries();

AliDielectron *fDiele=0x0;
Int_t          fCutDefinition=0;
Bool_t         fIsAOD=kFALSE;

AliDielectron* ConfigJpsi2ee(Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setup the instance of AliDielectron
  //
  
  fCutDefinition=cutDefinition;
  fIsAOD=isAOD;
  
  // create the actual framework object
  TString name=Form("%02d",fCutDefinition);
  if (fCutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(fCutDefinition)->GetName();
  }
  fDiele = new AliDielectron(Form("%s",name.Data()),
                             Form("Track cuts: %s",name.Data()));

  // cut setup
  SetupTrackCuts();
  SetupPairCuts();
  
  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // fDielelectron framework histograms will be filled
  //
  InitHistograms();

  // the last definition uses no cuts and only the QA histograms should be filled!
//   if (fCutDefinition<nDie-1)
  InitCF();

  return fDiele;
}

//______________________________________________________________________________________
void SetupTrackCuts()
{
  //
  // Setup the track cuts
  //
  
  //ESD quality cuts DielectronTrackCuts
  if (!fIsAOD) {
    fDiele->GetTrackFilter().AddCuts(SetupESDtrackCuts());
  } else {
    AliDielectronTrackCuts *trackCuts=new AliDielectronTrackCuts("trackCuts","trackCuts");
    if (fCutDefinition==0)
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    else if (fCutDefinition==1)
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetRequireITSRefit(kTRUE);
    fDiele->GetTrackFilter().AddCuts(trackCuts);
  }

  //Pt cut ----------------------------------------------------------
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");
  pt->AddCut(AliDielectronVarManager::kPt,0.6,1e30);

  //AOD additions since there are no AliESDtrackCuts -----------------
  //
  if (fIsAOD){
    // TPC #clusteres cut
    pt->AddCut(AliDielectronVarManager::kNclsTPC,90.,160.);
    pt->AddCut(AliDielectronVarManager::kEta,-0.8,0.8);
    //TODO: DCA cuts to be investigated!!!
//       pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
//       pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }
  fDiele->GetTrackFilter().AddCuts(pt);
    
  // PID cuts --------------------------------------------------------
  AliDielectronPID *pid = new AliDielectronPID("PID10","TPC nSigma 2.5<e<4. + |Pi|>3 + |P|>3");
  pid->SetDefaults(2);
  fDiele->GetTrackFilter().AddCuts(pid);
}

//______________________________________________________________________________________
void SetupPairCuts()
{
  //
  // Setup the pair cuts
  //
  
  //Invariant mass and rapidity selection
  AliDielectronVarCuts *pairCut=new AliDielectronVarCuts("2<M<4+|Y|<.8","2<M<4 + |Y|<.8");
  pairCut->AddCut(AliDielectronVarManager::kM,2.,4.);
  pairCut->AddCut(AliDielectronVarManager::kY,-0.8,0.8);
  fDiele->GetPairFilter().AddCuts(pairCut);
}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts()
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
  
  if (fCutDefinition==0)
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  else if (fCutDefinition==1)
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  
  return esdTrackCuts;
}


//______________________________________________________________________________________
void InitHistograms()
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(fDiele->GetName(),fDiele->GetTitle());
  
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
  
  fDiele->SetHistogramManager(histos);
}


void InitCF()
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(fDiele->GetName(),fDiele->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 1.4, 2.8, 4.2, 9.9, 100.0");
  
  cf->AddVariable(AliDielectronVarManager::kY,40,-2,2);
  cf->AddVariable(AliDielectronVarManager::kM,50,1.98,1.98+50*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 1.2, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 100, 120, 160",kTRUE);
  
  fDiele->SetCFManagerPair(cf);
  
}

