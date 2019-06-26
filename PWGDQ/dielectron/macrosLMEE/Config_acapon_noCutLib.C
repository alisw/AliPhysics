// Simple config file not requiring cutLlibrary.
// Used for testing. Would not recommend using as base for a new config.

void InitHistograms(AliDielectron* die, Int_t cutDefinition);
void SetupCuts(AliDielectron* die, Int_t cutDefinition);
AliESDtrackCuts* SetupESDtrackCuts(Int_t cutDefinition);
AliDielectronPID* SetPIDcuts(Int_t cutDefinition);
const AliDielectronEventCuts* GetEventCuts();
TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max);
TVectorD* GetVector(Int_t var);
enum {kMee = 0, kMee500, kPtee, kP2D, kRuns, kPhiV, kOpAng, kOpAng2, kEta2D, kEta3D, kSigmaEle, kSigmaOther, kTPCdEdx, kCent, kPhi2D};


AliDielectron* Config_acapon(TString name, Int_t cutDefinition = 0)
{
  // Setup the instance of AliDielectron
  AliDielectron* die = new AliDielectron(Form("%s", name.Data()), Form("Track cuts: %s", name.Data()));
  
  AliDielectronMixingHandler* mix = new AliDielectronMixingHandler;
  mix->SetMixType(AliDielectronMixingHandler::kAll);
  mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
  mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
  mix->SetDepth(30);
  die->SetMixingHandler(mix);

	// set track cuts
  SetupCuts(die,cutDefinition);

  // histogram setup
  InitHistograms(die,cutDefinition);

  die->SetNoPairing(kFALSE);
  
  return die;
}

//______________________________________________________________________________________
void SetupCuts(AliDielectron* die, Int_t cutDefinition)
{

  // Default is kTRUE...but will mess up results
  die->SetUseKF(kFALSE);
           
   AliDielectronCutGroup* allCuts  = new AliDielectronCutGroup("allCuts", "allCuts", AliDielectronCutGroup::kCompOR);

  // AND cut group to select low mass pairs with large opening angle
  AliDielectronCutGroup* convRejCut = new AliDielectronCutGroup("convRejCut", "convRejCut", AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts* convMassCut = new AliDielectronVarCuts("convMassCut", "convMassCut");
  AliDielectronVarCuts* convPhiVCut = new AliDielectronVarCuts("convPhiVCut", "convPhiVCut");
  convMassCut->AddCut(AliDielectronVarManager::kM, 0.00, 0.1);
  convPhiVCut->AddCut(AliDielectronVarManager::kPhivPair, 0., 2.);
  convRejCut->AddCut(convMassCut);
  convRejCut->AddCut(convPhiVCut);

  // Mass cut to include any pairs with mass greater than 0.1 GeV
  AliDielectronVarCuts* pairMassCut = new AliDielectronVarCuts("pairMassCut", "pairMassCut");
  pairMassCut->AddCut(AliDielectronVarManager::kM, 0.1, 5.0);

  allCuts->AddCut(convRejCut);
  allCuts->AddCut(pairMassCut);

  die->GetPairFilter().AddCuts(allCuts);
    
  die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));
  die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));

}
//______________________________________________________________________________________
//-----------------------------------pid------------------------------------------------

AliDielectronPID* SetPIDcuts(Int_t cutDefinition){
  
  AliDielectronPID* pid = new AliDielectronPID();

  if(cutDefinition == 0){
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE , AliDielectronPID::kRequire    , AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE, AliDielectronPID::kIfAvailable, AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS, AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE, AliDielectronPID::kRequire    , AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE, AliDielectronPID::kRequire    , AliDielectronVarManager::kPt);
  }	
 return pid;
}

//______________________________________________________________________________________
//-----------------------------------track cuts-----------------------------------------
AliESDtrackCuts* SetupESDtrackCuts(Int_t cutDefinition){

  AliESDtrackCuts* fesdTrackCuts = new AliESDtrackCuts();

  fesdTrackCuts->SetPtRange( 0.2 , 100. );
  fesdTrackCuts->SetEtaRange( -0.8 , 0.8 );
  fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fesdTrackCuts->SetDCAToVertex2D(kFALSE);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.);
  fesdTrackCuts->SetMaxDCAToVertexXY(1.);
 
  fesdTrackCuts->SetRequireTPCRefit(kTRUE);
  fesdTrackCuts->SetRequireITSRefit(kTRUE);

  if(cutDefinition == 0){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  return fesdTrackCuts;
}
//______________________________________________________________________________________
//-------------------------------Histogram definition-----------------------------------
void InitHistograms(AliDielectron* die, Int_t cutDefinition)
{
  
  // Setup histogram classes
  AliDielectronHistos* histos= new AliDielectronHistos(die->GetName(), die->GetTitle());
  
  // Initialise histogram classes
  histos->SetReservedWords("Track;Pair");

  // Event class
  histos->AddClass("Event");

  // Track classes
  for(Int_t i = 0; i < 2; ++i){
    histos->AddClass(Form("Track_%s", AliDielectron::TrackClassName(i)));
  }
  // Pair classes
  for(Int_t i = 0; i < 3; ++i){
    histos->AddClass(Form("Pair_%s", AliDielectron::PairClassName(i)));
    // Legs of final Pairs. Both charges together. No duplicate entries.
  }
  // Mixed event pairs if mixing handler present
  if(die->GetMixingHandler()){
    histos->AddClass(Form("Pair_%s", AliDielectron::PairClassName(3)));
    histos->AddClass(Form("Pair_%s", AliDielectron::PairClassName(4)));
    histos->AddClass(Form("Pair_%s", AliDielectron::PairClassName(6)));
    histos->AddClass(Form("Pair_%s", AliDielectron::PairClassName(7)));
  }

  TH1::AddDirectory(kFALSE);

  //---------------------------------------------------------------------------
  // Add histograms to event class
  histos->UserHistogram("Event","nEvents","",1,0.,1.,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","nESDTracks","",800,0,800,AliDielectronVarManager::kNTrk);
  histos->UserHistogram("Event","zVertexPrimary","",122,-11,11,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","NVtxContrib","Number of Vertex Contributor;N of Vertex Contributors;N of events",
                        200,-0.5,199.5,AliDielectronVarManager::kNVtxContrib);
  histos->UserHistogram("Event","CentralityV0A","Centrality V0;V0A percentile",   102, -1, 101, AliDielectronVarManager::kCentralityV0A);
  histos->UserHistogram("Event","CentralityV0C","Centrality V0;V0C percentile",   102, -1, 101, AliDielectronVarManager::kCentralityV0C);

  //---------------------------------------------------------------------------
  // Single track histograms
  // Kinematic variables
  histos->UserHistogram("Track","Pt",";Pt (GeV);#tracks",100,0,5.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Px",";Px (GeV);#tracks",100,0,5.,AliDielectronVarManager::kPx);
  histos->UserHistogram("Track","Py",";Py (GeV);#tracks",100,0,5.,AliDielectronVarManager::kPy);
  histos->UserHistogram("Track","Pz",";Pz (GeV);#tracks",100,0,5.,AliDielectronVarManager::kPz);
  histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)",100,0,10,AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track","Eta","",200,-2,2,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
  // DCA
  histos->UserHistogram("Track","dXY","",400,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ" ,"",600,-4.,4.,AliDielectronVarManager::kImpactParZ);
  // ITS
  histos->UserHistogram("Track","ITSnCls",";ITS number clusters;#tracks",6,-0.5,6.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","ITSchi2",";ITS chi2/Cl;#tracks",110,0.,11.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","nITSshared","#shared ITS clusters", 7, 0, 7, AliDielectronVarManager::kNclsSITS);
  histos->UserHistogram("Track","fracITSshared","frac. shared ITS clusters", 120, 0,  1.2, AliDielectronVarManager::kNclsSITS);
  // TPC
  histos->UserHistogram("Track","TPCnCls",";TPC number clusters;#tracks",170,-0.5,169.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","TPCchi2",";TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","nTPCshared","#shared TPC clusters", 170, 0, 170, AliDielectronVarManager::kNclsSTPC);
  histos->UserHistogram("Track","NclsSFracTPC",";TPC fraction of shared clusters;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCnCrossed",";TPC findable clusters;#tracks",170,-0.5,169.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable",";TPC crossed rows over findable clusters;#tracks",240,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  // PID plots (1D)
  histos->UserHistogram("Track","nSigITSeRaw",";raw n#sigma^{ITS}_{e};#tracks", GetVector(kSigmaEle),   AliDielectronVarManager::kITSnSigmaEleRaw);
  histos->UserHistogram("Track","nSigITSe",";n#sigma^{ITS}_{e};#tracks",        GetVector(kSigmaEle),   AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","nSigTPCeRaw",";raw n#sigma^{TPC}_{e};#tracks", GetVector(kSigmaEle),   AliDielectronVarManager::kTPCnSigmaEleRaw);
  histos->UserHistogram("Track","nSigTPCe",";n#sigma^{TPC}_{e};#tracks",        GetVector(kSigmaEle),   AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","nSigTOFeRaw",";raw n#sigma^{TOF}_{e};#tracks", GetVector(kSigmaEle),   AliDielectronVarManager::kTOFnSigmaEleRaw);
  histos->UserHistogram("Track","nSigTOFe",";n#sigma^{TOF}_{e};#tracks",        GetVector(kSigmaEle),   AliDielectronVarManager::kTOFnSigmaEle);
  histos->UserHistogram("Track","nSigTPCpi",";n#sigma^{TPC}_{#pi};#tracks",     GetVector(kSigmaOther), AliDielectronVarManager::kTPCnSigmaPio);

  //---------------------------------------------------------------------------
  // Pair histograms
  histos->UserHistogram("Pair","InvMass_PairPt_Rapdity",";Inv. Mass (GeV);Pair Pt (GeV);Y_{ee}",
                        GetVector(kMee), GetVector(kPtee), BinsToVector(200, -2, 2),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kY);
  histos->UserHistogram("Pair", "InvMass_PairPt_PhiV", ";Inv. Mass (GeV);Pair Pt (GeV);PhiV",
                        GetVector(kMee), GetVector(kPtee), GetVector(kPhiV),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
  histos->UserHistogram("Pair", "InvMass_PairPt_Centrality", ";Inv. Mass (GeV);Pair Pt (GeV);Centrality",
                        GetVector(kMee), GetVector(kPtee), GetVector(kCent),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kCentralityV0A);

  die->SetHistogramManager(histos);

}

const AliDielectronEventCuts* GetEventCuts(){

  AliDielectronEventCuts* eventCuts= new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1); 
  
  return eventCuts;
}

TVectorD* GetVector(Int_t var){

  switch(var){

    case kPhiV:   return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng:  return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng2: return AliDielectronHelper::MakeLinBinning( 50, 0., TMath::Pi()/2.);
    case kPhi2D:  return AliDielectronHelper::MakeLinBinning(100, 0, 2*TMath::Pi());
    case kEta2D:  return AliDielectronHelper::MakeLinBinning(100,-1,1);
    case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);

    case kSigmaEle:
      return AliDielectronHelper::MakeLinBinning( 50, -5., 5.);
    case kSigmaOther:
      return AliDielectronHelper::MakeLinBinning( 100,-10.,20.);
    case kTPCdEdx:
      return AliDielectronHelper::MakeLinBinning( 50, 50.,100.);

    case kMee:
      return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.47, 0.62, 0.70,"
                                                         "0.77, 0.80, 0.90, 0.95, 0.99, 1.02, 1.03, 1.10, 1.40, 1.70,"
                                                         "2.00, 2.30, 2.60, 2.80, 2.90, 3.00, 3.04, 3.08, 3.10, 3.12,"
                                                         "3.20, 3.50, 5.00");
    case kMee500:
      return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,"
                                                       "0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 0.50");
    case kPtee:
      return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,"
                                                       "0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,"
                                                       "1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90,"
                                                       "2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90,"
                                                       "3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90,"
                                                       "4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 5.00, 5.50, 6.00, 6.50,"
                                                       "7.00, 8.00, 10.0");

    case kP2D:
      return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,"
                                                       "0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,"
                                                       "1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90,"
                                                       "2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90,"
                                                       "3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90,"
                                                       "4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 5.00, 5.50, 6.00, 6.50,"
                                                       "7.00, 8.00, 10.0");

    case kCent:
      return AliDielectronHelper::MakeArbitraryBinning("0, 0.5, 5.0, 10, 20, 40, 60, 80, 100");

    default: std::cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << std::endl;
      break;
  }
  return 0x0;
}

TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max){

  return AliDielectronHelper::MakeLinBinning(nbins,min,max);

}
