// Simple config file not requiring cutLibrary.
// Used for testing. Would not recommend using as base for a new config.

void InitHistograms(AliDielectron* die, Int_t cutDefinition, Bool_t doPairing);
void SetupCuts(AliDielectron* die, Int_t cutDefinition);
AliESDtrackCuts* SetupESDtrackCuts(Int_t cutDefinition);
AliDielectronPID* SetPIDcuts(Int_t cutDefinition);
const AliDielectronEventCuts* GetEventCuts();
void SetSignalsMC(AliDielectron* die);
TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max);
TVectorD* GetVector(Int_t var);
enum {kMee = 0, kMee500, kPtee, kP2D, kPhiV, kOpAng, kEta2D, kEta3D, kSigmaEle, kSigmaOther, kCent, kDCA};


AliDielectron* Config_acapon(TString name,
                             Bool_t hasMC,
                             Bool_t doPairing,
                             Bool_t doMixing)
{
  // Setup the instance of AliDielectron
  AliDielectron* die = new AliDielectron(Form("%s", name.Data()), Form("Track cuts: %s", name.Data()));

  if(hasMC){
    die->SetHasMC(hasMC);
  }
  if(!doPairing){
    die->SetNoPairing();
  }

  Int_t cutDefinition = -99;
  // Standard set of cuts and output histograms
  if(name == "kTheoPID"){
    cutDefinition = 0; // Cuts for ESDs
  }
  // Compare LF and HF DCA distributions within a MC production
  // No track cuts applied. Very basic, so far.
  else if(name = "kDCAdists"){
    cutDefinition = 1; // Basic cuts for AODs
    SetSignalsMC(die);
  
  }
  else{
    std::cout << "Invalid cut set specified!!!!" << std::endl;
    return 0x0;
  }

  AliDielectronMixingHandler* mix = new AliDielectronMixingHandler;
  mix->SetMixType(AliDielectronMixingHandler::kAll);
  mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
  mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
  mix->SetDepth(30);
  if(doMixing){
    die->SetMixingHandler(mix);
  }

	// Set track cuts
  SetupCuts(die,cutDefinition);

  // Histogram setup
  InitHistograms(die,cutDefinition, doPairing);

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

  if(cutDefinition == 0){
    allCuts->AddCut(convRejCut);
    allCuts->AddCut(pairMassCut);

    die->GetPairFilter().AddCuts(allCuts);

    die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));
    die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));
  }else if(cutDefinition == 1){
    die->GetTrackFilter().AddCuts(SetupAODtrackCuts(cutDefinition));
  }

}

//______________________________________________________________________________________
//----------------------------------- PID ----------------------------------------------
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
//----------------------------------- Track Cuts ---------------------------------------
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

AliAnalysisCuts* SetupAODtrackCuts(Int_t cutDefinition){

  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(1<<4);
  trackCutsDiel->SetRequireITSRefit(kTRUE);
  trackCutsDiel->SetRequireTPCRefit(kTRUE);

  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -0.8, 0.8);
  trackCutsAOD->AddCut(AliDielectronVarManager::kPt, 0.2, 20.);

  AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("Trackcuts","Trackcuts",AliDielectronCutGroup::kCompAND);
  trackCuts->AddCut(trackCutsAOD);
  trackCuts->AddCut(trackCutsDiel);

  trackCuts->Print();

  return trackCuts;
}

//______________________________________________________________________________________
//------------------------------- Histogram definition ---------------------------------
void InitHistograms(AliDielectron* die, Int_t cutDefinition, Bool_t doPairing)
{

  // Setup histogram classes
  AliDielectronHistos* histos = new AliDielectronHistos(die->GetName(), die->GetTitle());

  // Initialise histogram classes
  histos->SetReservedWords("Track;Pair");

  // Event class
  histos->AddClass("Event");

  // Track classes
  for(Int_t i = 0; i < 2; ++i){
    histos->AddClass(Form("Track_%s", AliDielectron::TrackClassName(i)));
  }
  // Pair classes
  if(doPairing){
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
  }

  // Add MC signal histograms
  if(die->GetMCSignals()){
    for(Int_t i = 0; i < die->GetMCSignals()->GetEntriesFast(); ++i){
      histos->AddClass(Form("Track_%s_%s", AliDielectron::PairClassName(1), die->GetMCSignals()->At(i)->GetName()));
      if(doPairing){
        histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
      }
    }
  }

  TH1::AddDirectory(kFALSE);

  //---------------------------------------------------------------------------
  // Add histograms to event class
  histos->UserHistogram("Event", "nEvents",        "", 1,   0.,   1.,    AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event", "zVertexPrimary", "", 122, -11,  11,    AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event", "NVtxContrib",    "", 200, -0.5, 199.5, AliDielectronVarManager::kNVtxContrib);
  histos->UserHistogram("Event", "nESDTracks","",                    3000, 0, 300, AliDielectronVarManager::kNTrk);
  histos->UserHistogram("Event", "NumTrack(after cuts)","",          3000, 0, 300, AliDielectronVarManager::kTracks); //run2only
  histos->UserHistogram("Event", "Num. acc. tracks","",              3000, 0, 300, AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event", "Num. ac.. Tracklets (eta<0.9)","", 3000, 0, 300, AliDielectronVarManager::kNaccTrcklts09);
  histos->UserHistogram("Event", "Num. ac.. Tracklets (eta<1)","",   3000, 0, 300, AliDielectronVarManager::kNaccTrcklts10);
  histos->UserHistogram("Event", "Num. acc. Tracklets (eta<1.6)","", 3000, 0, 300, AliDielectronVarManager::kNaccTrcklts);
  histos->UserHistogram("Event", "RefMultTPConly","",                3000, 0, 300, AliDielectronVarManager::kRefMultTPConly); //run2only
  histos->UserHistogram("Event", "nTrckltsSPD05","",                 3000, 0, 300, AliDielectronVarManager::kNaccTrckltsEsd05);
  histos->UserHistogram("Event", "nTrckltsSPD10","",                 3001, 0, 300, AliDielectronVarManager::kNaccTrckltsEsd10);
  histos->UserHistogram("Event", "nTrItsPureESD05","",               3000, 0, 300, AliDielectronVarManager::kNaccItsPureEsd05);
  histos->UserHistogram("Event", "nTrItsPureESD10","",               3000, 0, 300, AliDielectronVarManager::kNaccItsPureEsd10);
  histos->UserHistogram("Event", "nTrITSTPC05","",                   3000, 0, 300, AliDielectronVarManager::kNaccItsTpcEsd05Corr);
  histos->UserHistogram("Event", "nTrITSTPC10","",                   3000, 0, 300, AliDielectronVarManager::kNaccItsTpcEsd10);

  //---------------------------------------------------------------------------
  // Single track histograms
  // Kinematic variables
  histos->UserHistogram("Track", "Pt",    "", 100, 0,  5.,             AliDielectronVarManager::kPt);
  histos->UserHistogram("Track", "Px",    "", 100, 0,  5.,             AliDielectronVarManager::kPx);
  histos->UserHistogram("Track", "Py",    "", 100, 0,  5.,             AliDielectronVarManager::kPy);
  histos->UserHistogram("Track", "Pz",    "", 100, 0,  5.,             AliDielectronVarManager::kPz);
  histos->UserHistogram("Track", "P_PIn", "", 100, 0,  10,             AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track", "Eta",   "", 200, -2, 2,              AliDielectronVarManager::kEta);
  histos->UserHistogram("Track", "Phi",   "", 120, 0., TMath::TwoPi(), AliDielectronVarManager::kPhi);
  // DCA
  histos->UserHistogram("Track", "dXY", "",  400, -2., 2., AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track", "dZ",   "", 600, -4., 4., AliDielectronVarManager::kImpactParZ);
  // ITS
  histos->UserHistogram("Track", "ITSnCls",       "", 6,   -0.5, 6.5, AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track", "ITSchi2",       "", 110, 0.,   11., AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track", "nITSshared",    "", 7,   0,    7,   AliDielectronVarManager::kNclsSITS);
  histos->UserHistogram("Track", "fracITSshared", "", 120, 0,    1.2, AliDielectronVarManager::kNclsSITS);
  // TPC
  histos->UserHistogram("Track", "TPCnCls",                    "", 170, -0.5, 169.5, AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track", "TPCchi2",                    "", 100, 0.,   10.,   AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track", "nTPCshared",                 "", 170, 0,    170,   AliDielectronVarManager::kNclsSTPC);
  histos->UserHistogram("Track", "NclsSFracTPC",               "", 200, 0,    10.,   AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track", "TPCnCrossed",                "", 170, -0.5, 169.5, AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track", "TPCcrossedRowsOverFindable", "", 240, 0.,   1.2,   AliDielectronVarManager::kNFclsTPCfCross);
  // PID plots (1D)
  histos->UserHistogram("Track", "nSigITSeRaw", "", GetVector(kSigmaEle),   AliDielectronVarManager::kITSnSigmaEleRaw);
  histos->UserHistogram("Track", "nSigITSe",    "", GetVector(kSigmaEle),   AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track", "nSigTPCeRaw", "", GetVector(kSigmaEle),   AliDielectronVarManager::kTPCnSigmaEleRaw);
  histos->UserHistogram("Track", "nSigTPCe",    "", GetVector(kSigmaEle),   AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track", "nSigTOFeRaw", "", GetVector(kSigmaEle),   AliDielectronVarManager::kTOFnSigmaEleRaw);
  histos->UserHistogram("Track", "nSigTOFe",    "", GetVector(kSigmaEle),   AliDielectronVarManager::kTOFnSigmaEle);
  histos->UserHistogram("Track", "nSigTPCpi",   "", GetVector(kSigmaOther), AliDielectronVarManager::kTPCnSigmaPio);

  histos->UserHistogram("Track", "pdgCode",       "", 4000, -2000, 2000, AliDielectronVarManager::kPdgCode);
  histos->UserHistogram("Track", "pdgCodeMother", "", 4000, -2000, 2000, AliDielectronVarManager::kPdgCodeMother);


  //---------------------------------------------------------------------------
  // Pair histograms
  if(doPairing){
    histos->UserHistogram("Pair", "InvMass",      "", 500, 0.,  5.,          AliDielectronVarManager::kM);
    histos->UserHistogram("Pair", "PairPt",       "", 160, 0.,  8.,          AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair", "Rapidity",     "", 200, -2., 2.,          AliDielectronVarManager::kY);
    histos->UserHistogram("Pair", "OpeningAngle", "", 240, 0.,  TMath::Pi(), AliDielectronVarManager::kOpeningAngle);

    histos->UserHistogram("Pair", "InvMass_PairPt_Rapdity", "", GetVector(kMee), GetVector(kPtee), BinsToVector(200, -2, 2),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kY);
    histos->UserHistogram("Pair", "InvMass_PairPt_PhiV", "", GetVector(kMee), GetVector(kPtee), GetVector(kPhiV),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair", "InvMass_PairPt_Centrality", "", GetVector(kMee), GetVector(kPtee), GetVector(kCent),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kCentralityV0A);
    histos->UserHistogram("Pair", "InvMass_pPt_PairDCA","", GetVector(kMee),GetVector(kPtee), GetVector(kDCA),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPairDCAsigXY);

    histos->UserHistogram("Pair", "PairDCAsigXY", "", GetVector(kDCA), AliDielectronVarManager::kPairDCAsigXY);
    histos->UserHistogram("Pair", "PairDCAsigZ",  "", GetVector(kDCA), AliDielectronVarManager::kPairDCAsigZ);
    histos->UserHistogram("Pair", "PairDCAabsXY", "", GetVector(kDCA), AliDielectronVarManager::kPairDCAabsXY);
    histos->UserHistogram("Pair", "PairDCAabsZ",  "", GetVector(kDCA), AliDielectronVarManager::kPairDCAabsZ);
  }

  die->SetHistogramManager(histos);

}

//______________________________________________________________________________________
//-------------------------------- Set Event Cuts --------------------------------------
const AliDielectronEventCuts* GetEventCuts(){

  AliDielectronEventCuts* eventCuts= new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD);
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1);
  eventCuts->Print();

  eventCuts->SetRequireAliEventCuts(kTRUE);

  return eventCuts;
}


//______________________________________________________________________________________
//----------------------------- Define MC Signals --------------------------------------
void SetSignalsMC(AliDielectron* die){

  // Used pdg codes (defined in AliDielectronMC::ComparePDG)
  // 401: open charm meson
  // 404: charged open charmed mesons NO s quark
  // 405: neutral open charmed mesons
  // 406: charged open charmed mesons with s quark
  // 501: open beauty mesons
  // 503: all beauty hadrons
  // 504: charged open beauty mesons NO s quark
  // 505: neutral open beauty mesons
  // 506: charged open beauty mesons with s quark
  // all D mesons

  // decay channels
  // (1) D -> e X
  // (1) B -> e X
  // (2) B -> D X -> e X Y
  // (3) B -> e D X -> ee X Y always produces ULS pair

  // Electrons from open beauty mesons and baryons
  AliDielectronSignalMC* eleFinalStateFromB = new AliDielectronSignalMC("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB->SetLegPDGs(11,-11);
  eleFinalStateFromB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB->SetMotherPDGs(502, 502);
  eleFinalStateFromB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  eleFinalStateFromB->SetCheckCorrelatedHF(kTRUE);
  die->AddSignalMC(eleFinalStateFromB);

  // Electrons from open charm mesons and baryons
  AliDielectronSignalMC* eleFinalStateFromD = new AliDielectronSignalMC("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD->SetLegPDGs(11,-11);
  eleFinalStateFromD->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD->SetMotherPDGs(402, 402);
  eleFinalStateFromD->SetCheckBothChargesMothers(kTRUE,kTRUE);
  eleFinalStateFromD->SetCheckCorrelatedHF(kTRUE);
  die->AddSignalMC(eleFinalStateFromD);

  /* // D+- meson (1)(1) */
  /* // Dielectrons originating from open charm hadrons */
  /* AliDielectronSignalMC* diEleOpenCharmCharged = new AliDielectronSignalMC("DmesonsCharged","di-electrons from open charm D+- mesons no B grandmother"); */
  /* diEleOpenCharmCharged->SetLegPDGs(11,-11); */
  /* diEleOpenCharmCharged->SetMotherPDGs(401,401); */
  /* diEleOpenCharmCharged->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); */
  /* diEleOpenCharmCharged->SetMothersRelation(AliDielectronSignalMC::kDifferent); */
  /* diEleOpenCharmCharged->SetCheckStackForPDG(kTRUE); */
  /* diEleOpenCharmCharged->SetPDGforStack(503); */
  /* diEleOpenCharmCharged->SetCheckBothChargesLegs(kTRUE,kTRUE); */
  /* diEleOpenCharmCharged->SetCheckBothChargesMothers(kTRUE,kTRUE); */
  /* die->AddSignalMC(diEleOpenCharmCharged); */

  /* // D0 meson (1)(1) */
  /* // Dielectrons originating from open charm hadrons */
  /* AliDielectronSignalMC* diEleOpenCharmNeutral = new AliDielectronSignalMC("DmesonsNeutral","di-electrons from open charm D0 mesons no B grandmother"); */
  /* diEleOpenCharmNeutral->SetLegPDGs(11,-11); */
  /* diEleOpenCharmNeutral->SetMotherPDGs(405,405); */
  /* diEleOpenCharmNeutral->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); */
  /* diEleOpenCharmNeutral->SetMothersRelation(AliDielectronSignalMC::kDifferent); */
  /* diEleOpenCharmNeutral->SetCheckStackForPDG(kTRUE); */
  /* diEleOpenCharmNeutral->SetPDGforStack(503); */
  /* diEleOpenCharmNeutral->SetCheckBothChargesLegs(kTRUE,kTRUE); */
  /* diEleOpenCharmNeutral->SetCheckBothChargesMothers(kTRUE,kTRUE); */
  /* die->AddSignalMC(diEleOpenCharmNeutral); */

  /* //B meson (3) */
  /* // Dielectrons originating from open charm hadrons */
  /* AliDielectronSignalMC* diEleOneOpenB = new AliDielectronSignalMC("B2ee","di-electrons from one B meson"); */
  /* diEleOneOpenB->SetLegPDGs(11,-11); */
  /* diEleOneOpenB->SetMotherPDGs(401,501); */
  /* diEleOneOpenB->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); */
  /* diEleOneOpenB->SetGrandMotherPDGs(501,0); */
  /* diEleOneOpenB->SetCheckMotherGrandmotherRelation(kTRUE,kTRUE); */
  /* diEleOneOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE); */
  /* diEleOneOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE); */
  /* diEleOneOpenB->SetCheckBothChargesGrandMothers(kTRUE,kTRUE); */
  /* die->AddSignalMC(diEleOneOpenB); */

  /* // B meson (1)(1) */
  /* // Dielectrons originating from open charm hadrons */
  /* AliDielectronSignalMC* diEleOpenB = new AliDielectronSignalMC("BMesons","di-electrons from B mesons"); */
  /* diEleOpenB->SetLegPDGs(11,-11); */
  /* diEleOpenB->SetMotherPDGs(501,501); */
  /* diEleOpenB->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); */
  /* diEleOpenB->SetMothersRelation(AliDielectronSignalMC::kDifferent); */
  /* diEleOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE); */
  /* diEleOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE); */
  /* die->AddSignalMC(diEleOpenB); */

  /* // B meson (2)(2) */
  /* // Dielectrons originating from open charm hadrons */
  /* AliDielectronSignalMC* diEleOpenBtoD = new AliDielectronSignalMC("B2D2ee","di-electrons from B->D-> e"); */
  /* diEleOpenBtoD->SetLegPDGs(11,-11); */
  /* diEleOpenBtoD->SetMotherPDGs(401,401); */
  /* diEleOpenBtoD->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); */
  /* diEleOpenBtoD->SetGrandMotherPDGs(501,501); */
  /* diEleOpenBtoD->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent); */
  /* diEleOpenBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE); */
  /* diEleOpenBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE); */
  /* diEleOpenBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE); */
  /* die->AddSignalMC(diEleOpenBtoD); */

  /* // B meson (1)(2) */
  /* // Dielectrons originating from open charm hadrons */
  /* AliDielectronSignalMC* diEleOpenBandBtoD = new AliDielectronSignalMC("B2eAndB2D2e","di-electrons from B->e and B->D->e"); */
  /* diEleOpenBandBtoD->SetLegPDGs(11,11); */
  /* diEleOpenBandBtoD->SetMotherPDGs(401,501); */
  /* diEleOpenBandBtoD->SetGrandMotherPDGs(501,0); */
  /* diEleOpenBandBtoD->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); */
  /* diEleOpenBandBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE); */
  /* diEleOpenBandBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE); */
  /* diEleOpenBandBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE); */
  /* diEleOpenBandBtoD->SetCheckMotherGrandmotherRelation(kTRUE,kFALSE); */
  /* die->AddSignalMC(diEleOpenBandBtoD); */

  // Dielectrons originating from dalitz decays
  AliDielectronSignalMC* PiDalitz = new AliDielectronSignalMC("Pi0","di-electrons from Pi0 dalitz");
  PiDalitz->SetLegPDGs(11,-11);
  PiDalitz->SetMotherPDGs(111,111);
  PiDalitz->SetMothersRelation(AliDielectronSignalMC::kSame);
  PiDalitz->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PiDalitz->SetCheckBothChargesLegs(kTRUE,kTRUE);
  PiDalitz->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(PiDalitz);

  AliDielectronSignalMC* PiDalitzNoFeedDown = new AliDielectronSignalMC("Pi0NoFeedDown","di-electrons from Pi0 dalitz no feeddown from Ks");  // dielectrons originating from dalitz decays
  PiDalitzNoFeedDown->SetLegPDGs(11,-11);
  PiDalitzNoFeedDown->SetMotherPDGs(111,111);
  PiDalitzNoFeedDown->SetMothersRelation(AliDielectronSignalMC::kSame);
  PiDalitzNoFeedDown->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PiDalitzNoFeedDown->SetCheckStackForPDG(kTRUE);
  PiDalitzNoFeedDown->SetPDGforStack(310);
  PiDalitzNoFeedDown->SetCheckBothChargesLegs(kTRUE,kTRUE);
  PiDalitzNoFeedDown->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(PiDalitzNoFeedDown);
}

TVectorD* GetVector(Int_t var){

  switch(var){

    case kPhiV:   return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng:  return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kEta2D:  return AliDielectronHelper::MakeLinBinning(100,-1,1);
    case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);

    case kSigmaEle:
      return AliDielectronHelper::MakeLinBinning( 50, -5., 5.);
    case kSigmaOther:
      return AliDielectronHelper::MakeLinBinning( 100,-10.,20.);

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
    case kDCA:
      return AliDielectronHelper::MakeLinBinning(50, 0., 20.);

    default: std::cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << std::endl;
      break;
  }
  return 0x0;
}

TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max){

  return AliDielectronHelper::MakeLinBinning(nbins,min,max);

}
