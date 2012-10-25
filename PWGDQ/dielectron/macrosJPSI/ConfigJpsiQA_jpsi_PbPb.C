void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);
void InitHF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

void AddMCSignals(AliDielectron *die);
void SetEtaCorrection();
TVectorD *GetRunNumbers();

TString names=("QA");
enum { kQA=0 };

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t  isESD = kTRUE;
Bool_t  hasMC = kFALSE;
TString list  = gSystem->Getenv("LIST");


AliDielectron* ConfigJpsiQA_jpsi_PbPb(Int_t cutDefinition, TString prod="")
{
  //
  // Setup the instance of AliDielectron
  //

  // find mc or not?
  if( list.IsNull()) list=prod;
  if( list.Contains("LHC10h")   || list.Contains("LHC11h")   ) hasMC=kFALSE;
  if( list.Contains("LHC11a10") || list.Contains("LHC12a17") ) hasMC=kTRUE;

  //ESD handler?
  isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());

  // switch off some configurations
  if(hasMC) { // MONTE CARLO
    switch(cutDefinition) {
      //    case kQA:        return 0x0;
    }
  } else { // COLLISION DATA
    switch(cutDefinition) {
      //    case kQA:  return 0x0;
    }
  }

  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("Track cuts: %s",name.Data()));
  die->SetHasMC(hasMC);

  printf(" Adding %s %s config %s for %s \n",(isESD?"ESD":"AOD"),(hasMC?"MC":""),name.Data(),list.Data());

  // Monte Carlo Signals and TRD efficiency tables
  if(hasMC) {
    AddMCSignals(die);
    printf(" Add %d MC signals \n",die->GetMCSignals()->GetEntriesFast());

    // trd tables
    /* if (list.Contains("LHC11a") ) {
       TString pidTab="$TRAIN_ROOT/util/dielectron/dielectron/TRDpidEff_eleProb07_TRDntr4_6.root";
       TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
       if (trainRoot.IsNull()) pidTab="$ALICE_ROOT/PWGDQ/dielectron/files/TRDpidEff_eleProb07_TRDntr4_6.root";

       if (gSystem->AccessPathName(gSystem->ExpandPathName(pidTab.Data())))
       Error("ConfigPbPb","PID table not found: %s",pidTab.Data());
       else
       die->SetTRDcorrectionFilename(pidTab.Data());
       }*/
  }

  // histogram setup
  InitHistograms(die,cutDefinition);
  printf(" Add %d classes to the manager \n",die->GetHistogramList()->GetEntries());

  // cut setup
  SetupTrackCuts(die,cutDefinition);
  //  SetupPairCuts(die,cutDefinition);


  // CF container setup, switched off
  //InitCF(die,cutDefinition);
  if(die->GetCFManagerPair())
    printf(" Add %d pair, %d leg vars, %p steps and %p bins to the container \n",
	   die->GetCFManagerPair()->GetNvarsPair(),die->GetCFManagerPair()->GetNvarsLeg(),
	   die->GetCFManagerPair()->GetContainer(), die->GetCFManagerPair()->GetContainer() );

  // HF arrays setup
  // InitHF(die,cutDefinition);

  /*
  // tpc event plane
  if(!hasMC) {
    // TPC event plane configurations
    Double_t gGap;
    switch(cutDefinition) {
    default: gGap=0.0;
    }

    AliDielectronVarCuts *poi = new AliDielectronVarCuts("PoI","PoI");
    poi->AddCut(AliDielectronVarManager::kM,2.92,3.20);     // particles of interest, jpsi mass window
    die->GetEventPlanePOIPreFilter().AddCuts(poi);

    AliDielectronVarCuts *etaGap = new AliDielectronVarCuts(AliDielectronVarManager::GetValueName(AliDielectronVarManager::kEta),"etaGap");
    etaGap->AddCut(AliDielectronVarManager::kEta,-1*gGap,gGap,kTRUE);
    die->GetEventPlanePreFilter().AddCuts(etaGap);
    // die->SetLikeSignSubEvents();

    die->SetPreFilterEventPlane();
  }
  */
  // prefilter settings
  die->SetNoPairing();
  //die->SetPreFilterUnlikeOnly();
  //die->SetPreFilterAllSigns();

  // setup eta correction
  //  if(isESD && list.Contains("LHC10h")) SetEtaCorrection();

  // VZERO calibration
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (!trainRoot.IsNull()) {
    die->SetVZEROCalibrationFilename("$TRAIN_ROOT/util/dielectron/dielectron/VzeroCalibrationLHC10h.root");
    die->SetVZERORecenteringFilename("$TRAIN_ROOT/util/dielectron/dielectron/VzeroRecenteringLHC10h.root");
  }

  return die;
}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  //  cuts->AddDielectron(die);
  die->GetTrackFilter().AddCuts(cuts);

  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  //  trkFilter->AddCutQA(die->GetHistoManager());
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  //  trkFilter->SetMinNCrossedRowsOverFindable(0.6);
  if(!isESD) cuts->AddCut(trkFilter);

  //Pt cut, should make execution a bit faster
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("PtCut","PtCut");
  //  pt->AddCutQA(die->GetHistoManager());
  pt->AddCut(AliDielectronVarManager::kPt,0.8,1e30);    //1.1
  cuts->AddCut(pt);
  pt->Print();

  // track cuts ESD and AOD
  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");
  //  varCuts->AddCutQA(die->GetHistoManager());
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  //varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     50.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.0);
  //varCuts->AddCut(AliDielectronVarManager::kTOFbeta,      0.2,   0.9, kTRUE);
  varCuts->Print();
  cuts->AddCut(varCuts);


  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  //  trkCuts->AddCutQA(die->GetHistoManager());
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 15); // ITS-4 = 1+2+4+8
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);
  cuts->AddCut(trkCuts);


  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronPID *pid = new AliDielectronPID("PID","PID");
  //  pid->AddCutQA(die->GetHistoManager());
  ////////////////////////////////// DATA
  if(!hasMC) {
    // TPC
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,   -4.,4.0);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -100.,4.0,0.,0.,kTRUE);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -100.,3.5,0.,0.,kTRUE);

    // TOF
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-5,5.,0.,0.,kFALSE,AliDielectronPID::kIfAvailable);
    //    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,AliDielectronPID::kIfAvailable);

    // TRD 1- or 2-dimensonal
    /*    pid->AddCut(AliDielectronPID::kTRDeleEff,AliPID::kElectron,.8,1.,3.5.,6.,kFALSE,
	  AliDielectronPID::kIfAvailable,AliDielectronVarManager::kTRDpidQuality);
	  pid->AddCut(AliDielectronPID::kTRDeleEff2D,AliPID::kElectron,.8,1.,3.5.,6.,kFALSE,
	  AliDielectronPID::kIfAvailable,AliDielectronVarManager::kTRDpidQuality);
    */
  }

  ////////////////////////////////// MC
  if(hasMC) {

    // electron
    Double_t nSigmaPi = 3.5; Double_t nSigmaP = 3.5;
    Double_t resolution=0.0549;
    Double_t BBpro[5] = {0};
    Double_t BBpio[5] = {0};

    for(Int_t icent=0; icent<8; icent++) {

      switch (icent) {
      case 0:  // 0-10%
	BBpro[0] = 0.031555;  BBpro[1] = 26.0595; BBpro[2] = 3.02422e-11;  BBpro[3] = 2.05594; BBpro[4] = 5.99848;
	BBpio[0] = 0.0252122; BBpio[1] = 38.8991; BBpio[2] = 4.0901e-11;   BBpio[3] = 5.27988; BBpio[4] = 4.3108;
	break;
      case 1:  // 10-20%
	BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
	BBpio[0] = 0.0252127; BBpio[1] = 33.8617; BBpio[2] = 3.56866e-11;  BBpio[3] = 5.24831; BBpio[4] = 4.31093;
	break;
      case 2:  // 20-30%
	BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
	BBpio[0] = 0.0263205; BBpio[1] = 37.9307; BBpio[2] = 4.29724e-11;  BBpio[3] = 5.74458; BBpio[4] = 4.32459;
	break;
      case 3:  // 30-40%
	BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
	BBpio[0] = 0.026294;  BBpio[1] = 39.0346; BBpio[2] = 4.12261e-11;  BBpio[3] = 5.28808; BBpio[4] = 4.31301;
	break;
      case 4:  // 40-50%
	BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
	BBpio[0] = 0.0263134; BBpio[1] = 38.2084; BBpio[2] = 3.75159e-11;  BBpio[3] = 5.78125; BBpio[4] = 4.31363;
	break;
      case 5:  // 50-60%
	BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
	BBpio[0] = 0.0263134; BBpio[1] = 38.2084; BBpio[2] = 3.75159e-11;  BBpio[3] = 5.78125; BBpio[4] = 4.31363;
	break;
      case 6:  // 60-70%
	BBpro[0] = 0.031555;  BBpro[1] = 26.0595; BBpro[2] = 3.02422e-11;  BBpro[3] = 2.05594; BBpro[4] = 5.99848;
	BBpio[0] = 0.026302;  BBpio[1] = 38.6888; BBpio[2] = 3.56792e-11;  BBpio[3] = 5.2465;  BBpio[4] = 4.31094;
	break;
      case 7:  // 70-80%
	BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
	BBpio[0] = 0.0263134; BBpio[1] = 38.2084; BBpio[2] = 3.75159e-11;  BBpio[3] = 5.78125; BBpio[4] = 4.31363;
	break;
      case 8:  // 80-90%
	BBpro[0] = 0.0313438; BBpro[1] = 25.8666; BBpro[2] = 4.5457e-11;   BBpro[3] = 2.07912; BBpro[4] = 5.99986;
	BBpio[0] = 0.0252127; BBpio[1] = 33.8617; BBpio[2] = 3.56866e-11;  BBpio[3] = 5.24831; BBpio[4] = 4.31093;
	break;
      case 9:  // 90-100%
	BBpro[0] = 0.0319126; BBpro[1] = 36.8784; BBpro[2] = 3.4274e-11;   BBpro[3] = 3.2431;  BBpro[4] = 5.93388;
	BBpio[0] = 0.027079;  BBpio[1] = 67.5936; BBpio[2] = 9.72548e-11;  BBpio[3] = 9.61382; BBpio[4] = 5.99372;
	break;
      }


      TF1 *ffPro=new TF1(Form("fBethe%d_c%d",AliPID::kProton,icent), Form("(%f*%f+(AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])-AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])))/%f", nSigmaP,resolution, AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kElectron), resolution), 0.05,200.);

      TF1 *ffPio=new TF1(Form("fBethe%d_c%d",AliPID::kPion,icent), Form("(%f*%f+(AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])-AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])))/%f", nSigmaPi,resolution, AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kElectron), resolution), 0.05,200.);

      //LHC11a10b
      if (list.Contains("LHC11a10b")) {
        printf(" LHC11a10b parameters\n");
        ffPro->SetParameters(BBpro[0],BBpro[1],BBpro[2],BBpro[3],BBpro[4]);
        ffPio->SetParameters(BBpio[0],BBpio[1],BBpio[2],BBpio[3],BBpio[4]);

        // proton cut
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,ffPro,10,((double)icent)*10.,((double)icent+1)*10,
                    kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kCentrality);
        // pion cut
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,ffPio,10,((double)icent)*10.,((double)icent+1)*10,
                    kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kCentrality);
      }
    }

    // shifts for the nSigma electrons
    TGraph* nSigmaCorrection = new TGraph();
    // LHC11a10b
    if (list.Contains("LHC11a10b")) {
      nSigmaCorrection->SetPoint(0, 137161., -0.50-(0.28));
      nSigmaCorrection->SetPoint(1, 139510., -0.50-(0.28));
      pid->SetCorrGraph(nSigmaCorrection);
    }

  } //hasMC

  cuts->AddCut(pid);
  /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ PID CUTS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/


  // exclude conversion electrons selected by the tender
  AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  cuts->AddCut(noconv);

}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //

  // conversion rejection
  Double_t gCut;
  switch(cutDefinition) {
  case kQA:    gCut=0.0;  break;
  default: gCut=0.05;             // default
  }

  AliDielectronVarCuts *gammaCuts = new AliDielectronVarCuts("GammaCuts","GammaCuts");
  //  gammaCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.1,  kTRUE);
  //  gammaCuts->AddCut(AliDielectronVarManager::kLegDist,      0.0,   0.25, kTRUE);
  //  gammaCuts->AddCut(AliDielectronVarManager::kR,            3.0,   90.0, kTRUE);
  //  gammaCuts->AddCut(AliDielectronVarManager::kPsiPair,      0.0,   0.05, kTRUE);
  //  gammaCuts->AddCut(AliDielectronVarManager::kChi2NDF,      0.0,   10.0, kTRUE);
  gammaCuts->AddCut(AliDielectronVarManager::kM,            0.0,   gCut);
  gammaCuts->Print();
  die->GetPairPreFilter().AddCuts(gammaCuts);


  // rapidity selection
  //  AliDielectronVarCuts *rapCut=new AliDielectronVarCuts("|Y|<.9","|Y|<.9");
  // rapCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  // die->GetPairFilter().AddCuts(rapCut);

}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //

  // booleans for histo selection
  Bool_t bHistEvtQA=kFALSE, bHistTrackQA=kFALSE, bHistPairQA = kFALSE;
  switch (cutDefinition) {
  case kQA:       bHistEvtQA=kTRUE; bHistTrackQA=kTRUE; break;
  }


  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());


  ////// EVENT HISTOS /////
  if(bHistEvtQA) {

    //add histograms to event class
    histos->AddClass("Event");

    histos->UserHistogram("Event","RunNumber_Centrality","",
			  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(80,0.,80.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kCentrality);
    histos->UserProfile("Event","RunNumber-Centrality","", AliDielectronVarManager::kCentrality,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","RunNumber-CentralitySPD","", AliDielectronVarManager::kCentralitySPD,
    			GetRunNumbers(), AliDielectronVarManager::kRunNumber);

    histos->UserHistogram("Event","RunNumber_VtxZ","",
			  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(150,-15.,15.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kZvPrim);
    histos->UserProfile("Event","RunNumber-VtxZ","", AliDielectronVarManager::kZvPrim,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","RunNumber-sigVtxZ","", AliDielectronVarManager::kZvPrim,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber, "s;-10;10");
    histos->UserProfile("Event","Centrality-VtxZ","", AliDielectronVarManager::kZvPrim,
			AliDielectronHelper::MakeLinBinning(80,0.,80.), AliDielectronVarManager::kCentrality);

    histos->UserHistogram("Event","RunNumber_v0ArpH2","",
			  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(100,-2.,+2.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kv0ArpH2);
    histos->UserHistogram("Event","RunNumber_v0CrpH2","",
			  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(100,-2.,+2.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kv0CrpH2);
    //     histos->UserHistogram("Event","RunNumber_TPCrpH2",";run;#Psi_{2}^{TPC} (rad.)",
    // 			  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(100,-2.,+2.),
    // 			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kTPCrpH2);

  }

  ////// EVENT HISTOS /////
  if(bHistTrackQA) {

    //add histograms to track class
    histos->SetReservedWords("Track");
    for (Int_t i=0; i<2; ++i) histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));

    histos->UserHistogram("Track","Pt", "", 400, 0.,20., AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","Eta","", 200,-1.,+1., AliDielectronVarManager::kEta);
    histos->UserHistogram("Track","Phi","", 180,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    histos->UserProfile("Track","RunNumber-Pt","", AliDielectronVarManager::kPt,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","RunNumber-Eta","", AliDielectronVarManager::kEta,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","RunNumber-Phi","", AliDielectronVarManager::kPhi,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);

    histos->UserHistogram("Track","ImpactParXY", "", 400,-1.,+1., AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","ImpactParZ",  "", 600,-3.,+3., AliDielectronVarManager::kImpactParZ);
    histos->UserProfile("Track","Eta_Phi-ImpactParXY", "", AliDielectronVarManager::kImpactParXY,
			200,-1.,+1., 180,0.,TMath::TwoPi(),
			AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    histos->UserProfile("Track","Eta_Phi-ImpactParZ", "", AliDielectronVarManager::kImpactParZ,
			200,-1.,+1., 180,0.,TMath::TwoPi(),
			AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    histos->UserProfile("Track","RunNumber-ImpactParXY","", AliDielectronVarManager::kImpactParXY,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","RunNumber-ImpactParZ","", AliDielectronVarManager::kImpactParZ,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);


    // TPC
    histos->UserHistogram("Track","NclsTPC", "", 160,0.,160., AliDielectronVarManager::kNclsTPC);
    histos->UserProfile("Track","RunNumber-NclsTPC","", AliDielectronVarManager::kNclsTPC,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","Cent-NclsTPC","", AliDielectronVarManager::kNclsTPC,
			80, 0.,80., AliDielectronVarManager::kCentrality);
    histos->UserProfile("Track","Eta_Phi-NclsTPC", "", AliDielectronVarManager::kNclsTPC,
			200,-1.,+1., 180,0.,TMath::TwoPi(),
			AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);

    histos->UserHistogram("Track","TPCclsSegments", "", 8,0.,8., AliDielectronVarManager::kTPCclsSegments);
    histos->UserProfile("Track","RunNumber-TPCclsSegments","", AliDielectronVarManager::kTPCclsSegments,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","Cent-TPCclsSegments","", AliDielectronVarManager::kTPCclsSegments,
			80, 0.,80., AliDielectronVarManager::kCentrality);
    histos->UserProfile("Track","Eta_Phi-TPCclsSegments", "", AliDielectronVarManager::kTPCclsSegments,
			200,-1.,+1., 180,0.,TMath::TwoPi(),
			AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    
    // ITS
    histos->UserHistogram("Track","NclsITS", "", 7,0.,7., AliDielectronVarManager::kNclsITS);
    histos->UserProfile("Track","RunNumber-NclsITS","", AliDielectronVarManager::kNclsITS,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","Cent-NclsITS","", AliDielectronVarManager::kNclsITS,
			80, 0.,80., AliDielectronVarManager::kCentrality);
    histos->UserProfile("Track","Eta_Phi-NclsITS", "", AliDielectronVarManager::kNclsITS,
			200,-1.,+1., 180,0.,TMath::TwoPi(),
			AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    
    // TPC PID
    histos->UserProfile("Track","RunNumber-TPCnSigmaEle","", AliDielectronVarManager::kTPCnSigmaEle,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","RunNumber_Cent_TPCnSigmaEle","",
			  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(16, 0.,80.), AliDielectronHelper::MakeLinBinning(100, -5.,+5),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","RunNumber_Eta_TPCnSigmaEle","",
			  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(20, -1.,1.), AliDielectronHelper::MakeLinBinning(100, -5.,+5),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","RunNumber_Phi_TPCnSigmaEle","",
			  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(18, 0.,TMath::TwoPi()), AliDielectronHelper::MakeLinBinning(100, -5.,+5),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kPhi, AliDielectronVarManager::kTPCnSigmaEle);

    // TOF PID
    histos->UserProfile("Track","RunNumber-TOFnSigmaEle","", AliDielectronVarManager::kTOFnSigmaEle,
			GetRunNumbers(), AliDielectronVarManager::kRunNumber, "h;-10;+10");
    histos->UserProfile("Track","RunNumber_Cent-TOFnSigmaEle","", AliDielectronVarManager::kTOFnSigmaEle,
			GetRunNumbers(), AliDielectronHelper::MakeLinBinning(16, 0.,80.),
			AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kCentrality,"h;-10;+10");
    histos->UserProfile("Track","Eta_Phi-TOFnSigmaEle","", AliDielectronVarManager::kTOFnSigmaEle,
			AliDielectronHelper::MakeLinBinning(20, -1.,+1.), AliDielectronHelper::MakeLinBinning(18, 0.,TMath::TwoPi()),
			AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi, "h;-10;+10");
    histos->UserHistogram("Track","TOFnSigmaEle", "", 100,-10.,+10., AliDielectronVarManager::kTOFnSigmaEle);
  }

/*


    histos->UserHistogram("Event","RunNumber","Events per run;run;events",
    GetRunNumbers(),
    AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Event","RunNumber_VtxZ",";run;z_{vtx} (cm)",
    GetRunNumbers(), AliDielectronHelper::MakeLinBinning(150,-15.,15.),
    AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","RunNumber_Multiplicity",";run;multiplicity V0",
    GetRunNumbers(), AliDielectronHelper::MakeLinBinning(250,0.,25000.),
    AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","RunNumber_v0ArpH2",";run;#Psi_{2}^{V0A} (rad.)",
    GetRunNumbers(), AliDielectronHelper::MakeLinBinning(100,-2.,+2.),
    AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kv0ArpH2);
    histos->UserHistogram("Event","RunNumber_v0CrpH2",";run;#Psi_{2}^{V0C} (rad.)",
    GetRunNumbers(), AliDielectronHelper::MakeLinBinning(100,-2.,+2.),
    AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kv0CrpH2);
    histos->UserHistogram("Event","RunNumber_TPCrpH2",";run;#Psi_{2}^{TPC} (rad.)",
    GetRunNumbers(), AliDielectronHelper::MakeLinBinning(100,-2.,+2.),
    AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kTPCrpH2);

    histos->UserHistogram("Event","VtxZ","Vertex Z;z_{vtx} (cm)", 300,-15.,15.,
    AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","Multiplicity","Multiplicity V0;Multiplicity V0;events", 500,0.,25000.,
    AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","Cent_Mult","Centrality vs. Multiplicity;centrality (%);Multiplicity V0",
    100,0.,100., 250,0.,25000.,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kMultV0);
    histos->UserProfile("Event","Cent_Nacc", "accepted tracks;centrality (%)",
    AliDielectronVarManager::kNacc,
    "0.,5.,10.,20.,40.,50.,60.,80.,100.",
    AliDielectronVarManager::kCentrality);
    histos->UserProfile("Event","Cent_NVtxContrib", "number of vertex contributors;centrality (%)",
    AliDielectronVarManager::kNVtxContrib,
    100,0.,100.,
    AliDielectronVarManager::kCentrality);
    } //hist: event


    ////// FLOW //////
    if(bHistFlow) {

    // RP angles versus centrality
    histos->UserHistogram("Event","Cent_TPCrpH2","TPC RP;centrality (%);#Psi^{TPC} (rad.)",
    10,0.,100.,100,-2.,2.,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCrpH2);
    histos->UserHistogram("Event","Cent_TPCsub1rpH2","TPC-1 RP;centrality (%);#Psi^{sub1} (rad.)",
    10,0.,100.,100,-2.,2.,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsub1rpH2);
    histos->UserHistogram("Event","Cent_TPCsub2rpH2","TPC-2 RP;centrality (%);#Psi^{sub2} (rad.)",
    10,0.,100.,100,-2.,2.,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsub2rpH2);

    histos->UserHistogram("Event","Cent_v0ACrpH2","VZERO-AC RP;centrality (%);#Psi_{2}^{V0AC} (rad.)",
    10,0.,100.,100,-2.0,2.0,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0ACrpH2);
    histos->UserHistogram("Event","Cent_v0ArpH2","VZERO-A RP;centrality (%);#Psi_{2}^{V0A} (rad.)",
    10,0.,100.,100,-2.0,2.0,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0ArpH2);
    histos->UserHistogram("Event","Cent_v0CrpH2","VZERO-C RP;centrality (%);#Psi_{2}^{V0C} (rad.)",
    10,0.,100.,100,-2.0,2.0,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0CrpH2);

    histos->UserHistogram("Event","Cent_v0A0rpH2","VZERO-A RP;centrality (%);#Psi_{2}^{V0A0} (rad.)",
    10,0.,100.,100,-2.0,2.0,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0A0rpH2);
    histos->UserHistogram("Event","Cent_v0A3rpH2","VZERO-A RP;centrality (%);#Psi_{2}^{V0A3} (rad.)",
    10,0.,100.,100,-2.0,2.0,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0A3rpH2);
    histos->UserHistogram("Event","Cent_v0C0rpH2","VZERO-C RP;centrality (%);#Psi_{2}^{V0C0} (rad.)",
    10,0.,100.,100,-2.0,2.0,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0C0rpH2);
    histos->UserHistogram("Event","Cent_v0C3rpH2","VZERO-C RP;centrality (%);#Psi_{2}^{V0C3} (rad.)",
    10,0.,100.,100,-2.0,2.0,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0C3rpH2);
    } // hist: flow

    if(bHistFlowQA) {
    // TPC event plane
    histos->UserHistogram("Event","TPCxH2","TPC Qx component;TPCxH2",
    100,-1500.,1500.,
    AliDielectronVarManager::kTPCxH2);
    histos->UserHistogram("Event","TPCyH2","TPC Qy component;TPCyH2",
    100,-1500.,1500.,
    AliDielectronVarManager::kTPCyH2);
    histos->UserHistogram("Event","TPCrpH2","TPC reaction plane; #Psi^{TPC}",
    100,-2.,2.,
    AliDielectronVarManager::kTPCrpH2);
    histos->UserHistogram("Event","TPCsub1rpH2","TPC reaction plane sub1; #Psi^{sub1}",
    100,-2.,2.,
    AliDielectronVarManager::kTPCsub1rpH2);
    histos->UserHistogram("Event","TPCsub2rpH2","TPC reaction plane sub2; #Psi^{sub2}",
    100,-2.,2.,
    AliDielectronVarManager::kTPCsub2rpH2);
    histos->UserHistogram("Event","TPCsub12DiffH2","TPC reaction plane diff; cos(2(#Psi^{sub1}-#Psi^{sub2}))",
    100,-1.,1.,
    AliDielectronVarManager::kTPCsub12DiffH2);
    
  */

  /* // uncorrected eventplane
     histos->UserHistogram("Event","TPCxH2uc","TPC Qx component;TPCxH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCxH2uc);
     histos->UserHistogram("Event","TPCyH2uc","TPC Qy component;TPCyH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCyH2uc);
     histos->UserHistogram("Event","TPCrpH2uc","TPC reaction plane;TPCrpH2uc",
     100,-2.,2.,
     AliDielectronVarManager::kTPCrpH2uc);
     histos->UserHistogram("Event","TPCsub1xH2uc","TPC Qx component sub1;TPCsub1xH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCsub1xH2uc);
     histos->UserHistogram("Event","TPCsub1yH2uc","TPC Qy component sub1;TPCsub1yH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCsub1yH2uc);
     histos->UserHistogram("Event","TPCsub1rpH2uc","TPC reaction plane sub1;TPCsub1rpH2uc",
     100,-2.,2.,
     AliDielectronVarManager::kTPCsub1rpH2uc);
     histos->UserHistogram("Event","TPCsub2xH2uc","TPC Qx component sub2;TPCsub2xH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCsub2xH2uc);
     histos->UserHistogram("Event","TPCsub2yH2uc","TPC Qy component sub2;TPCsub2yH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCsub2yH2uc);
     histos->UserHistogram("Event","TPCsub2rpH2uc","TPC reaction plane sub2;TPCsub2rpH2uc",
     100,-2.,2.,
     AliDielectronVarManager::kTPCsub2rpH2uc);
     histos->UserHistogram("Event","TPCsub12DiffH2uc","TPC reaction plane difference;TPCsub12DiffH2uc",
     100,-1.,1.,
     AliDielectronVarManager::kTPCsub12DiffH2uc);
  */
  /*
  // EP resolution calculation
  histos->UserHistogram("Event","Cent_v0ATPCDiffH2","VZERO-A TPC diff;centrality (%);cos(2(#Psi^{V0A}-#Psi^{TPC}))",
  10,0.,100.,300,-1.0,1.0,
  AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0ATPCDiffH2);
  histos->UserHistogram("Event","Cent_v0CTPCDiffH2","VZERO-C TPC diff;centrality (%);cos(2(#Psi^{V0C}-#Psi^{TPC}))",
  10,0.,100.,300,-1.0,1.0,
  AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0CTPCDiffH2);
  histos->UserHistogram("Event","Cent_v0Av0CDiffH2","VZERO-A VZERO-C diff;centrality (%);cos(2(#Psi^{V0A}-#Psi^{V0C}))",
  10,0.,100.,300,-1.0,1.0,
  AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0Av0CDiffH2);
  histos->UserHistogram("Event","Cent_TPCsub12DiffH2","TPC-sub1 TPC-sub2 diff;centrality (%);cos(2(#Psi^{sub1}-#Psi^{sub2}))",
  10,0.,100.,300,-1.0,1.0,
  AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsub12DiffH2);
  // detector effects
  histos->UserHistogram("Event","Cent_TPCsub12DiffH2Sin","TPC-sub1 TPC-sub2 diff;centrality (%);sin(2(#Psi^{sub1}-#Psi^{sub2}))",
  10,0.,100.,300,-1.0,1.0,
  AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsub12DiffH2Sin);
  // recentering stuff
  histos->UserProfile("Pair","TPCxH2-Cent-RunNumber", ";centrality (%);run;#LTQ_{x}#GT",
  AliDielectronVarManager::kTPCxH2,
  AliDielectronHelper::MakeLinBinning(10, 0.,100.), GetRunNumbers(),
  AliDielectronVarManager::kCentrality, AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Pair","TPCyH2-Cent-RunNumber", ";centrality (%);run;#LTQ_{x}#GT",
  AliDielectronVarManager::kTPCyH2,
  AliDielectronHelper::MakeLinBinning(10, 0.,100.), GetRunNumbers(),
  AliDielectronVarManager::kCentrality, AliDielectronVarManager::kRunNumber);

  } //hist: flowQA


  // RUN QA
  if(bHistTrackQA) {

  // Event QA
  histos->UserHistogram("Event","RunNumber","Events per run;run;#events",
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Event","RunNumber-Centrality",";run;#LTcentrality#GT (%)", AliDielectronVarManager::kCentrality,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Event","RunNumber-CentralitySPD",";run;#LTcentrality_{SPD}#GT (%)", AliDielectronVarManager::kCentralitySPD,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Event","RunNumber-VtxZ",";run;#LTz_{vtx}#GT (cm)", AliDielectronVarManager::kZvPrim,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Event","RunNumber-sigVtxZ",";run;#LTRMS(z_{vtx})#GT (cm)", AliDielectronVarManager::kZvPrim,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber, "s;-10;10");

  //Track classes
  histos->SetReservedWords("Track");
  for (Int_t i=0; i<2; ++i) histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));

  histos->UserProfile("Track","RunNumber-Pt",";run;#LTp_{T}#GT (GeV/c)", AliDielectronVarManager::kPt,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Track","RunNumber-Eta",";run;#LT#eta#GT", AliDielectronVarManager::kEta,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Track","RunNumber-Phi",";run;#LT#varphi#GT", AliDielectronVarManager::kPhi,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Track","RunNumber-ImpactParXY",";run;#LTdXY#GT (cm)", AliDielectronVarManager::kImpactParXY,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Track","RunNumber-ImpactParZ",";run;#LTdZ#GT (cm)", AliDielectronVarManager::kImpactParZ,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  // TPC
  histos->UserProfile("Track","RunNumber-NclsTPC",";run;#LTN_{cls}^{TPC}#GT", AliDielectronVarManager::kNclsTPC,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Track","RunNumber_Cent-NclsTPC",";run;centrality (%);#LTN_{cls}^{TPC}#GT", AliDielectronVarManager::kNclsTPC,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(32, 0.,80.),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kCentrality);
  histos->UserProfile("Track","RunNumber_Eta-NclsTPC",";run;#eta;#LTN_{cls}^{TPC}#GT", AliDielectronVarManager::kNclsTPC,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(20, -1.,+1.),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kEta);
  histos->UserProfile("Track","RunNumber_Phi-NclsTPC",";run;#varphi;#LTN_{cls}^{TPC}#GT", AliDielectronVarManager::kNclsTPC,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(18, 0.,TMath::Pi()),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kPhi);

  // ITS
  histos->UserProfile("Track","RunNumber-NclsITS",";run;#LTN_{cls}^{ITS}#GT", AliDielectronVarManager::kNclsITS,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Track","RunNumber_Cent-NclsITS",";run;centrality (%);#LTN_{cls}^{ITS}#GT", AliDielectronVarManager::kNclsITS,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(32, 0.,80.),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kCentrality);
  histos->UserProfile("Track","RunNumber_Eta-NclsITS",";run;#eta;#LTN_{cls}^{ITS}#GT", AliDielectronVarManager::kNclsITS,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(20, -1.,+1.),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kEta);
  histos->UserProfile("Track","RunNumber_Phi-NclsITS",";run;#varphi;#LTN_{cls}^{ITS}#GT", AliDielectronVarManager::kNclsITS,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(18, 0.,TMath::Pi()),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kPhi);

  // TPC PID
  histos->UserProfile("Track","RunNumber-TPCnSigmaEle",";run;#LTn#sigma_{ele}^{TPC}#GT", AliDielectronVarManager::kTPCnSigmaEle,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Track","RunNumber_Cent-TPCnSigmaEle",";run;centrality (%);#LTn#sigma_{ele}^{TPC}#GT", AliDielectronVarManager::kTPCnSigmaEle,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(32, 0.,80.),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kCentrality);
  histos->UserProfile("Track","RunNumber_Eta-TPCnSigmaEle",";run;#eta;#LTn#sigma_{ele}^{TPC}#GT", AliDielectronVarManager::kTPCnSigmaEle,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(20, -1.,+1.),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kEta);
  histos->UserProfile("Track","RunNumber_Phi-TPCnSigmaEle",";run;#varphi;#LTn#sigma_{ele}^{TPC}#GT", AliDielectronVarManager::kTPCnSigmaEle,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(18, 0.,TMath::Pi()),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kPhi);
  // TOF PID
  histos->UserProfile("Track","RunNumber-TOFnSigmaEle",";run;#LTn#sigma_{ele}^{TOF}#GT", AliDielectronVarManager::kTOFnSigmaEle,
  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  histos->UserProfile("Track","RunNumber_Cent-TOFnSigmaEle",";run;centrality (%);#LTn#sigma_{ele}^{TOF}#GT", AliDielectronVarManager::kTOFnSigmaEle,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(32, 0.,80.),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kCentrality);
  histos->UserProfile("Track","RunNumber_Eta-TOFnSigmaEle",";run;#eta;#LTn#sigma_{ele}^{TOF}#GT", AliDielectronVarManager::kTOFnSigmaEle,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(20, -1.,+1.),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kEta);
  histos->UserProfile("Track","RunNumber_Phi-TOFnSigmaEle",";run;#varphi;#LTn#sigma_{ele}^{TOF}#GT", AliDielectronVarManager::kTOFnSigmaEle,
  GetRunNumbers(), AliDielectronHelper::MakeLinBinning(18, 0.,TMath::Pi()),
  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kPhi);
  }


  if(bHistPair) {
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");

  //Pair classes
  // to fill also mixed event histograms loop until 7 or 10
  for (Int_t i=0; i<(bHistPairME ? 8 : 3); ++i){
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
  }
  //legs from pair (fill SE)
  for (Int_t i=0; i<3; ++i){
  histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
  }
  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
  histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
  //track rotation
  //   histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
  //   histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));

  // Vertex
  histos->UserHistogram("Track","ImpactParXY",";dXY (cm);#tracks", 500,-1.,1.,
  AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","ImpactParZ",";dZ (cm);#tracks", 600,-3.,3.,
  AliDielectronVarManager::kImpactParZ);

  // Kinematics
  histos->UserHistogram("Track","pt",";p_{T} (GeV/c);#tracks", 400,0,20.,
  AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Eta_Phi",";#eta;#varphi;#tracks", 200,-1,1,200,0,6.285,
  AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  // TPC
  histos->UserHistogram("Track","TPCnCls",";N_{cls}^{TPC};#tracks", 160,-0.5,159.5,
  AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","TPCsignalN",";N_{cls}^{TPC};#tracks", 160,-0.5,159.5,
  AliDielectronVarManager::kTPCsignalN);
  histos->UserHistogram("Track","NclFr",";N_{max.cls}^{TPC};#tracks", 160,-0.5,159.5,
  AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","Ncl_NclFr",";N_{cls}^{TPC};N_{max.cls}^{TPC};#tracks",
  160,-0.5,159.5,160,-0.5,159.5,
  AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);

  // TRD
  histos->UserHistogram("Track","NtrklTRD",";N_{trkl}^{TRD};#tracks",
  8,-0.5,7.5,
  AliDielectronVarManager::kTRDpidQuality);

  // PID
  histos->UserHistogram("Track","dEdx_P",";p (GeV/c);TPC signal (arb units);#tracks",
  400,0.2,20.,200,0.,200.,
  AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaEle_P","p (GeV/c);n#sigma_{ele}^{TPC};#tracks",
  400,0.2,20.,200,-10.,10.,
  AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TOFbeta_P",";p (GeV/c);#beta;#tracks",
  250,0.0,5.0,300,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta,kTRUE);

  ///// add histograms to Pair classes /////
  histos->UserHistogram("Pair","InvMass",";m_{ee} (GeV/c^{2});#pairs",
  300,.0,300*0.04, AliDielectronVarManager::kM); // 40MeV bins, 12GeV/c2
  histos->UserHistogram("Pair","Rapidity",";y;#pairs",
  100,-1.,1., AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","Pt",";p_{T} (GeV/c);#pairs",
  400,0,20., AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","OpeningAngle",";opening angle (rad.);#pairs",
  100,0.,3.15, AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","Chi2NDF",";#chi^{2}/NDF;#pairs",
  100,0.,20, AliDielectronVarManager::kChi2NDF);
  histos->UserHistogram("Pair","PsiPair",";#psi;#pairs",
  100,0.,3.15, AliDielectronVarManager::kPsiPair);
  histos->UserHistogram("Pair","R",";r (cm)",
  200,0.,100., AliDielectronVarManager::kR);
  histos->UserHistogram("Pair","LegDist",";dca (cm)",
  50,0.,5., AliDielectronVarManager::kLegDist);
  histos->UserHistogram("Pair","LegDistXY",";dca_{xy} (cm)",
  50,0.,5., AliDielectronVarManager::kLegDistXY);

  //// FLOW results use tprofiles
  if(bHistFlow) {
  histos->UserProfile("Pair","M_Cent_Pt_v0ACrpH2FlowV2",
  "cos(2(#varphi-#Psi^{V0AC}));mass (GeV/c^{2});centrality (%);p_{T} (GeV/c)",
  AliDielectronVarManager::kv0ACrpH2FlowV2,
  125,0.,125*.04, 10, 0.,100., 200,0.,100.,
  AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPt);
  histos->UserProfile("Pair","M_Cent_Pt_v0ArpH2FlowV2",
  "cos(2(#varphi-#Psi^{V0A}));mass (GeV/c^{2});centrality (%);p_{T} (GeV/c)",
  AliDielectronVarManager::kv0ArpH2FlowV2,
  125,0.,125*.04, 10, 0.,100., 200,0.,100.,
  AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPt);
  histos->UserProfile("Pair","M_Cent_Pt_v0CrpH2FlowV2",
  "cos(2(#varphi-#Psi^{V0C}));mass (GeV/c^{2});centrality (%);p_{T} (GeV/c)",
  AliDielectronVarManager::kv0CrpH2FlowV2,
  125,0.,125*.04, 10, 0.,100., 200,0.,100.,
  AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPt);
  } //hist: flow
  } //hist: pair
  */

  ////// MONTE CARLO //////
  /*
    if(cutDefinition == kTOFTRD && hasMC) {
    histos->AddClass("MCEvent");
    histos->UserHistogram("MCEvent","Cent_NJPsis","Centrality vs. generated incl. J/#psi per event;centrality (%);N_{J/#psi}",
    10,0.,100., 21,-0.5,20.5,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kNumberOfJPsis);
    }
  */

  //  histos->UserHistogram("Track","TRDprobEle",";P_{ele}^{TRD};#tracks", 
  //			100,0.,1.,            AliDielectronVarManager::kTRDprobEle);
  
  die->SetHistogramManager(histos);
}

void InitHF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the HF arrays
  //

  AliDielectronHF *hf=new AliDielectronHF(die->GetName(),die->GetTitle());
  //  if(hasMC) hf->SetStepForMCGenerated();
  hf->SetPairTypes(AliDielectronHF::kAll);
  hf->SetVariable(AliDielectronVarManager::kM, 125, 0.0, 0.04*125);

  hf->AddCutVariable(AliDielectronVarManager::kCentrality,  "0.,5.,10.,20.,40.,50.,60.,80."  );
  hf->AddCutVariable(AliDielectronVarManager::kPt,          "0.,2.5,5.,100."                 );
  hf->AddCutVariable(AliDielectronVarManager::kDeltaPhiv0ArpH2, 8,-1.*TMath::Pi(),TMath::Pi());
  //  hf->AddCutVariable(AliDielectronVarManager::kY,           1, -0.9, 0.9                     );
  //  hf->AddCutVariable(AliDielectronVarManager::kPt,          "0.8, 1.0, 1.1, 1.2, 1.5, 100.0", kTRUE, AliDielectronHF::kBinToMax);
  // hf->AddCutVariable(AliDielectronVarManager::kNclsTPC,     "70,90,100,120,160",              kTRUE, AliDielectronHF::kBinToMax);
  //  hf->AddCutVariable(AliDielectronVarManager::kTPCnSigmaEle,"-4,-3,-2.5,-2,2,2.5,3,4",        kTRUE, AliDielectronHF::kSymBin);
  //hf->AddCutVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.,3.5,4.,100.",                 kTRUE, AliDielectronHF::kBinToMax);
  //hf->AddCutVariable(AliDielectronVarManager::kITSLayerFirstCls,4,0.,4.,              kFALSE, kTRUE, AliDielectronHF::kBinFromMin);
  //hf->AddCutVariable(AliDielectronVarManager::kNclsITS,         5,2.,7.,              kFALSE, kTRUE, AliDielectronHF::kBinToMax);
  //hf->AddCutVariable(AliDielectronVarManager::kRunNumber,  GetRunNumbers());

  die->SetHistogramArray(hf);
}

void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the CF Manager if needed
  //

  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  // event variables
  cf->AddVariable(AliDielectronVarManager::kCentrality,"0.,5.,10.,20.,40.,50.,60.,80.");
  //    if(!hasMC) cf->AddVariable(AliDielectronVarManager::kZvPrim,20, -10., 10.);
  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kNacc,20,0.,3000.0);
  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kNVtxContrib,20,0.,4000.);
  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kRunNumber, GetRunNumbers() );

  // pair variables
  //    cf->AddVariable(AliDielectronVarManager::kY,"-0.8,0.8");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);
  cf->AddVariable(AliDielectronVarManager::kPt,"0., 1., 2.5, 5., 100.0");
  if(hasMC) cf->AddVariable(AliDielectronVarManager::kTRDpidEffPair,101,0.0,1.01);
  //    if(hasMC) cf->AddVariable(AliDielectronVarManager::kThetaCS,15,-1.,1.);

  // flow variables
  cf->AddVariable(AliDielectronVarManager::kDeltaPhiv0ArpH2,4,-1.*TMath::Pi(),TMath::Pi());
  cf->AddVariable(AliDielectronVarManager::kDeltaPhiv0CrpH2,4,-1.*TMath::Pi(),TMath::Pi());

  // leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.8, 1.0, 1.1, 1.2, 1.5, 100.0",kTRUE);
  if(hasMC) cf->AddVariable(AliDielectronVarManager::kEta,"-0.9,0.9",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,7,-1.5,5.5,kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kNclsITS,"1,2,3,4,5,6",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-3,-2.5,-2,2,2.5,3",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"2.5,3.0,3.5,4.0,4.5,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"70, 90, 100, 120, 160",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.5,4.0,4.5,5.0,100",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kTOFnSigmaEle,"-3,-2,2,3",kTRUE); break;
  //    cf->AddVariable(AliDielectronVarManager::kTRDpidQuality,"3.5, 4.5, 5.5, 6.5",kTRUE);
  //    if(!hasMC && isESD) cf->AddVariable(AliDielectronVarManager::kTRDchi2,"-1.,0.,2.,4.",kTRUE);

  // mc steps
  if(hasMC) {
    if(cutDefinition==kTOFTRD) cf->SetStepForMCtruth();
    cf->SetStepsForMCtruthOnly();
    // cf->SetStepsForBackground();
  }

  die->SetCFManagerPair(cf);
}

void AddMCSignals(AliDielectron *die){
  //Do we have an MC handler?
  if (!hasMC) return;
  
  AliDielectronSignalMC* inclusiveJpsi = new AliDielectronSignalMC("inclusiveJpsi","Inclusive J/psi");
  inclusiveJpsi->SetLegPDGs(11,-11);
  inclusiveJpsi->SetMotherPDGs(443,443);
  inclusiveJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  inclusiveJpsi->SetFillPureMCStep(kTRUE);
  inclusiveJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  inclusiveJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(inclusiveJpsi);
  
  AliDielectronSignalMC* promptJpsi = new AliDielectronSignalMC("promptJpsi","Prompt J/psi");   // prompt J/psi (not from beauty decays)
  promptJpsi->SetLegPDGs(11,-11);
  promptJpsi->SetMotherPDGs(443,443);
  promptJpsi->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsi->SetFillPureMCStep(kTRUE);
  promptJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(promptJpsi);
  
  AliDielectronSignalMC* beautyJpsi = new AliDielectronSignalMC("beautyJpsi","Beauty J/psi");
  beautyJpsi->SetLegPDGs(11,-11);
  beautyJpsi->SetMotherPDGs(443,443);
  beautyJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  beautyJpsi->SetGrandMotherPDGs(500,500);
  beautyJpsi->SetFillPureMCStep(kTRUE);
  beautyJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(beautyJpsi);
  
  AliDielectronSignalMC* directJpsi = new AliDielectronSignalMC("directJpsi","Direct J/psi");   // embedded J/psi
  directJpsi->SetLegPDGs(11,-11);
  directJpsi->SetMotherPDGs(443,443);
  directJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  directJpsi->SetFillPureMCStep(kTRUE);
  directJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  directJpsi->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);
  directJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  directJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(directJpsi);
  
  AliDielectronSignalMC* conversionElePairs = new AliDielectronSignalMC("conversionElePairs","conversion electron pairs");  // pairs made from conversion (may be also from 2 different conversions)                                                                                                                                                                                             
  conversionElePairs->SetLegPDGs(11,-11);
  conversionElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  conversionElePairs->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  conversionElePairs->SetMotherPDGs(22,22);
  //   die->AddSignalMC(conversionElePairs);
}

void SetEtaCorrection()
{
  if (AliDielectronPID::GetEtaCorrFunction()) return;
  
  TString etaMap="$TRAIN_ROOT/jpsi_JPSI/EtaCorrMaps.root";
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) etaMap="$ALICE_ROOT/PWGDQ/dielectron/files/EtaCorrMaps.root";
  if (gSystem->AccessPathName(gSystem->ExpandPathName(etaMap.Data()))){
    Error("ConfigPbPb","Eta map not found: %s",etaMap.Data());
    return;
  }

  TFile f(etaMap.Data());
  if (!f.IsOpen()) return;
  TList *keys=f.GetListOfKeys();

  for (Int_t i=0; i<keys->GetEntries(); ++i){
    TString kName=keys->At(i)->GetName();
    TPRegexp reg(kName);
    if (reg.MatchB(list)){
      printf(" Using Eta Correction Function: %s\n",kName.Data());
      AliDielectronPID::SetEtaCorrFunction((TF1*)f.Get(kName.Data()));
    }
  }
}

TVectorD *GetRunNumbers() {
  
  Double_t runLHC10h[] = { // all good runs based on RCT 29.Mai
    139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161, 137135
  };
  
  Double_t runLHC11h[] = { // all good runs based on RCT 29.Mai
    170593, 170572, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170027, 169965, 169923, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169591, 169590, 169588, 169587, 169586, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168311, 168310, 168115, 168108, 168107, 168105, 168076, 168069, 167988, 167987, 167985, 167920, 167915
  };
  
  // selection via environement variable (works only for gsi trains)

  
  if(list.Contains("LHC10h") || list.Contains("LHC11a10")) {
    Int_t size = (int) (sizeof(runLHC10h)/sizeof(Double_t));
    TVectorD *vec = new TVectorD(size+1);
    
    (*vec)[size] = runLHC10h[0] + 1;
    for (int i = 0; i < size; i++) {
      (*vec)[i] = runLHC10h[size-1-i];
    }
    //    vec->Print("");
    return vec;
  }

  if( list.Contains("LHC11h") || list.Contains("LHC12a17") ) {
    
    Int_t size = (int) (sizeof(runLHC11h)/sizeof(Double_t));
    TVectorD *vec = new TVectorD(size+1);
    
    (*vec)[size] = runLHC11h[0] + 1;
    for (int i = 0; i < size; i++) {
      (*vec)[i] = runLHC11h[size-1-i];
    }
    //   vec->Print("");
    return vec;
  }

  TVectorD *vec = new TVectorD(2);
  (*vec)[0] = 0;
  (*vec)[0] = 1;
  return vec;
     
}
