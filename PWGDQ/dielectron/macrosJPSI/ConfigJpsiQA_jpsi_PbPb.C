void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupEventCuts(AliDielectron *die, ULong64_t triggers);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die,  Int_t cutDefinition);

void AddMCSignals(AliDielectron *die);

void SetEtaCorrection();
TVectorD *GetRunNumbers();

TString names=("QA");
enum { kQA=0 };

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t  isESD      = kTRUE;
Bool_t  hasMC      = kFALSE;
TString list       = gSystem->Getenv("LIST");

AliDielectron* ConfigJpsiQA_jpsi_PbPb(Int_t cutDefinition, TString prod="", ULong64_t triggers=AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB)
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
  printf(" Adding %s %s config %s for %s \n",(isESD?"ESD":"AOD"),(hasMC?"MC":""),name.Data(),list.Data());
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("Track cuts: %s",name.Data()));
  die->SetHasMC(hasMC);

  // cut setup
  SetupEventCuts(die,triggers);
  SetupTrackCuts(die,cutDefinition);
  SetupPairCuts(die,cutDefinition);

  // MC signals
  if(hasMC) {
    AddMCSignals(die);
    printf(" Add %d MC signals \n",die->GetMCSignals()->GetEntriesFast());
  }
  // histogram setup
  InitHistograms(die,cutDefinition);
  printf(" Add %d classes to the manager \n",die->GetHistogramList()->GetEntries());

  // CF container setup
  InitCF(die,cutDefinition);

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
  //  die->SetNoPairing();
  die->SetPreFilterUnlikeOnly();
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
void SetupEventCuts(AliDielectron *die, ULong64_t triggers)
{
  //
  // Setup the event cuts
  //

  // trigger specific centrality cuts (reject trigger inefficiencies)
  Double_t minCent=0.0, maxCent=80.;
  if(!hasMC) {
    switch(triggers) {
    case AliVEvent::kCentral:     minCent= 0.; maxCent= 9.; break;
    case AliVEvent::kSemiCentral: minCent=12.; maxCent=53.; break;
    case AliVEvent::kMB:          minCent= 0.; maxCent=80.; break;
    default:                      minCent= 0.; maxCent=80.; break;
    }
  }

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","eventCuts");
  if(!isESD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
  //eventCuts->SetCentralityRange(minCent,maxCent);

  /*
  TF1 *fMean  = new TF1("fMean", "pol1",               0,25e+3);
  fMean->SetParameters(631.301, 1.49836);
  TF1 *fSigma = new TF1("fSigma","[0]+sqrt([1]*x+[2])",0,25e+3);
  fSigma->SetParameters(-25.7843, 32.8055, 35275.7);
  eventCuts->SetCutOnV0MultipicityNTrks(fMean, fSigma, 4.0);
  */

  eventCuts->Print();
  die->GetEventFilter().AddCuts(eventCuts);

}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);

  //  AliDielectronVarCuts *mc = new AliDielectronVarCuts("MCCut","MCCut");
  //mc->AddCut(AliDielectronVarManager::kHasCocktailMother,      0); // exclude enhanced signals
  //mc->AddCut(AliDielectronVarManager::kHasCocktailGrandMother, 0); // exclude enhanced signals
  //  mc->AddCut(AliDielectronVarManager::kPdgCode, -6,  6, kTRUE); // exclude quarks to speed up
  //  mc->AddCut(AliDielectronVarManager::kPdgCode, 21, 21, kTRUE); // exclude gluons to speed up
  //if(hasMC) {
    //    cuts->AddCut(mc);
    //    mc->Print();
  //}

  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  //  trkFilter->SetMinNCrossedRowsOverFindable(0.6);
  //if(!isESD) cuts->AddCut(trkFilter);

  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 15); // ITS-4 = 1+2+4+8
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);
  cuts->AddCut(trkCuts);

  //Pt cut, should make execution a bit faster
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("PtCut","PtCut");
  pt->AddCut(AliDielectronVarManager::kPt,0.8,1e30);    //1.1
  cuts->AddCut(pt);
  pt->Print();

  // track cuts ESD and AOD
  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");
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

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronPID *pid = new AliDielectronPID("PID","PID");
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

    // TPC
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,   -4.,4.0);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -100.,4.0,0.,0.,kTRUE);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -100.,3.5,0.,0.,kTRUE);

    // TOF
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-5,5.,0.,0.,kFALSE,AliDielectronPID::kIfAvailable);

  } //hasMC

  //  cuts->AddCut(pid);
  /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ PID CUTS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/


  // exclude conversion electrons selected by the tender
  AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  //cuts->AddCut(noconv);

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
  case kQA:    gCut = 0.05;  break;
  default:     gCut = 0.05;
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
  //AliDielectronVarCuts *rapCut=new AliDielectronVarCuts("|Y|<.9","|Y|<.9");
  //rapCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  //die->GetPairFilter().AddCuts(rapCut);

  // minv cut (for better jpsi candidate to mc comparison) 
  AliDielectronVarCuts *minvCut=new AliDielectronVarCuts("PairCut","PairCut");
  minvCut->AddCut(AliDielectronVarManager::kM,2.92,3,.16);
  if(!hasMC) {
    //    minvCut->Print();
    //    die->GetPairFilter().AddCuts(minvCut);
  }
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
  case kQA:       bHistEvtQA=kTRUE; bHistTrackQA=kTRUE; bHistPairQA=kTRUE; break;
  }


  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());


  ////// EVENT HISTOS /////
  if(bHistEvtQA) {

    //add histograms to event class
    histos->AddClass("Event");
    Int_t bins[]={125,100,100,55}; Double_t min[]={0.,0.,0.,0.}; Double_t max[]={25000.,20000.,4000.,1.1}; 
    UInt_t var[]={AliDielectronVarManager::kMultV0,AliDielectronVarManager::kNTrk,AliDielectronVarManager::kNacc,AliDielectronVarManager::kMatchEffITSTPC};
    histos->UserSparse("Event", 4, bins, min, max, var);
    //histos->UserHistogram("Event", 3, bins, min, max, var);
    histos->UserHistogram("Event","","", 80,0.,80., AliDielectronVarManager::kCentrality);
    histos->UserHistogram("Event","","", GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Event","","", GetRunNumbers(), AliDielectronHelper::MakeLinBinning(80,0.,80.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kCentrality);
    histos->UserProfile("Event","","", AliDielectronVarManager::kCentrality,    GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","","", AliDielectronVarManager::kCentralitySPD, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","","", AliDielectronVarManager::kZvPrim,        GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","","", AliDielectronVarManager::kZvPrim,	GetRunNumbers(), AliDielectronVarManager::kRunNumber, "s;-10;10");
    histos->UserProfile("Event","","", AliDielectronVarManager::kZvPrim,	AliDielectronHelper::MakeLinBinning(80,0.,80.), AliDielectronVarManager::kCentrality);
    histos->UserProfile("Event","","", AliDielectronVarManager::kZvPrim,        250,0.,25000.,   AliDielectronVarManager::kMultV0A);
    histos->UserProfile("Event","","", AliDielectronVarManager::kZvPrim,        250,0.,25000.,   AliDielectronVarManager::kMultV0C);
    histos->UserHistogram("Event","","", GetRunNumbers(), AliDielectronHelper::MakeLinBinning(150,-15.,15.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","","", 80.,0.,80., 150,-15.,15.,
			  AliDielectronVarManager::kCentrality, AliDielectronVarManager::kZvPrim);

    histos->UserHistogram("Event","","", 300.,0.,6000., 300,0.,6000.,
			  AliDielectronVarManager::kNVtxContrib, AliDielectronVarManager::kNVtxContribTPC);
    histos->UserHistogram("Event","","", 200.,0.,4000., 200,0.,4000.,
			  AliDielectronVarManager::kNVtxContrib, AliDielectronVarManager::kNacc);


    histos->UserHistogram("Event","","", GetRunNumbers(), AliDielectronHelper::MakeLinBinning(250,0.,25000.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","","", 80.,0.,80., 250,0.,25000.,
			  AliDielectronVarManager::kCentrality, AliDielectronVarManager::kMultV0);
    histos->UserProfile("Event","","", AliDielectronVarManager::kMultV0,        80.,0.,80., AliDielectronVarManager::kCentrality);
    histos->UserProfile("Event","","", AliDielectronVarManager::kMultV0,        80.,0.,80., AliDielectronVarManager::kCentralitySPD);
    histos->UserProfile("Event","","", AliDielectronVarManager::kMultV0,        GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Event","","", 250,0.,25000., AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","","", 250,0.,25000., AliDielectronVarManager::kMultV0A);
    histos->UserHistogram("Event","","", 250,0.,25000., AliDielectronVarManager::kMultV0C);
    histos->UserHistogram("Event","","", 200,0.,20000., 250,0.,25000., AliDielectronVarManager::kNTrk, AliDielectronVarManager::kMultV0 );
    histos->UserHistogram("Event","","", 200,0.,20000., 250,0.,25000., AliDielectronVarManager::kNTrk, AliDielectronVarManager::kMultV0A);
    histos->UserHistogram("Event","","", 200,0.,20000., 250,0.,25000., AliDielectronVarManager::kNTrk, AliDielectronVarManager::kMultV0C);

    

    histos->UserHistogram("Event","","", 
			  AliDielectronHelper::MakeLinBinning(200,0.,20000.),
			  AliDielectronHelper::MakeLinBinning(250,0.,25000.),
			  GetRunNumbers(), 
			  AliDielectronVarManager::kNTrk, AliDielectronVarManager::kMultV0A,  AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Event","","", 
			  AliDielectronHelper::MakeLinBinning(200,0.,20000.),
			  AliDielectronHelper::MakeLinBinning(250,0.,25000.),
			  GetRunNumbers(), 
			  AliDielectronVarManager::kNTrk, AliDielectronVarManager::kMultV0C,  AliDielectronVarManager::kRunNumber);

    histos->UserHistogram("Event","","", 110,0.,1.1, AliDielectronVarManager::kMatchEffITSTPC);
    histos->UserHistogram("Event","","", 200,0.,20000., 110,0.,1.1,
			  AliDielectronVarManager::kNTrk, AliDielectronVarManager::kMatchEffITSTPC);
    histos->UserHistogram("Event","","", 80.,0.,80., 80.,0.,80., AliDielectronVarManager::kCentrality, AliDielectronVarManager::kCentralitySPD);


    histos->UserHistogram("Event","","", GetRunNumbers(), AliDielectronHelper::MakeLinBinning(100,-2.,+2.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kv0ArpH2);
    histos->UserHistogram("Event","","", GetRunNumbers(), AliDielectronHelper::MakeLinBinning(100,-2.,+2.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kv0CrpH2);
    histos->UserHistogram("Event","","", GetRunNumbers(), AliDielectronHelper::MakeLinBinning(100,-2.,+2.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kTPCrpH2uc);
    histos->UserHistogram("Event","","", 80,0.,80., 100,-2.,+2.,
			  AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0CrpH2);
    histos->UserHistogram("Event","","", 80,0.,80., 100,-2.,+2.,
			  AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0ArpH2);
    histos->UserHistogram("Event","","", 80,0.,80., 100,-2.,+2.,
			  AliDielectronVarManager::kCentrality, AliDielectronVarManager::kTPCrpH2uc);
    histos->UserHistogram("Event","","", 80,0.,80., 100,0.,250.,
			  AliDielectronVarManager::kCentrality, AliDielectronVarManager::kTPCmagH2uc);

  }

  ////// PAIR HISTOS /////
  if(bHistPairQA) {

    //add histograms to track class
    histos->SetReservedWords("Pair");
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PM)));

    //add MC signal histograms to track class
    if(die->GetMCSignals()) {
      for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i)
	histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
    }

    histos->UserHistogram("Pair","","",  20, 0.,10.,           AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","","", 200,-1.,+1.,           AliDielectronVarManager::kEta);
    histos->UserHistogram("Pair","","", 180,-1.*TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kPhi);
    histos->UserHistogram("Pair","","", 180, 0.,TMath::Pi(),   AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","","", 300, 0.,300*0.04,      AliDielectronVarManager::kM);
    histos->UserHistogram("Pair","","", 180, 0.,TMath::Pi(),   AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","","", 200, 0.,1.,            AliDielectronVarManager::kCosPointingAngle);
    histos->UserHistogram("Pair","","", 200, 0.,20.,           AliDielectronVarManager::kR);
    histos->UserHistogram("Pair","","", 100,-1.,+1.,           AliDielectronVarManager::kThetaCS);

    if(hasMC) histos->UserHistogram("Pair","","",10000,-5000.5,4999.5, 300, 0.,300*0.04, AliDielectronVarManager::kPdgCode, AliDielectronVarManager::kM);

    histos->UserProfile("Pair","","",  AliDielectronVarManager::kPt,  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Pair","","",  AliDielectronVarManager::kEta, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Pair","","",  AliDielectronVarManager::kPhi, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Pair","","",  AliDielectronVarManager::kOpeningAngle, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Pair","","",  AliDielectronVarManager::kM,   GetRunNumbers(), AliDielectronVarManager::kRunNumber);
  }

  ////// TRACK HISTOS /////
  if(bHistTrackQA) {

    //add histograms to track class
    histos->SetReservedWords("Track");
    for (Int_t i=0; i<2; ++i) histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
    //legs from pair (fill SE PM)
    for (Int_t i=1; i<2; ++i) histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
    
    //add MC signal histograms to track class
    if(die->GetMCSignals()) {
      for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i)
	histos->AddClass(Form("Track_Legs_%s",die->GetMCSignals()->At(i)->GetName()));
    }

    histos->UserHistogram("Track","","", 400, 0.,20., AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","","", 200,-1.,+1., AliDielectronVarManager::kEta);
    histos->UserHistogram("Track","","", 180,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    histos->UserProfile("Track","","", AliDielectronVarManager::kPt,  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","","", AliDielectronVarManager::kEta, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","","", AliDielectronVarManager::kPhi, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","","", AliDielectronVarManager::kPt,  80,0.,80., AliDielectronVarManager::kCentrality);
    histos->UserProfile("Track","","", AliDielectronVarManager::kEta, 80,0.,80., AliDielectronVarManager::kCentrality);
    histos->UserProfile("Track","","", AliDielectronVarManager::kPhi, 80,0.,80., AliDielectronVarManager::kCentrality);

    histos->UserHistogram("Track","","", 400,-1.,+1., AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","","", 600,-3.,+3., AliDielectronVarManager::kImpactParZ);
    histos->UserProfile("Track","","", AliDielectronVarManager::kImpactParXY,
			200,-1.,+1., 180,0.,TMath::TwoPi(), AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    histos->UserProfile("Track","","", AliDielectronVarManager::kImpactParZ,
			200,-1.,+1., 180,0.,TMath::TwoPi(), AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    histos->UserProfile("Track","","", AliDielectronVarManager::kImpactParXY, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","","", AliDielectronVarManager::kImpactParZ,  GetRunNumbers(), AliDielectronVarManager::kRunNumber);


    // TPC
    histos->UserHistogram("Track","","", 160,0.,160., AliDielectronVarManager::kNclsTPC);
    histos->UserProfile("Track","","", AliDielectronVarManager::kNclsTPC,  GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","","", AliDielectronVarManager::kNclsTPC,  80, 0.,80.,      AliDielectronVarManager::kCentrality);
    histos->UserProfile("Track","","", AliDielectronVarManager::kNclsTPC,  40, 0.,20.,      AliDielectronVarManager::kPt);
    histos->UserProfile("Track","","", AliDielectronVarManager::kNclsTPC,
			200,-1.,+1., 180,0.,TMath::TwoPi(), AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);


    histos->UserHistogram("Track","","", 8,0.,8., AliDielectronVarManager::kTPCclsSegments);
    histos->UserProfile("Track","","", AliDielectronVarManager::kTPCclsSegments, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","","", AliDielectronVarManager::kTPCclsSegments, 80, 0.,80.,      AliDielectronVarManager::kCentrality);
    histos->UserProfile("Track","", "", AliDielectronVarManager::kTPCclsSegments,
			200,-1.,+1., 180,0.,TMath::TwoPi(), AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);

    histos->UserHistogram("Track","","", 110,0.,1.1, AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserProfile("Track","","", AliDielectronVarManager::kNFclsTPCfCross, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","","", AliDielectronVarManager::kNFclsTPCfCross, 80, 0.,80.,      AliDielectronVarManager::kCentrality);
    histos->UserProfile("Track","", "", AliDielectronVarManager::kNFclsTPCfCross,
			200,-1.,+1., 180,0.,TMath::TwoPi(), AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);

    // ITS
    histos->UserHistogram("Track","","", 7,0.,7., AliDielectronVarManager::kNclsITS);
    histos->UserProfile("Track","","", AliDielectronVarManager::kNclsITS, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Track","","", AliDielectronVarManager::kNclsITS, 80, 0.,80.,      AliDielectronVarManager::kCentrality);
    histos->UserProfile("Track","","", AliDielectronVarManager::kNclsITS, 40, 0.,20.,      AliDielectronVarManager::kPt);
    histos->UserProfile("Track","", "", AliDielectronVarManager::kNclsITS,
			200,-1.,+1., 180,0.,TMath::TwoPi(), AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);

    histos->UserHistogram("Track","","", 7,0.,7., AliDielectronVarManager::kITSLayerFirstCls);
    histos->UserProfile("Track","", "", AliDielectronVarManager::kITSLayerFirstCls,
			200,-1.,+1., 180,0.,TMath::TwoPi(), AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);

    // TRD
    histos->UserHistogram("Track","", "",   7, 0., 7., AliDielectronVarManager::kTRDpidQuality);
    histos->UserHistogram("Track","", "", 105,-1.,20., AliDielectronVarManager::kTRDchi2);

    // TPC PID
    histos->UserProfile("Track","","", AliDielectronVarManager::kTPCnSigmaEle, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","","", GetRunNumbers(), AliDielectronHelper::MakeLinBinning(16, 0.,80.), AliDielectronHelper::MakeLinBinning(100, -5.,+5),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","","", GetRunNumbers(), AliDielectronHelper::MakeLinBinning(20, -1.,1.), AliDielectronHelper::MakeLinBinning(100, -5.,+5),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","","", GetRunNumbers(), AliDielectronHelper::MakeLinBinning(18, 0.,TMath::TwoPi()), AliDielectronHelper::MakeLinBinning(100, -5.,+5),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kPhi, AliDielectronVarManager::kTPCnSigmaEle);

    // TOF PID
    histos->UserHistogram("Track","", "", 100,-10.,+10., AliDielectronVarManager::kTOFnSigmaEle);
    histos->UserHistogram("Track","", "", 250,0.0,5., 100,-10.,+10., AliDielectronVarManager::kPIn, AliDielectronVarManager::kTOFnSigmaEle);
    histos->UserProfile("Track","","", AliDielectronVarManager::kTOFnSigmaEle, GetRunNumbers(), AliDielectronVarManager::kRunNumber, "h;-10;+10");
    histos->UserProfile("Track","","", AliDielectronVarManager::kTOFnSigmaEle, 
			GetRunNumbers(), AliDielectronHelper::MakeLinBinning(16, 0.,80.), AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kCentrality,"h;-10;+10");
    histos->UserProfile("Track","","", AliDielectronVarManager::kTOFnSigmaEle,
			AliDielectronHelper::MakeLinBinning(20, -1.,+1.), AliDielectronHelper::MakeLinBinning(18, 0.,TMath::TwoPi()),
			AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi, "h;-10;+10");

    // MC pdg codes
    if(hasMC) {
      //      histos->UserHistogram("Event","","",   25,    0.,   25.,  AliDielectronVarManager::kNumberOfJPsisPrompt);
      histos->UserHistogram("Track","","",10000,-5000.5,4999.5, AliDielectronVarManager::kPdgCodeGrandMother);
      histos->UserHistogram("Track","","",10000,-5000.5,4999.5, AliDielectronVarManager::kPdgCodeMother);
      histos->UserHistogram("Track","","",10000,-5000.5,4999.5, AliDielectronVarManager::kPdgCode);
    }

  }

  die->SetHistogramManager(histos);
}

void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the CF Manager
  //

  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  // event variables
  cf->AddVariable(AliDielectronVarManager::kMultV0,         125,0.,25000. );
  cf->AddVariable(AliDielectronVarManager::kNTrk,           100,0.,20000. );
  cf->AddVariable(AliDielectronVarManager::kNacc,           100,0., 4000. );
  cf->AddVariable(AliDielectronVarManager::kMatchEffITSTPC,  55,0.,    1.1);

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

  AliDielectronSignalMC* gammaConversion = new AliDielectronSignalMC("gammaConversion","gamma conversions");
  gammaConversion->SetLegPDGs(11,-11);
  gammaConversion->SetCheckBothChargesLegs(kTRUE,kTRUE);
  gammaConversion->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  gammaConversion->SetMotherPDGs(22,22);
  gammaConversion->SetMothersRelation(AliDielectronSignalMC::kSame);
  die->AddSignalMC(gammaConversion);

  /*
  AliDielectronSignalMC* gammaConversionHI = new AliDielectronSignalMC("gammaConversionHI","gamma conversions Hijing");
  gammaConversionHI->SetLegPDGs(11,-11);
  gammaConversionHI->SetCheckBothChargesLegs(kTRUE,kTRUE);
  gammaConversionHI->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  gammaConversionHI->SetMotherSources(AliDielectronSignalMC::kNoCocktail, AliDielectronSignalMC::kNoCocktail);
  gammaConversionHI->SetGrandMotherSources(AliDielectronSignalMC::kNoCocktail, AliDielectronSignalMC::kNoCocktail);
  gammaConversionHI->SetMotherPDGs(22,22);
  gammaConversionHI->SetMothersRelation(AliDielectronSignalMC::kSame);
  die->AddSignalMC(gammaConversionHI);
  */

  AliDielectronSignalMC* conversionElePairs = new AliDielectronSignalMC("conversionEle","conversion electron pairs");  // pairs made from conversion (may be also from 2 different conversions)
  conversionElePairs->SetLegPDGs(11,-11);
  conversionElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  conversionElePairs->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  conversionElePairs->SetMotherPDGs(22,22);
  die->AddSignalMC(conversionElePairs);

  AliDielectronSignalMC* inclJpsiElePairs = new AliDielectronSignalMC("inclJpsiEle","Inclusive J/psi electron pairs"); // pairs made from Jpsi electron + any different electron (no real jpsis)
  inclJpsiElePairs->SetLegPDGs(11,-11);
  inclJpsiElePairs->SetMotherPDGs(443,0);
  inclJpsiElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  inclJpsiElePairs->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  die->AddSignalMC(inclJpsiElePairs);

  // prompt J/psi radiative channel
  AliDielectronSignalMC* promptJpsiRad = new AliDielectronSignalMC("promptJpsiRad","Prompt J/psi Radiative");   // prompt J/psi (not from beauty decays)
  promptJpsiRad->SetLegPDGs(11,-11);
  promptJpsiRad->SetMotherPDGs(443,443);
  promptJpsiRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsiRad->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsiRad->SetFillPureMCStep(kTRUE);
  promptJpsiRad->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsiRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsiRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsiRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  promptJpsiRad->SetJpsiRadiative(AliDielectronSignalMC::kIsRadiative);
  //  die->AddSignalMC(promptJpsiRad);

  // prompt J/psi Non radiative channel
  AliDielectronSignalMC* promptJpsiNonRad = new AliDielectronSignalMC("promptJpsiNonRad","Prompt J/psi non-Radiative");   // prompt J/psi (not from beauty decays)
  promptJpsiNonRad->SetLegPDGs(11,-11);
  promptJpsiNonRad->SetMotherPDGs(443,443);
  promptJpsiNonRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsiNonRad->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsiNonRad->SetFillPureMCStep(kTRUE);
  promptJpsiNonRad->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsiNonRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsiNonRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsiNonRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  promptJpsiNonRad->SetJpsiRadiative(AliDielectronSignalMC::kIsNotRadiative);
  //  die->AddSignalMC(promptJpsiNonRad);
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
  (*vec)[1] = 1;
  return vec;
     
}
