void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);
void InitHF(AliDielectron* die, Int_t cutDefinition);

void SetupEventCuts(AliDielectron *die, ULong64_t triggers, Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupV0Cuts( AliDielectron *die,  Int_t cutDefinition);
void SetupPairCuts( AliDielectron *die,  Int_t cutDefinition);

void ConfigEvtPlane(AliDielectron *die,  Int_t cutDefinition);
void ConfigBgrd(    AliDielectron *die,  Int_t cutDefinition);

void AddMCSignals(AliDielectron *die,  Int_t cutDefinition);
void SetEtaCorrection(AliDielectron *die);
TVectorD *GetRunNumbers2011();
TVectorD *GetDeltaPhiBins();

TString names=("Event;Rec;TPC;TOF;TRD;Ionut;TPC-TOF-TRD;NOPID;NOPIDeta;avgPt;LegEff");
enum { kEvent, kRec, kTPC, kTOF, kTRD, kIonut, kTPCTOFTRD, kNoPID, kNoPIDEta, kAvgPt, kLegEff };

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t  isESD = kTRUE;
TString list  = gSystem->Getenv("LIST");

AliDielectron* ConfigJpsi_jb_PbPb(Int_t cutDefinition, Bool_t hasMC=kFALSE, ULong64_t triggers=AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB)
{
  //
  // Setup the instance of AliDielectron
  //

  // gsi train?
  TString trainRoot = gSystem->Getenv("TRAIN_ROOT");
  Bool_t isGSItrain = (trainRoot.IsNull()?kFALSE:kTRUE); 
  if(isGSItrain) {
    // find mc or not?
    if( list.IsNull()) list=prod;
    if( list.Contains("LHC10h")   || list.Contains("LHC11h")   ) hasMC=kFALSE;
    if( list.Contains("LHC11a10") || list.Contains("LHC12a17") ) hasMC=kTRUE;
  }

  //ESD handler?
  isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());

  // switch configurations ON and OFF
  if(hasMC) { // MONTE CARLO
    switch(cutDefinition) {
    case kEvent:   /* */ break;
      //    case kRec: /* */ break;
      //    case kTPC:   /* */ break;
      //    case kTOF:   /* */ break;
      //    case kTRD:   /* */ break;
    case kIonut: /* */ break;
    case kLegEff: /* */ break;
      //    case kNoPID: /* */ break;
    default:         return 0x0;
    }
  } else { // COLLISION DATA
    switch(cutDefinition) {
    case kEvent:   /* */ break;
      //    case kNoPID: /* */ break;
      //case kNoPIDEta: /* */ break;
      //    case kTPC:   /* */ break;
      //    case kTOF:   /* */ break;
      //    case kTRD:   /* */ break;
    case kIonut:       /* */ break;
      //    case kTPCTOFTRD:   /* */ break;
      //    case kAvgPt:       /* */ break;
    default:         return 0x0;
    }
  }

  // task name
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast())  name=arrNames->At(cutDefinition)->GetName();
  printf(" Adding %s%s config %s for %s \n",(isESD?"ESD":"AOD"),(hasMC?" MC":""),name.Data(),list.Data());

  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("%s",name.Data()));
  die->SetHasMC(hasMC);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  SetupEventCuts(die,triggers,cutDefinition);
  SetupTrackCuts(die,cutDefinition);
  SetupV0Cuts(die,cutDefinition);
  SetupPairCuts(die,cutDefinition);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MISC vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // Monte Carlo Signals
  if(hasMC)  AddMCSignals(die, cutDefinition);
  // PID eta correction
  if(cutDefinition!=kNoPID && cutDefinition!=kEvent) SetEtaCorrection(die);
  // bgrd estimators
  if(!hasMC) ConfigBgrd(die,cutDefinition);
  // tpc event plane configuration
  if(!hasMC) ConfigEvtPlane(die,cutDefinition);
  // prefilter settings
  if(cutDefinition==kEvent || cutDefinition==kNoPID || cutDefinition==kNoPIDEta)
    die->SetNoPairing();
  else
    die->SetPreFilterAllSigns();
  //die->SetPreFilterUnlikeOnly();
  // load single electron effieciency map ATTENTION
  if(!hasMC && 0) AliDielectronVarManager::InitEffMap("/home/aliceexp/jbook/gsi-train-trunk/jbook_jpsi/effMap.root");
  if(!hasMC && 1) AliDielectronVarManager::InitEffMap("/sw/alice/gsi-svn/testtrain/effMap.root");

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv OUTPUT vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // histogram setup
  InitHistograms(die,cutDefinition);
  // histogram grid setup
  /*if(!hasMC)*/
  if(cutDefinition!=kEvent && cutDefinition!=kNoPID && cutDefinition!=kNoPIDEta && cutDefinition!=kLegEff)
    InitHF(die,cutDefinition);
  // CF container setup
  if(hasMC && cutDefinition!=kEvent && cutDefinition!=kLegEff)
    InitCF(die,cutDefinition);

  // cut QA
  die->SetCutQA();

  return die;
}

//______________________________________________________________________________________
void SetupEventCuts(AliDielectron *die, ULong64_t triggers, Int_t cutDefinition)
{
  //
  // Setup the event cuts
  //

  // trigger specific centrality cuts (reject trigger inefficiencies)
  Double_t minCent=0.0, maxCent=100.;
  if(!die->GetHasMC()) {
    switch(triggers) {
    // case AliVEvent::kCentral:     minCent= 0.; maxCent=10.; break; //0-9
    // case AliVEvent::kSemiCentral: minCent=10.; maxCent=50.; break; //12-53
    // case AliVEvent::kMB:          minCent=10.; maxCent=90.; break;
    default:                      minCent= 0.; maxCent=90.; break;
    }
  }
  //  if(cutDefinition >= kEtaGap01) {minCent=20.; maxCent=50.;} // v2 analysis

  // VZERO multiplicity vs. number ob global tracks cut
  TF1 *fMean  = new TF1("fMean", "pol1",               0,25e+3);
  fMean->SetParameters(691.633, 1.4892);
  TF1 *fSigma = new TF1("fSigma","[0]+sqrt([1]*x+[2])",0,25e+3);
  fSigma->SetParameters(-83.6599, 36.7677, 69530.7);

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","eventCuts");
  if(!isESD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
  eventCuts->SetCentralityRange(minCent,maxCent);
  //  eventCuts->SetCutOnV0MultipicityNTrks(fMean, fSigma, 4.0);
  eventCuts->Print();
  die->GetEventFilter().AddCuts(eventCuts);

}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //
  Bool_t hasMC=die->GetHasMC();

  // Quality cuts
  //  die->GetTrackFilter().AddCuts(cuts);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  //  trkFilter->SetMinNCrossedRowsOverFindable(0.6);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varAccCuts   = new AliDielectronVarCuts("acc","acc");
  AliDielectronCutGroup  *grpRecCuts = new AliDielectronCutGroup("rec","rec",AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts   *varRecCuts = new AliDielectronVarCuts("VarRecCuts","VarRecCuts");
  AliDielectronTrackCuts *trkRecCuts = new AliDielectronTrackCuts("TrkRecCuts","TrkRecCuts");

  // config specific cuts
  switch(cutDefinition) {
  case kEvent:
  case kRec:
  case kIonut:
  case kTRD:
  case kTPCTOFTRD:
  case kAvgPt:
  case kNoPID:
  case kNoPIDEta:
    varAccCuts->AddCut(AliDielectronVarManager::kPt,           0.85, 20./*1e30*/); ///ATTENTION
    varAccCuts->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8); ////ATTENTION
  case kLegEff:
    varRecCuts->AddCut(AliDielectronVarManager::kTPCclsSegments,7.,   8.0);
    trkRecCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
    //if(hasMC)  varAccCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
    //else
    ///       
    break;
  case kTOF:
  case kTPC:
    varAccCuts->AddCut(AliDielectronVarManager::kPt,           0.95, 1e30);    //0.8
    varAccCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
    varRecCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
    trkRecCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 7); // ITS-3 = 1+2+4
    break;
  }

  // standrad reconstruction cuts
  varRecCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varRecCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varRecCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varRecCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.0);
  //  varRecCuts->AddCut(AliDielectronVarManager::kV0Index0,     0.0);
  trkRecCuts->SetRequireITSRefit(kTRUE);
  trkRecCuts->SetRequireTPCRefit(kTRUE);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronCutGroup *grpPIDCuts = new AliDielectronCutGroup("PID","PID",AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts  *pidSelCuts = new AliDielectronVarCuts("selPIDCuts","selPIDCuts");
  pidSelCuts->AddCut(AliDielectronVarManager::kTRDpidQuality,      4.0,   6.0);
  pidSelCuts->AddCut(AliDielectronVarManager::kTRDchi2,            0.0,   2.0);
  AliDielectronVarCuts *pidVarCuts = new AliDielectronVarCuts("varPIDCuts","varPIDCuts");
  AliDielectronVarCuts *pidMCCuts  = new AliDielectronVarCuts("mcPIDCuts","mcPIDCuts");
  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");

  switch(cutDefinition) {
  case kNoPIDEta:
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-10.,10.);
    break;
  case kTOF:
    pidVarCuts->AddCut(AliDielectronVarManager::kTOFbeta,      0.2,   0.9, kTRUE);
    pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,
		    AliDielectronPID::kIfAvailable);
  case kTPC:
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,4.);
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,4.0,0.,0.,kTRUE);
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-100.,3.5,0.,0.,kTRUE);
    break;
  case kTRD:
    pidCuts->AddCut(AliDielectronPID::kTRDeleEff,AliPID::kElectron,.95,1.,pidSelCuts,
		    kFALSE, AliDielectronPID::kIfAvailable);
  case kEvent:
  case kIonut:
  case kLegEff:
  case kAvgPt:
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.5.,3.);
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-0.9,3.0,-0.5,+0.5,kFALSE,
		    AliDielectronPID::kRequire,AliDielectronVarManager::kEta);
    // pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-0.65,3.0,-0.3,+0.3,kFALSE,
    // 		AliDielectronPID::kRequire,AliDielectronVarManager::kEta);
    // pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-0.4,3.0,-0.1,+0.1,kFALSE,
    // 		AliDielectronPID::kRequire,AliDielectronVarManager::kEta);
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-100.,4.0,0.,0.,kTRUE);
    // tof heavy particle exclusion
    //    pidVarCuts->AddCut(AliDielectronVarManager::kTOFbeta,      0.3,   0.7, kTRUE);
    //pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,
     //		    AliDielectronPID::kIfAvailable);
    break;
  case kTPCTOFTRD:
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.3.,3.);
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-0.9,3.0,-0.5,+0.5,kFALSE,
		    AliDielectronPID::kRequire,AliDielectronVarManager::kEta);
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-0.65,3.0,-0.3,+0.3,kFALSE,
		    AliDielectronPID::kRequire,AliDielectronVarManager::kEta);
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-0.4,3.0,-0.1,+0.1,kFALSE,
		    AliDielectronPID::kRequire,AliDielectronVarManager::kEta);
    pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-100.,4.0,0.,0.,kTRUE);
    pidVarCuts->AddCut(AliDielectronVarManager::kTOFbeta,      0.2,   0.9, kTRUE);
    pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,
		    AliDielectronPID::kIfAvailable);
    pidCuts->AddCut(AliDielectronPID::kTRDeleEff,AliPID::kElectron,.95,1.,pidSelCuts,
		    kFALSE, AliDielectronPID::kIfAvailable);
    break;
  }

  // mc identification
  if(cutDefinition==kLegEff && 0) {
    pidMCCuts->SetCutType(AliDielectronVarCuts::kAny);//only apply any of the two cuts
    pidMCCuts->AddCut(AliDielectronVarManager::kPdgCode,-11.5,-10.5 );
    pidMCCuts->AddCut(AliDielectronVarManager::kPdgCode,10.5,11.5 );
  }


  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TENDER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // exclude conversion electrons selected by the tender
  AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);

  // activate the cut sets (order might be CPU timewise important)
  switch(cutDefinition) {
  case kNoPID:
  case kNoPIDEta:
  case kRec:
    die->GetTrackFilter().AddCuts(varAccCuts);
    grpRecCuts->AddCut(trkRecCuts);
    grpRecCuts->AddCut(varRecCuts);
    die->GetTrackFilter().AddCuts(grpRecCuts);
    break;
  case kEvent:
  case kIonut:
  case kLegEff:
  case kTRD:
  case kTPCTOFTRD:
  case kAvgPt:
  case kTOF:
  case kTPC:
    die->GetTrackFilter().AddCuts(varAccCuts);
    grpRecCuts->AddCut(trkRecCuts);
    grpRecCuts->AddCut(varRecCuts);
    die->GetTrackFilter().AddCuts(grpRecCuts);
    grpPIDCuts->AddCut(pidCuts);
    grpPIDCuts->AddCut(pidVarCuts);
    grpPIDCuts->AddCut(pidMCCuts);
    die->GetTrackFilter().AddCuts(grpPIDCuts);
    //  if(!isESD) cuts->AddCut(trkFilter);
    //cuts->AddCut(noconv);
  }
  //  cuts->Print();

}

//______________________________________________________________________________________
void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the V0 cuts
  //

  // Quality cuts (add the gamma filter to the cut group)
  // TIter next(die->GetTrackFilter().GetCuts());
  // AliAnalysisCuts *cuts;
  // while((cuts = (AliAnalysisCuts*)next())) {
  //   if(cuts->IsA() == AliDielectronCutGroup::Class())  break;
  // }

  AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("V0","V0");
  gammaV0Cuts->SetPdgCodes(22,11,11);
  gammaV0Cuts->SetDefaultPID(16);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
  //  gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle,              0.0,   0.1, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  //  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // not sure if it works as expected
  gammaV0Cuts->SetExcludeTracks(kTRUE);
  gammaV0Cuts->Print();

 //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
 //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; &&  const Double_t cutQTG2 < 0.04;

  // if(cuts)
  //   ((AliDielectronCutGroup*)cuts)->AddCut(gammaV0Cuts);
  // else
  die->GetTrackFilter().AddCuts(gammaV0Cuts);
}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //
  Bool_t hasMC=die->GetHasMC();

  // rap and mass rejection
  Double_t gCut=0.05, yCut=0.9;
  switch(cutDefinition) {
  case kEvent:      yCut=0.0;  break;
  case kRec:        yCut=0.9;  break;
  case kNoPID:      return;
  case kNoPIDEta:   return;
  case kTPC:        yCut=0.9;  break;
  case kTOF:        yCut=0.9;  break;
  case kTRD:        yCut=0.8;  break;
  case kIonut:      yCut=0.8;  break;
  case kLegEff:     return;
  case kTPCTOFTRD:  yCut=0.8;  break;
  case kAvgPt:      yCut=0.8;  break;
    //  default: gCut=0.05;       // default
  }

  // MC
  //if(hasMC) yCut=0.9;

  // rapidity selection
  AliDielectronVarCuts *rapCut=new AliDielectronVarCuts(Form("|Y|<%.1f",yCut),Form("|Y|<%.1f",yCut));
  rapCut->AddCut(AliDielectronVarManager::kY,-1.*yCut,yCut);
  die->GetPairFilter().AddCuts(rapCut);

  // gamma rejection
  AliDielectronVarCuts *gammaCuts = new AliDielectronVarCuts("GammaCuts","GammaCuts");
  gammaCuts->AddCut(AliDielectronVarManager::kM,            0.0,   gCut);
  die->GetPairPreFilter().AddCuts(gammaCuts);

}

//______________________________________________________________________________________
void ConfigBgrd(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Configurate the background estimators
  //

  // skip event QA
  switch(cutDefinition) {
  case kEvent:
  case kNoPID:
  case kNoPIDEta:
  case kRec:
    return;
  }

  // add track rotations
  AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
  rot->SetIterations(10);
  rot->SetConeAnglePhi(TMath::Pi());
  rot->SetStartAnglePhi(TMath::Pi());
  //      die->SetTrackRotator(rot);

  // add mixed events
  AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
  switch(cutDefinition) {
  // case kIonut:
  //   mix->AddVariable(AliDielectronVarManager::kZvPrim,      20, -10., 10.);
  //   mix->AddVariable(AliDielectronVarManager::kCentrality,  36,   0., 90.);
  //   mix->AddVariable(AliDielectronVarManager::kv0ACrpH2,    10,   0., TMath::Pi());
  //   break;
  default: 
    mix->AddVariable(AliDielectronVarManager::kZvPrim,      "-10.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,10.");
    mix->AddVariable(AliDielectronVarManager::kCentrality,  9,  0.,90.);
    mix->AddVariable(AliDielectronVarManager::kTPCrpH2,     10, TMath::Pi()/-2, TMath::Pi()/2); // max res: 10%->10bins // 8bins
    //    mix->AddVariable(AliDielectronVarManager::kTPCmagH2,    "0.,20.,50.,80.,110.,150.,500.");
    break;
  }
  mix->SetSkipFirstEvent(kTRUE); // needed for flow analysis
  mix->SetMixType(AliDielectronMixingHandler::kAll);
  mix->SetDepth(150);
  die->SetMixingHandler(mix);

}

//______________________________________________________________________________________
void ConfigEvtPlane(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Configurate the TPC event plane 
  //

  switch(cutDefinition) {
  case kNoPID:
  case kNoPIDEta:
  case kRec:
    return;
  }

  //   Double_t gGap;
  //   switch(cutDefinition) {
  //   case kEtaGap01:   gGap=0.1;   break;
  //   case kEtaGap02:   gGap=0.2;   break;
  //   case kEtaGap03:   gGap=0.3;   break;
  //   case kEtaGap04:   gGap=0.4;   break;
  //   case kEtaGap05:   gGap=0.5;   break;
  //   default: gGap=0.0;
  //   }
  // eta gap in tpc event plane
  //   AliDielectronVarCuts *etaGap = new AliDielectronVarCuts(AliDielectronVarManager::GetValueName(AliDielectronVarManager::kEta),"etaGap");
  //   etaGap->AddCut(AliDielectronVarManager::kEta,-1*gGap,gGap,kTRUE);
  //   die->GetEventPlanePreFilter().AddCuts(etaGap);

  AliDielectronVarCuts *poi = new AliDielectronVarCuts("PoI","PoI");
  poi->AddCut(AliDielectronVarManager::kM,2.92,3.20);     // particles of interest, jpsi mass window
  die->GetEventPlanePOIPreFilter().AddCuts(poi);

  // die->SetLikeSignSubEvents();
  die->SetPreFilterEventPlane();
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  Bool_t hasMC=die->GetHasMC();

  // booleans for histo selection
  Bool_t bHistTrackQA=kFALSE, bHistEvts = kFALSE, bHistPair = kFALSE, bHistPairME = kFALSE, bHistFlow = kFALSE, bHistFlowQA=kFALSE, bHistPID=kFALSE;
  switch (cutDefinition) {
  case kEvent: bHistEvts=kTRUE; break;
  case kRec:
  case kNoPID:
  case kNoPIDEta:
  case kTPC:
  case kTOF:   //bHistEvts=kTRUE; //bHistFlow=kTRUE;
  case kTRD:   ///bHistEvts=kTRUE; //bHistFlow=kTRUE;
  case kIonut: ///*bHistEvts=kTRUE;*/ bHistPair=kTRUE; break; //bHistPairME=kTRUE;
  case kTPCTOFTRD: bHistPair=kTRUE; bHistPID=kTRUE; break; //bHistPairME=kTRUE;
  case kAvgPt: /* */ break;
  }
  if(hasMC) {
    bHistPID=kFALSE;
  }

  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());

  //add histograms to event class
  histos->AddClass("Event");
  Int_t maxMixBins = (die->GetMixingHandler() ? die->GetMixingHandler()->GetNumberOfBins() : 0);
  histos->UserHistogram("Event","","", 100, 0.0, 100.0,   AliDielectronVarManager::kCentrality);
  if(!hasMC && die->GetMixingHandler() )
    histos->UserHistogram("Event","","", maxMixBins, 0, maxMixBins, AliDielectronVarManager::kMixingBin);

  ////// EVENT HISTOS /////
  if(bHistEvts) {
    histos->UserHistogram("Event","","", 100,0.,100.,	  AliDielectronVarManager::kNTrk);
    //    histos->UserHistogram("Event","","", 100,0.,100.,	  AliDielectronVarManager::kNacc);
    histos->UserHistogram("Event","","", GetRunNumbers2011(), AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Event","","", 300,-15.,15.,	  AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","","", GetRunNumbers2011(), AliDielectronHelper::MakeLinBinning(3000, 0., 3000.),
			  AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kNacc);
    histos->UserHistogram("Event","","", AliDielectronHelper::MakeLinBinning(90, 0., 90.), AliDielectronHelper::MakeLinBinning(3000, 0., 3000.),
			  AliDielectronVarManager::kCentrality, AliDielectronVarManager::kNacc);
    histos->UserProfile(  "Event","","", AliDielectronVarManager::kNacc,        90, 0., 90.,  AliDielectronVarManager::kCentrality);
    histos->UserProfile(  "Event","","", AliDielectronVarManager::kNVtxContrib, 90, 0., 90.,  AliDielectronVarManager::kCentrality);

    if(!hasMC) {
      // event plane resolutions
      TObjArray *limits  = new TObjArray();
      limits->Add(AliDielectronHelper::MakeLinBinning( 18,  0.,    90.));
      //    limits->Add(GetRunNumbers2011());
      limits->Add(AliDielectronHelper::MakeLinBinning(100,  -1.,    1.));
      limits->Add(AliDielectronHelper::MakeLinBinning(100,  -1.,    1.));
      limits->Add(AliDielectronHelper::MakeLinBinning(100,  -1.,    1.));

      UInt_t var[4]={AliDielectronVarManager::kCentrality,
		     //		   AliDielectronVarManager::kRunNumber,
		     AliDielectronVarManager::kv0ATPCDiffH2,
		     AliDielectronVarManager::kv0CTPCDiffH2,
		     AliDielectronVarManager::kv0Av0CDiffH2 };
      histos->UserSparse("Event", 4, limits, var);

      // // event plane flattening
      TObjArray *limits2  = new TObjArray();
      limits2->Add(AliDielectronHelper::MakeLinBinning( 18,   0.,    90.));
      //limits2->Add(GetRunNumbers2011());
      limits2->Add(AliDielectronHelper::MakeLinBinning(100,  -1.6,   1.6));
      limits2->Add(AliDielectronHelper::MakeLinBinning(100,  -1.6,   1.6));
      limits2->Add(AliDielectronHelper::MakeLinBinning(100,  -1.6.,  1.6));

      UInt_t var2[4]={AliDielectronVarManager::kCentrality,
		      // 		    AliDielectronVarManager::kRunNumber,
		      AliDielectronVarManager::kv0ArpH2,
		      AliDielectronVarManager::kv0CrpH2,
		      AliDielectronVarManager::kTPCrpH2 };
      //    histos->UserSparse("Event", 4, limits2, var2);
    } //end: hasMC

  } //hist: event

  ////// FLOW //////
  if(bHistFlow) {
    histos->UserHistogram("Event","","", 100,-1500.,1500., AliDielectronVarManager::kTPCxH2);
    histos->UserHistogram("Event","","", 100,-1500.,1500., AliDielectronVarManager::kTPCyH2);
    histos->UserHistogram("Event","","", 100,-2.,2., AliDielectronVarManager::kTPCrpH2);
    histos->UserHistogram("Event","","", 100,-2.,2., AliDielectronVarManager::kv0ArpH2);
    histos->UserHistogram("Event","","", 100,-2.,2., AliDielectronVarManager::kv0CrpH2);
    histos->UserHistogram("Event","","", 100,-2.,2., AliDielectronVarManager::kv0ACrpH2);
    histos->UserHistogram("Event","","", 100,-2.,2., AliDielectronVarManager::kv0A0rpH2);
    histos->UserHistogram("Event","","", 100,-2.,2., AliDielectronVarManager::kv0C0rpH2);
    histos->UserHistogram("Event","","", 100,-2.,2., AliDielectronVarManager::kv0A3rpH2);
    histos->UserHistogram("Event","","", 100,-2.,2., AliDielectronVarManager::kv0C3rpH2);
    // EP resolution calculation centrality dependence
    // histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0ATPCDiffH2);
    // histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0CTPCDiffH2);
    // histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Av0CDiffH2);
  } // hist: flow

  // DEBUG
  if(bHistFlowQA) {
    // TPC event plane
    histos->UserHistogram("Event","","", 100,-1500.,1500., AliDielectronVarManager::kTPCxH2uc);
    histos->UserHistogram("Event","","", 100,-1500.,1500., AliDielectronVarManager::kTPCyH2uc);
    histos->UserHistogram("Event","","", 100,   -2.,   2., AliDielectronVarManager::kTPCrpH2uc);
    histos->UserHistogram("Event","","", 100,   -2.,   2., AliDielectronVarManager::kTPCsub1rpH2);
    histos->UserHistogram("Event","","", 100,   -2.,   2., AliDielectronVarManager::kTPCsub2rpH2);
    histos->UserHistogram("Event","","", 100,   -1.,   1., AliDielectronVarManager::kTPCsub12DiffH2);
    // detector effects
    histos->UserHistogram("Event","","", 10,0.,100., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kTPCsub12DiffH2Sin);
    // recentering stuff
    histos->UserProfile("Event","","", AliDielectronVarManager::kTPCxH2,
			AliDielectronHelper::MakeLinBinning(8, 0.,80.), GetRunNumbers2011(),
			AliDielectronVarManager::kCentrality, AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","","", AliDielectronVarManager::kTPCyH2,
     			AliDielectronHelper::MakeLinBinning(8, 0.,80.), GetRunNumbers2011(),
     			AliDielectronVarManager::kCentrality, AliDielectronVarManager::kRunNumber);
    // EP resolution calculation
    histos->UserHistogram("Event","","", 200,-1.0,1.0,  AliDielectronVarManager::kv0ATPCDiffH2);
    histos->UserHistogram("Event","","", 200,-1.0,1.0,  AliDielectronVarManager::kv0CTPCDiffH2);
    histos->UserHistogram("Event","","", 200,-1.0,1.0,  AliDielectronVarManager::kv0Av0CDiffH2);
    histos->UserHistogram("Event","","", 200,-1.0,1.0,  AliDielectronVarManager::kTPCsub12DiffH2);
    histos->UserHistogram("Event","","", 200,-1.0,1.0,  AliDielectronVarManager::kv0Av0C0DiffH2);
    histos->UserHistogram("Event","","", 200,-1.0,1.0,  AliDielectronVarManager::kv0Av0C3DiffH2);
    histos->UserHistogram("Event","","", 200,-1.0,1.0,  AliDielectronVarManager::kv0Cv0A0DiffH2);
    histos->UserHistogram("Event","","", 200,-1.0,1.0,  AliDielectronVarManager::kv0Cv0A3DiffH2);
    histos->UserHistogram("Event","","", 200,-1.0,1.0,  AliDielectronVarManager::kv0A0v0A3DiffH2);
    histos->UserHistogram("Event","","", 200,-1.0,1.0,  AliDielectronVarManager::kv0C0v0C3DiffH2);
    // EP resolution calculation centrality dependence
    histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0ATPCDiffH2);
    histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0CTPCDiffH2);
    histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Av0CDiffH2);
    // histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kTPCsub12DiffH2);
    // histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Av0C0DiffH2);
    // histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Av0C3DiffH2);
    // histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Cv0A0DiffH2);
    // histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Cv0A3DiffH2);
    // histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0A0v0A3DiffH2);
    // histos->UserHistogram("Event","","", 18,0.,90., 100,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0C0v0C3DiffH2);
    // EP angle correlation
    histos->UserHistogram("Event","","", 320,-3.2.,3.2, 320,-3.2.,3.2, 
			  AliDielectronVarManager::kTPCrpH2uc,AliDielectronVarManager::kv0ACrpH2);
    // EP Qvector magnitudes // TODO move to QA
    histos->UserHistogram("Event","","", 200,0.,200., AliDielectronVarManager::kTPCmagH2uc);
    histos->UserHistogram("Event","","", 200,0.,800., AliDielectronVarManager::kv0ACmagH2);
    histos->UserHistogram("Event","","", 200,0.,800., AliDielectronVarManager::kv0AmagH2);
    histos->UserHistogram("Event","","", 200,0.,800., AliDielectronVarManager::kv0CmagH2);
  } //hist: flowQA


  ///// PAIR HISTOS /////
  if(bHistPair) {
    //Initialise histogram classes
    histos->SetReservedWords("Track;Pair");

    //Pair classes
    // to fill also mixed event histograms loop until 7 or 10
    for (Int_t i=0; i<(bHistPairME ? 8 : 3); ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
    }
    //add MC signal histograms to pair class
    if(die->GetMCSignals()) {
      for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i)
	histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
    }

    //Track classes
    //legs from pair (fill SE)
    for (Int_t i=1; i<2; ++i){
      histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
    }
    //to fill also track info from 2nd event loop until 2
    //    for (Int_t i=0; i<2; ++i) histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
    histos->AddClass(Form("Track_%s",     AliDielectron::PairClassName(AliDielectron::kEv1PM)));


    ////// TRACK HISTOS /////
    // Vertex
    // histos->UserHistogram("Track","","", 500,-1.,1., AliDielectronVarManager::kImpactParXY);
    // histos->UserHistogram("Track","","", 600,-3.,3., AliDielectronVarManager::kImpactParZ);
    // Kinematics
    histos->UserHistogram("Track","","", 400,0,20.,  AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","","", 200,-1,1, 200,0,6.285, AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    // TPC
    // histos->UserHistogram("Track","","", 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC);
    // histos->UserHistogram("Track","","", 160,-0.5,159.5, AliDielectronVarManager::kTPCsignalN);
    // histos->UserHistogram("Track","","", 160,-0.5,159.5, AliDielectronVarManager::kNFclsTPCr);
    // histos->UserHistogram("Track","","", 160,-0.5,159.5, 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
    // TRD
    histos->UserHistogram("Track","","",   8,-0.5,  7.5, AliDielectronVarManager::kTRDpidQuality);
    // PID
    // histos->UserHistogram("Track","","", 400,0.2,20.,200,0.,200., AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
    // histos->UserHistogram("Track","","", 400,0.2,20.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
    histos->UserHistogram("Track","","", 250,0.0,5.0,300,0.,1.2,  AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta,kTRUE);
    //    histos->UserHistogram("Track","","", 100,-1.,+1.,200,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);

    histos->UserHistogram("Track","","", 100,0.2,10.,100,0.,200.,
			  AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
    histos->UserHistogram("Track","","", 100,0.2,10.,100,-5.,+5.,
			  AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
    histos->UserHistogram("Track","","", 100,0.,4000.,100,0.,200.,
			  AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","","", 100,0.,4000.,100,-5.,+5.,
			  AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","","", 100,0.,100.,100,0.,200.,
			  AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","","", 100,0.,100.,100,-5.,+5.,
			  AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","","", 100,-1.,+1.,100,0.,200.,
			  AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","","", 100,-1.,+1.,100,-5.,+5.,
			  AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","","", 100,-1.,+1.,100,-5.,+5.,
			  AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw);
    
    if(!hasMC) histos->UserHistogram("Track","","", 150,0.0,1.5, AliDielectronVarManager::kLegEff);
    if(!hasMC) histos->UserProfile("Track","","", AliDielectronVarManager::kLegEff, 18,0.0,90.0, AliDielectronVarManager::kCentrality);
    if(!hasMC) histos->UserProfile("Track","","", AliDielectronVarManager::kLegEff, 100,0.0,20.0, AliDielectronVarManager::kPt);
    if(!hasMC) histos->UserProfile("Track","","", AliDielectronVarManager::kLegEff, 20,0.0,TMath::TwoPi(), AliDielectronVarManager::kPhi);
    if(!hasMC) histos->UserProfile("Track","","", AliDielectronVarManager::kLegEff, 20,-1.,+1., AliDielectronVarManager::kEta);

    // TPC PID (detailed) //
    if(bHistPID) {
      // arbitrary binning for variables
      TVectorD *vcen = AliDielectronHelper::MakeLinBinning( 21,  0.,    52.5);
      (*vcen)[21] = 90.;

      TObjArray *limits  = new TObjArray();
      limits->Add(AliDielectronHelper::MakeLinBinning( 75,-10.,     5.));
      //  limits->Add(AliDielectronHelper::MakeLinBinning(100,  0.,   200.));
      limits->Add(AliDielectronHelper::MakeLinBinning(60,   0.,    1.2));
      limits->Add(AliDielectronHelper::MakeLinBinning(100,  0.,  4000.));
      limits->Add(AliDielectronHelper::MakeLinBinning( 50,  0.,    10.));
      limits->Add(vcen);
      //    limits->Add(AliDielectronHelper::MakeLinBinning( 36,  0.,    90.));
      limits->Add(AliDielectronHelper::MakeLinBinning( 20, -1.,     1.));
      //  limits->Add(GetRunNumbers());
      UInt_t var[6]={AliDielectronVarManager::kTPCnSigmaEle,
		     AliDielectronVarManager::kTOFbeta,
		     //              AliDielectronVarManager::kTPCsignal,
		     AliDielectronVarManager::kRefMultTPConly,
		     AliDielectronVarManager::kPIn,
		     AliDielectronVarManager::kCentrality,
		     AliDielectronVarManager::kEta
		     //      AliDielectronVarManager::kRunNumber
      };
      histos->UserSparse(Form("Track_%s",     AliDielectron::PairClassName(AliDielectron::kEv1PM)), 6, limits, var);      
    }

    ///// Pair classes /////
    histos->UserHistogram("Pair","","",  125,.0,125*0.04, AliDielectronVarManager::kM); // 40MeV bins, 12GeV/c2
    histos->UserHistogram("Pair","","",  100,-1.,1.,      AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","","",  400,0,20.,       AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","","",  100,0.,3.15,     AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","","",  100,0.,20,       AliDielectronVarManager::kChi2NDF);
    histos->UserHistogram("Pair","","",  100,0.,3.15,     AliDielectronVarManager::kPsiPair);
    histos->UserHistogram("Pair","","",  200,0.,100.,     AliDielectronVarManager::kR);
    histos->UserHistogram("Pair","","",   50,0.,5.,       AliDielectronVarManager::kLegDist);
    histos->UserHistogram("Pair","","",   50,0.,5.,       AliDielectronVarManager::kLegDistXY);
    histos->UserHistogram("Pair","","", 100,-1*TMath::Pi(),+1*TMath::Pi(), AliDielectronVarManager::kDeltaPhiv0ArpH2);
    histos->UserHistogram("Pair","","", 100,-1*TMath::Pi(),+1*TMath::Pi(), AliDielectronVarManager::kDeltaPhiv0CrpH2);
    histos->UserHistogram("Pair","","", 100,-1*TMath::Pi(),+1*TMath::Pi(), AliDielectronVarManager::kDeltaPhiTPCrpH2);
    //    histos->UserHistogram("Pair","","", 210,-1.05,1.05, 100,0.,2.5, AliDielectronVarManager::kArmAlpha,AliDielectronVarManager::kArmPt);
    if(!hasMC) histos->UserProfile("Pair","","", AliDielectronVarManager::kPairEff, 125,.0,125*0.04, AliDielectronVarManager::kM);
    if(!hasMC) histos->UserProfile("Pair","","", AliDielectronVarManager::kPairEff, 100,.0,15.0,     AliDielectronVarManager::kPt);

    if(!hasMC && die->GetMixingHandler() )
      histos->UserHistogram("Pair","","", maxMixBins, 0, maxMixBins, AliDielectronVarManager::kMixingBin);


  } //hist: pair


  ////// MONTE CARLO //////
  if(hasMC && cutDefinition==kLegEff && 1) { ////////ATTENTION

    //Initialise histogram classes
    histos->SetReservedWords("Track");

    //add MC signal histograms to track class
    if(die->GetMCSignals()) {
      TString className = Form("Track_%s_%s",AliDielectron::PairClassName(AliDielectron::kEv1PM),die->GetMCSignals()->Last()->GetName());
      histos->AddClass(Form("%s_MCtruth",className.Data()));
      histos->AddClass(className.Data());

      //workaround for AN-09 dummy
      for(Int_t isig=0; isig<die->GetMCSignals()->GetEntriesFast(); isig++) {
	histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(isig)->GetName()));
	histos->AddClass(Form("Pair_%s_MCtruth",die->GetMCSignals()->At(isig)->GetName()));
	// histos->AddClass(Form("Track_Legs_%s",die->GetMCSignals()->At(isig)->GetName()));
	// histos->AddClass(Form("Track_Legs_%s_MCtruth",die->GetMCSignals()->At(isig)->GetName()));
	histos->AddClass(Form("Track_%s_%s",AliDielectron::PairClassName(AliDielectron::kEv1PM),die->GetMCSignals()->At(isig)->GetName()));
	histos->AddClass(Form("Track_%s_%s_MCtruth",AliDielectron::PairClassName(AliDielectron::kEv1PM),die->GetMCSignals()->At(isig)->GetName()));
      }

      // single leg efficiency //
      // arbitrary binning for variables
      // TVectorD *vpt = AliDielectronHelper::MakeLinBinning( 41,  0.0,    10.25);
      // (*vpt)[41] = 20.;

      TObjArray *limEpm  = new TObjArray();
      //    limEpm->Add(AliDielectronHelper::MakeLinBinning( 75,-10.,     5.));
      limEpm->Add(AliDielectronHelper::MakeLinBinning( 18,  0.,    90.));
      limEpm->Add(AliDielectronHelper::MakeLinBinning( 80,  0.,    20.));
      //      limEpm->Add(vpt);
      limEpm->Add(AliDielectronHelper::MakeLinBinning( 20,  0.,   TMath::TwoPi()));
      limEpm->Add(AliDielectronHelper::MakeLinBinning( 20, -1.,    +1.));
      //      limEpm->Add(AliDielectronHelper::MakeLinBinning( 20, -2.,    +2.));
      //  limEpm->Add(GetRunNumbers());
      UInt_t varEpm[]={//AliDielectronVarManager::kTPCnSigmaEle,
	AliDielectronVarManager::kCentrality
	,AliDielectronVarManager::kPt
	,AliDielectronVarManager::kPhi
	,AliDielectronVarManager::kEta
	//	,AliDielectronVarManager::kY
	//      ,AliDielectronVarManager::kRunNumber
      };

      histos->UserSparse("Track", limEpm->GetEntriesFast(), limEpm, varEpm);
    }
  } //end: hasMC



  /*
    if(cutDefinition == kTOFTRD && hasMC) {
    histos->AddClass("MCEvent");
    histos->UserHistogram("MCEvent","Cent_NJPsis","Centrality vs. generated incl. J/#psi per event;centrality (%);N_{J/#psi}",
    10,0.,100., 21,-0.5,20.5,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kNumberOfJPsis);
    }
  */

  die->SetHistogramManager(histos);


  ////// LOG //////
  TIter nextClass(histos->GetHistogramList());
  THashList *l=0;
  while ( (l=static_cast<THashList*>(nextClass())) ) {
    //printf(" [D] HistogramManger: Class %s: Histograms: %04d \n", l->GetName(), l->GetEntries());
  }


}

void InitHF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the HF arrays
  //
  Bool_t hasMC=die->GetHasMC();

  // booleans for histo selection
  Bool_t bHistFlow = kTRUE, bHistPt = kTRUE, bHistEff = kTRUE; //ATTENTION
  switch (cutDefinition) {
  case kAvgPt:  bHistFlow=kFALSE;; break;
    //  default:      bHistPt=kTRUE; bHistFlow = kTRUE; bHistEff=kTRUE; break;
  }

  AliDielectronHF *hf=new AliDielectronHF(die->GetName(),die->GetTitle());
  // define pair types and sources
  if(hasMC) hf->SetStepForMCGenerated();
  //  else
  if(hasMC)  hf->SetPairTypes(AliDielectronHF::kOSandMIX);
    //  hf->SetPairTypes(AliDielectronHF::kOSonly);
  if(!hasMC) hf->SetPairTypes(AliDielectronHF::kAll);

  //// define output objects ////

  // defaults
  hf->UserHistogram("Pair", AliDielectronHelper::MakeLinBinning(125,.0,125*0.04), AliDielectronVarManager::kM);
  // for mixed event weighting
  Int_t maxMixBins = (die->GetMixingHandler() ? die->GetMixingHandler()->GetNumberOfBins() : 0);
  //  if(maxMixBins)  hf->UserHistogram("Pair", AliDielectronHelper::MakeLinBinning(maxMixBins,0,maxMixBins), AliDielectronVarManager::kMixingBin);

  // mean pt analysis
  if(bHistPt) {
    hf->UserProfile("Pair", AliDielectronVarManager::kPt,
		    AliDielectronHelper::MakeLinBinning(125,.0,125*0.04), AliDielectronVarManager::kM);
  }

  // on the fly efficienies
  if(bHistEff && !hasMC) {
    hf->UserProfile("Pair", AliDielectronVarManager::kPairEff, 
		    AliDielectronHelper::MakeLinBinning(125,.0,125*0.04), AliDielectronVarManager::kM);
    hf->UserProfile("Pair", AliDielectronVarManager::kPt,
		    AliDielectronHelper::MakeLinBinning(125,.0,125*0.04), AliDielectronVarManager::kM, 
		    "", AliDielectronVarManager::kOneOverPairEff);
    hf->UserHistogram("Pair", AliDielectronHelper::MakeLinBinning(125,.0,125*0.04), AliDielectronVarManager::kM,
		      AliDielectronVarManager::kOneOverPairEff);
  }

  // flow analysis
  if(bHistFlow && !hasMC) {  
    // flow versus minv
    hf->UserProfile("Pair", AliDielectronVarManager::kv0ArpH2FlowV2,
		    AliDielectronHelper::MakeLinBinning(125,.0,125*0.04), AliDielectronVarManager::kM);
    hf->UserProfile("Pair", AliDielectronVarManager::kv0CrpH2FlowV2,
		    AliDielectronHelper::MakeLinBinning(125,.0,125*0.04), AliDielectronVarManager::kM);
    hf->UserProfile("Pair", AliDielectronVarManager::kTPCrpH2FlowV2,
		    AliDielectronHelper::MakeLinBinning(125,.0,125*0.04), AliDielectronVarManager::kM);
    // detector effects
    hf->UserProfile("Pair", AliDielectronVarManager::kCosTPCrpH2,
		    AliDielectronHelper::MakeLinBinning(1,.0,125*0.04), AliDielectronVarManager::kM);
    hf->UserProfile("Pair", AliDielectronVarManager::kSinTPCrpH2,
		    AliDielectronHelper::MakeLinBinning(1,.0,125*0.04), AliDielectronVarManager::kM);
    hf->UserProfile("Pair", AliDielectronVarManager::kCosPhiH2,
		    AliDielectronHelper::MakeLinBinning(1,.0,125*0.04), AliDielectronVarManager::kM);
    hf->UserProfile("Pair", AliDielectronVarManager::kSinPhiH2,
		    AliDielectronHelper::MakeLinBinning(1,.0,125*0.04), AliDielectronVarManager::kM);
    hf->UserProfile("Pair", AliDielectronVarManager::kTPCrpH2FlowV2Sin,
		    AliDielectronHelper::MakeLinBinning(125,.0,125*0.04), AliDielectronVarManager::kM);
  }

  //// define the grid size and granularity /////
  // event variables //
  //  hf->AddCutVariable(AliDielectronVarManager::kCentrality, AliDielectronHelper::MakeArbitraryBinning("0,10,50,90")); // flow only
  if(hasMC) hf->AddCutVariable(AliDielectronVarManager::kCentrality,               18, 0., 90.);
  else      hf->AddCutVariable(AliDielectronVarManager::kCentrality,               9, 0., 90.);
  //  if(hasMC)  hf->AddCutVariable(AliDielectronVarManager::kRunNumber, GetRunNumbers2011() );
  //  if(hasMC)  hf->AddCutVariable(AliDielectronVarManager::kNacc,        3000,0.,3000.);
  //  if(hasMC)  hf->AddCutVariable(AliDielectronVarManager::kNVtxContrib,   20,0.,4000.);

  // pair variables //ATTENTION
  if(hasMC && 0) hf->AddCutVariable(AliDielectronVarManager::kY,  18, -0.9,  0.9);
  //  hf->AddCutVariable(AliDielectronVarManager::kPt, 10,  0.0, 10.0);
  if(!hasMC)
     hf->AddCutVariable(AliDielectronVarManager::kPt, AliDielectronHelper::MakeArbitraryBinning("0,1,2,3,4,5,6,7,8,9,10,100"));
  //  hf->AddCutVariable(AliDielectronVarManager::kPt, 20,  0.0, 10.0);
  else
    hf->AddCutVariable(AliDielectronVarManager::kPt, 50,  0.0, 10.0);

  //    if(hasMC) hf->AddCutVariable(AliDielectronVarManager::kTRDpidEffPair,101,0.0,1.01);
  //    if(hasMC) hf->AddCutVariable(AliDielectronVarManager::kThetaCS,15,-1.,1.);

  // flow variables //
  //  if(!hasMC) hf->AddCutVariable(AliDielectronVarManager::kDeltaPhiv0ArpH2, GetDeltaPhiBins());
  //  if(!hasMC) hf->AddCutVariable(AliDielectronVarManager::kDeltaPhiv0CrpH2, GetDeltaPhiBins());
  if(!hasMC && bHistFlow) hf->AddCutVariable(AliDielectronVarManager::kDeltaPhiTPCrpH2, GetDeltaPhiBins());

  // leg variables // NOTE: switched off in HF??
  //  if(hasMC) hf->AddCutVariable(AliDielectronVarManager::kPt, "0.85, 0.95, 1.1, 100.0", kTRUE, AliDielectronHF::kBinToMax);
  // if(hasMC) hf->AddCutVariable(AliDielectronVarManager::kEta,"-0.9,-0.8,0.8,0.9",      kTRUE, AliDielectronHF::kSymBin);
  //  hf->AddCutVariable(AliDielectronVarManager::kY,           1, -0.9, 0.9                     );
  //  hf->AddCutVariable(AliDielectronVarManager::kPt,          "0.8, 1.0, 1.1, 1.2, 1.5, 100.0", kTRUE, AliDielectronHF::kBinToMax);
  //  hf->AddCutVariable(AliDielectronVarManager::kNclsTPC,     "70,90,100,120,160",              kTRUE, AliDielectronHF::kBinToMax);
  //  hf->AddCutVariable(AliDielectronVarManager::kTPCnSigmaEle,"-4,-3,-2.5,-2,2,2.5,3,4",        kTRUE, AliDielectronHF::kSymBin);
  //hf->AddCutVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.,3.5,4.,100.",                 kTRUE, AliDielectronHF::kBinToMax);
  //hf->AddCutVariable(AliDielectronVarManager::kITSLayerFirstCls,4,0.,4.,              kFALSE, kTRUE, AliDielectronHF::kBinFromMin);
  //hf->AddCutVariable(AliDielectronVarManager::kNclsITS,         5,2.,7.,              kFALSE, kTRUE, AliDielectronHF::kBinToMax);
  //hf->AddCutVariable(AliDielectronVarManager::kRunNumber,  GetRunNumbers2011());

  die->SetHistogramArray(hf);
}

void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the CF Manager if needed
  //
  Bool_t hasMC=die->GetHasMC();

  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  //event variables
  cf->AddVariable(AliDielectronVarManager::kCentrality,"0.,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,70.0,80.,90.");
  cf->AddVariable(AliDielectronVarManager::kRunNumber, GetRunNumbers2011() );

  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,100,0.0,100*0.2);
  
  //  cf->AddVariable(AliDielectronVarManager::kY,"-5,-1,-0.9,-0.8,-0.5,0.5,0.8,0.9,1.0,5");
  cf->AddVariable(AliDielectronVarManager::kY,36,-0.9,+0.9);
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,1,1,1);
  //  cf->AddVariable(AliDielectronVarManager::kThetaCS,20,-1.,+1.);
  //  cf->AddVariable(AliDielectronVarManager::kThetaHE,20,-1.,+1.);
  //  cf->AddVariable(AliDielectronVarManager::kPhiCS,20,-3.2,+3.2);
  //  cf->AddVariable(AliDielectronVarManager::kPhiHE,20,-3.2,+3.2);
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.75, 0.80, 0.85, 1.0, 1.1, 100.0",kTRUE);
  //  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 70, 80, 90, 100, 120, 160",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-0.9,-0.85,-0.8,-0.75,-0.70,0.70,0.75,0.8,0.85,0.9",kTRUE);
  /*
  // event variables
  cf->AddVariable(AliDielectronVarManager::kCentrality,               18,0.,  90.);
  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kRunNumber, GetRunNumbers2011() );
  //  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kNacc,        3000,0.,3000.);
  //  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kNVtxContrib,   20,0.,4000.);

  // pair variables
  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kPairType,1,1,1);
  else       cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);
  cf->AddVariable(AliDielectronVarManager::kM, 125,  0.0,  5.0); //40Mev Steps
  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kY,  18, -0.9,  0.9);
  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kPt, 100,  0.0, 20.0);
  else       cf->AddVariable(AliDielectronVarManager::kPt, AliDielectronHelper::MakeArbitraryBinning("0,1,2,3,4,5,6,7,8,9,10,20"));
  //    if(hasMC) cf->AddVariable(AliDielectronVarManager::kTRDpidEffPair,101,0.0,1.01);
  //    if(hasMC) cf->AddVariable(AliDielectronVarManager::kThetaCS,15,-1.,1.);

  // flow variables
  //if(!hasMC) cf->AddVariable(AliDielectronVarManager::kDeltaPhiv0ArpH2, GetDeltaPhiBins());
  //if(!hasMC) cf->AddVariable(AliDielectronVarManager::kDeltaPhiv0CrpH2, GetDeltaPhiBins());

  // leg variables
  if(hasMC) cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.85, 0.95, 1.0, 1.1, 100.0",kTRUE);
  if(hasMC) cf->AddVariable(AliDielectronVarManager::kEta,"-0.9,-0.8,-0.7, 0.7, 0.8, 0.9",  kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,7,-1.5,5.5,kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kNclsITS,"1,2,3,4,5,6",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-3,-2.5,-2,2,2.5,3",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"2.5,3.0,3.5,4.0,4.5,100",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kNclsTPC,"70, 90, 100, 120, 160",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.5,4.0,4.5,5.0,100",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kTOFnSigmaEle,"-3,-2,2,3",kTRUE); break;
  //    cf->AddVariable(AliDielectronVarManager::kTRDpidQuality,"3.5, 4.5, 5.5, 6.5",kTRUE);
  //    if(!hasMC && isESD) cf->AddVariable(AliDielectronVarManager::kTRDchi2,"-1.,0.,2.,4.",kTRUE);
  */
  // mc steps
  if(hasMC) {
    //if(cutDefinition==kRec) 
    cf->SetStepForMCtruth();
    //    cf->SetStepForNoCutsMCmotherPid();
    //    cf->SetStepForAfterAllCuts();
    //    cf->SetStepsForEachCut();
    //    cf->SetStepsForSignal();
    //    cf->SetStepsForBackground(); 
    cf->SetStepsForMCtruthOnly();
  }

  die->SetCFManagerPair(cf);
}

void AddMCSignals(AliDielectron *die, Int_t cutDefinition){
  //Do we have an MC handler?
  if (!die->GetHasMC()) return;

  AliDielectronSignalMC* inclusiveJpsi = new AliDielectronSignalMC("inclusiveJpsi","Inclusive");
  inclusiveJpsi->SetLegPDGs(11,-11);
  inclusiveJpsi->SetMotherPDGs(443,443);
  inclusiveJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  inclusiveJpsi->SetFillPureMCStep(kTRUE);
  inclusiveJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  inclusiveJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);

  AliDielectronSignalMC* beautyJpsi = new AliDielectronSignalMC("beautyJpsi","Beauty");
  beautyJpsi->SetLegPDGs(11,-11);
  beautyJpsi->SetMotherPDGs(443,443);
  beautyJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  beautyJpsi->SetGrandMotherPDGs(500,500);
  beautyJpsi->SetFillPureMCStep(kTRUE);
  beautyJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);

  AliDielectronSignalMC* promptJpsi = new AliDielectronSignalMC("promptJpsi","Prompt");   // prompt J/psi (not from beauty decays)
  promptJpsi->SetLegPDGs(11,-11);
  promptJpsi->SetMotherPDGs(443,443);
  promptJpsi->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsi->SetFillPureMCStep(kTRUE);
  promptJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);

  // prompt J/psi radiative channel
  AliDielectronSignalMC* promptJpsiRad = new AliDielectronSignalMC("promptJpsiRad","PromptRadiative");   // prompt J/psi (not from beauty decays)
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

  // prompt J/psi Non radiative channel
  AliDielectronSignalMC* promptJpsiNonRad = new AliDielectronSignalMC("promptJpsiNonRad","PromptNonRadiative");   // prompt J/psi (not from beauty decays)
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

  AliDielectronSignalMC* directJpsi = new AliDielectronSignalMC("directJpsi","Direct");   // embedded J/psi
  directJpsi->SetLegPDGs(11,-11);
  directJpsi->SetMotherPDGs(443,443);
  directJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  directJpsi->SetFillPureMCStep(kTRUE);
  directJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  directJpsi->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);
  directJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  directJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);

  AliDielectronSignalMC* gammaConversion = new AliDielectronSignalMC("gammaConversion","gamma conversions");
  gammaConversion->SetLegPDGs(11,-11);
  gammaConversion->SetCheckBothChargesLegs(kTRUE,kTRUE);
  gammaConversion->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  gammaConversion->SetMotherPDGs(22,22);
  gammaConversion->SetMothersRelation(AliDielectronSignalMC::kSame);


  AliDielectronSignalMC* electrons = new AliDielectronSignalMC("electrons","electrons");
  electrons->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  electrons->SetCheckBothChargesLegs(kTRUE,kTRUE);
  electrons->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  electrons->SetFillPureMCStep(kTRUE);
  //  electrons->SetMothersRelation(AliDielectronSignalMC::kSame);

  /*
  AliDielectronSignalMC* eleHijing = new AliDielectronSignalMC("eleHijing","eleHijing");
  eleHijing->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleHijing->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleHijing->SetLegSources(AliDielectronSignalMC::kNoCocktail, AliDielectronSignalMC::kNoCocktail);
  eleHijing->SetFillPureMCStep(kTRUE);
  */
  /*
  AliDielectronSignalMC* electrons = new AliDielectronSignalMC("electrons","electrons");
  electrons->SetLegPDGs(11,-11);  //dummy second leg (never MCtrue)
  electrons->SetCheckBothChargesLegs(kTRUE,kTRUE);
  electrons->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //  electrons->SetMotherPDGs(111,111,kTRUE,kTRUE);   // not from pi0
  //  electrons->SetMothersRelation(AliDielectronSignalMC::kSame);
  electrons->SetFillPureMCStep(kTRUE);
  */

  // selection
  switch(cutDefinition) {
  case kRec:
    die->AddSignalMC(inclusiveJpsi);
    die->AddSignalMC(directJpsi);
    break;
  case kIonut:
    die->AddSignalMC(inclusiveJpsi);
    //  die->AddSignalMC(beautyJpsi);
    //die->AddSignalMC(promptJpsi);
    //die->AddSignalMC(promptJpsiRad);
    //die->AddSignalMC(promptJpsiNonRad);
    die->AddSignalMC(directJpsi);
    //  die->AddSignalMC(gammaConversion);
    break;
  case kLegEff:
    die->AddSignalMC(directJpsi);
    die->AddSignalMC(inclusiveJpsi);
    die->AddSignalMC(electrons);
    break;
  default:         return;
  }


}

void SetEtaCorrection(AliDielectron *die) {

  if (AliDielectronPID::GetCentroidCorrFunction()) return;

  TF2 *fCntrdCorr=0x0;
  TF1 *fWdthCorr=0x0;
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv DATA vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // either data or MC with tune on data option
  if( !die->GetHasMC() || ((AliPIDResponse*)AliAnalysisManager::GetAnalysisManager()->GetTasks()->At(0))->IsTunedOnData()) {
    // 2-dimensional eta correction for the centroid of electron sigmas
    fCntrdCorr = new TF2("fCntrdCorr", "[0] + [1]*y + [2]*y*y + [3]*TMath::Power(y,3) + [4]*TMath::Power(y,4) + [5]*TMath::Power(y,5) + [6]*TMath::Power(y,6) + [7]*x",
			 //			      0.0, 3000.0, -0.9, +0.9);
			 0.0, 90.0, -0.9, +0.9);
    //    fCntrdCorr->SetParameters(0.723106, 0.23958, -6.31221, -0.687976, 15.912, 0.579609, -11.6901, -0.000354381); //Nacc dep.
    fCntrdCorr->SetParameters(+0.149002, +0.214644 , -6.034930, -0.529588, +14.97902, +0.402640, -10.890027, +0.011248); //Cent dep.
    
    // 1-dimensional eta correction for the width of electron sigmas
    fWdthCorr = new TF1("fWdthCorr", "pol2", 0.0, 90.0);
    //    fWdthCorr->SetParameters(1.06108, 0.000217804,-5.80291e-08); //Nacc dep.
    fWdthCorr->SetParameters(+1.290755, -0.005261, +0.000021); //Cent dep.

    // apply corrections
    // AliDielectronPID::SetCentroidCorrFunction(fCntrdCorr,AliDielectronVarManager::kNacc,AliDielectronVarManager::kEta);
    // AliDielectronPID::SetWidthCorrFunction(fWdthCorr,AliDielectronVarManager::kNacc);
    AliDielectronPID::SetCentroidCorrFunction(fCntrdCorr,AliDielectronVarManager::kCentrality,AliDielectronVarManager::kEta);
    AliDielectronPID::SetWidthCorrFunction(fWdthCorr,AliDielectronVarManager::kCentrality);
    printf(" DATA PID correction loaded!!!\n");
  }
  else  {
    /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MONTE CARLO vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
    // 2-dimensional eta correction for the centroid of electron sigmas
    fCntrdCorr = new TF2("fCntrdCorr", "[0] + [1]*y + [2]*y*y + [3]*TMath::Power(y,3) + [4]*TMath::Power(y,4) + [5]*TMath::Power(y,5) + [6]*TMath::Power(y,6) + [7]*x",
			      0.0, 3000.0, -0.9, +0.9);
    fCntrdCorr->SetParameters(+0.378611, -0.070831, -3.076778, +0.121977, +8.576097, +0.113009, -5.001368, -0.000181);
    // 1-dimensional eta correction for the width of electron sigmas
    fWdthCorr = new TF1("fWdthCorr", "pol1", 0.0, 3000.0);
    fWdthCorr->SetParameters(+0.881894, +0.000053);

    // apply corrections
    AliDielectronPID::SetCentroidCorrFunction(fCntrdCorr,AliDielectronVarManager::kNacc,AliDielectronVarManager::kEta);
    AliDielectronPID::SetWidthCorrFunction(fWdthCorr,AliDielectronVarManager::kNacc);
    // AliDielectronPID::SetCentroidCorrFunction(fCntrdCorr,AliDielectronVarManager::kCentrality,AliDielectronVarManager::kEta);
    // AliDielectronPID::SetWidthCorrFunction(fWdthCorr,AliDielectronVarManager::kCentrality);
  }


}

TVectorD *GetRunNumbers2011() {
  
  Double_t runLHC10h[] = { // all good runs based on RCT 29.Mai
    139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161, 137135
  };
  Double_t runLHC11h[] = { // all good runs based on RCT 29.Mai
    167915, 167920, 167985, 167987, 167988, 168069, 168076, 168105, 168107, 168108, 168115, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169965, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593, 
    170593+1
  };
  /*
  if(list.Contains("LHC10h") || list.Contains("LHC11a10")) {
  if( list.Contains("LHC11h") || list.Contains("LHC12a17") ) {
  */  
  Int_t size = (int) (sizeof(runLHC11h)/sizeof(Double_t));
  TVectorD *vec = new TVectorD(size,runLHC11h);
  //vec->Print("");
  return vec;
}

TVectorD *GetDeltaPhiBins() {
  //
  // for in and out of event plane bins
  //
  Double_t pi = TMath::Pi();
  TVectorD *deltaPhi = new TVectorD(6);
  (*deltaPhi)[0] = -1.    *pi;
  (*deltaPhi)[1] = -3./4. *pi;
  (*deltaPhi)[2] = -1./4. *pi;
  (*deltaPhi)[3] = +1./4. *pi;
  (*deltaPhi)[4] = +3./4. *pi;
  (*deltaPhi)[5] = +1.    *pi;
  return deltaPhi;
}
