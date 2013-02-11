void InitHistograms(AliDielectron *die, Int_t cutDefinition);

void SetupEventCuts(AliDielectron *die, ULong64_t triggers, Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts( AliDielectron *die,  Int_t cutDefinition);

void ConfigEvtPlane(AliDielectron *die,  Int_t cutDefinition);
void ConfigBgrd(    AliDielectron *die,  Int_t cutDefinition);

TString names=("NoBins;Zvtx;ZvtxNoKF;ZvtxCent;ZvtxNcontr;ZvtxNcontrPE;ZvtxNcontrPEepTPC;ZvtxNcontrPEepTPCmag");
enum { kNoBins=0, kZvtx, kZvtxNoKF, kZvtxCent, kZvtxNcontr, kZvtxNcontrPE, kZvtxNcontrPEepTPC, kZvtxNcontrPEepTPCmag};

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t  isESD = kTRUE;
Bool_t  hasMC = kFALSE;
TString list  = gSystem->Getenv("LIST");

AliDielectron* ConfigJpsiME_jpsi_PbPb(Int_t cutDefinition, TString prod="", ULong64_t triggers=AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB)
{
  //
  // Setup the instance of AliDielectron
  //

  // gsi train?
  TString trainRoot = gSystem->Getenv("TRAIN_ROOT");
  Bool_t isGSItrain = (trainRoot.IsNull()?kFALSE:kTRUE); 

  // find mc or not?
  if( list.IsNull()) list=prod;
  if( list.Contains("LHC10h")   || list.Contains("LHC11h")   ) hasMC=kFALSE;
  if( list.Contains("LHC11a10") || list.Contains("LHC12a17") ) hasMC=kTRUE;

  //ESD handler?
  isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());

  // task name
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast())  name=arrNames->At(cutDefinition)->GetName();
  printf(" Adding %s%s config %s for %s \n",(isESD?"ESD":"AOD"),(hasMC?" MC":""),name.Data(),list.Data());

  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("ME config: %s",name.Data()));
  die->SetHasMC(hasMC);

  // cut setup
  SetupEventCuts(die,triggers,cutDefinition);
  SetupTrackCuts(die,cutDefinition);
  SetupPairCuts(die,cutDefinition);

  // histogram setup
  InitHistograms(die,cutDefinition);
  printf(" Add %d class types to the histo manager \n",die->GetHistogramList()->GetEntries());

  // bgrd estimators
  ConfigBgrd(die,cutDefinition);

  // tpc event plane configuration
  ConfigEvtPlane(die,cutDefinition);

  // prefilter settings
  die->SetPreFilterUnlikeOnly();
  //die->SetPreFilterAllSigns();
  //die->SetNoPairing();

  // KF usgae
  if(cutDefinition==kZvtxNoKF) die->SetUseKF(kFALSE);

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
  if(!hasMC) {
    switch(triggers) {
    case AliVEvent::kCentral:     minCent= 0.; maxCent= 9.; break;
    case AliVEvent::kSemiCentral: minCent=12.; maxCent=53.; break;
    case AliVEvent::kMB:          minCent= 0.; maxCent=80.; break;
    default:                      minCent= 0.; maxCent=80.; break;
    }
  }

  // VZERO multiplicity vs. number ob global tracks cut
  TF1 *fMean  = new TF1("fMean", "pol1",               0,25e+3);
  fMean->SetParameters(691.633, 1.4892);
  TF1 *fSigma = new TF1("fSigma","[0]+sqrt([1]*x+[2])",0,25e+3);
  fSigma->SetParameters(-83.6599, 36.7677, 69530.7);

  // number of vertex contributors TPC vs. global cut
  TF1* vtxContribUp = new TF1("vtxContribUp","pol1",0.,20000.);
  vtxContribUp->SetParameters(0.,1.38);      //  --> strong cut, removes about 40% of events
  TF1* vtxContribLow= new TF1("vtxContribLow","pol1",0.,20000.);
  vtxContribLow->SetParameters(-100.,1.2);

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","eventCuts");
  if(!isESD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
  eventCuts->SetCentralityRange(minCent,maxCent);
  
  // apply pile-up event (PE) rejection
  switch(cutDefinition) {
  case kZvtxNcontrPE:
  case kZvtxNcontrPEepTPC:
  case kZvtxNcontrPEepTPCmag:
    eventCuts->SetCutOnV0MultipicityNTrks(fMean, fSigma, 4.0);
    //    eventCuts->SetCutOnNVtxContributorsGloablTPC(vtxContribLow, vtxContribUp);
  }
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

  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  //  trkFilter->SetMinNCrossedRowsOverFindable(0.6);
  if(!isESD) cuts->AddCut(trkFilter);

  //Pt cut, should make execution a bit faster
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("PtCut","PtCut");
  pt->AddCut(AliDielectronVarManager::kPt,1.1,1e30);    //0.8
  cuts->AddCut(pt);

  // track cuts ESD and AOD
  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9); // -0.9, 0.9
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.0);
  varCuts->AddCut(AliDielectronVarManager::kTOFbeta,      0.2,   0.9, kTRUE);
  cuts->AddCut(varCuts);
  varCuts->Print();

  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // ITS-4 = 1+2+4+8
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);
  cuts->AddCut(trkCuts);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronPID *pid = new AliDielectronPID("PID","PID");
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,4.0,0.,0.,kTRUE);
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-100.,3.5,0.,0.,kTRUE);
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);
  pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,AliDielectronPID::kIfAvailable);
  cuts->AddCut(pid);
  /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ PID CUTS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  // exclude conversion electrons selected by the tender
  AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  //  cuts->AddCut(noconv);

}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //

  // conversion rejection
  Double_t gCut = 0.05;
  AliDielectronVarCuts *gammaCuts = new AliDielectronVarCuts("GammaCuts","GammaCuts");
//  gammaCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.1,  kTRUE);
//  gammaCuts->AddCut(AliDielectronVarManager::kLegDist,      0.0,   0.25, kTRUE);
//  gammaCuts->AddCut(AliDielectronVarManager::kR,            3.0,   90.0, kTRUE);
//  gammaCuts->AddCut(AliDielectronVarManager::kPsiPair,      0.0,   0.05, kTRUE);
//  gammaCuts->AddCut(AliDielectronVarManager::kChi2NDF,      0.0,   10.0, kTRUE);
  gammaCuts->AddCut(AliDielectronVarManager::kM,            0.0,   gCut);
  die->GetPairPreFilter().AddCuts(gammaCuts);

  // rapidity selection
  //  AliDielectronVarCuts *rapCut=new AliDielectronVarCuts("|Y|<.9","|Y|<.9");
  // rapCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  // die->GetPairFilter().AddCuts(rapCut);

}

//______________________________________________________________________________________
void ConfigBgrd(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Configurate the background estimators
  //

  // add mixed events
  AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;

  // mixing handler
  switch(cutDefinition) {
  case kNoBins: /* */ break;
  case kZvtx:
  case kZvtxNoKF:
    mix->AddVariable(AliDielectronVarManager::kZvPrim, "-10.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,10.");
    break;
  case kZvtxCent:
    mix->AddVariable(AliDielectronVarManager::kZvPrim, "-10.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,10.");
    mix->AddVariable(AliDielectronVarManager::kCentrality,  8,  0.,80.);
    break;
  case kZvtxNcontrPEepTPCmag:
    mix->AddVariable(AliDielectronVarManager::kTPCmagH2,    "0.,20.,50.,80.,110.,150.,500.");
  case kZvtxNcontrPEepTPC:
    mix->AddVariable(AliDielectronVarManager::kTPCrpH2,     8,  TMath::Pi()/-2., TMath::Pi()/2.);
  case kZvtxNcontrPE:
  case kZvtxNcontr:
    mix->AddVariable(AliDielectronVarManager::kNVtxContrib, 32, 0.,3200.);
    mix->AddVariable(AliDielectronVarManager::kZvPrim, "-10.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,10.");
    break;
  }
  mix->SetMixType(AliDielectronMixingHandler::kOSonly);
  mix->SetDepth(150);
  mix->Print();

  die->SetMixingHandler(mix);
}

//______________________________________________________________________________________
void ConfigEvtPlane(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Configurate the TPC event plane 
  //

  if(cutDefinition!=kZvtxNcontrPEepTPC && cutDefinition!=kZvtxNcontrPEepTPCmag) return;

  Double_t gGap = 0.0;
  AliDielectronVarCuts *poi = new AliDielectronVarCuts("PoI","PoI");
  poi->AddCut(AliDielectronVarManager::kM,2.92,3.20);     // particles of interest, jpsi mass window
  die->GetEventPlanePOIPreFilter().AddCuts(poi);

  // eta gap in tpc event plane
  //AliDielectronVarCuts *etaGap = new AliDielectronVarCuts(AliDielectronVarManager::GetValueName(AliDielectronVarManager::kEta),"etaGap");
  //etaGap->AddCut(AliDielectronVarManager::kEta,-1*gGap,gGap,kTRUE);
  //die->GetEventPlanePreFilter().AddCuts(etaGap);
  //if(cutDefinition==kSubLS) die->SetLikeSignSubEvents();

  die->SetPreFilterEventPlane();
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //

  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());

  //add histograms to event class
  histos->AddClass("Event");


  switch(cutDefinition) {
  case kNoBins: /* */ break;
  case kZvtxCent:
    histos->UserHistogram("Event","","", 200,-10.,   10.,    AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","","", 100,  0.0, 100.0,   AliDielectronVarManager::kCentrality);
    break;
  case kZvtxNcontrPEepTPCmag:
    histos->UserHistogram("Event","","", 250,  0.,  500.,    AliDielectronVarManager::kTPCmagH2);
  case kZvtxNcontrPEepTPC:
    histos->UserHistogram("Event","","", 16, TMath::Pi()/-2., TMath::Pi()/2.,    AliDielectronVarManager::kTPCrpH2);
  case kZvtxNcontrPE:
  case kZvtxNcontr:
    histos->UserHistogram("Event","","", 200,  0., 4000.,    AliDielectronVarManager::kNVtxContrib);
  case kZvtx:
  case kZvtxNoKF:
    histos->UserHistogram("Event","","", 200,-10.,   10.,    AliDielectronVarManager::kZvPrim);
    break;
  }

  //Initialise histogram classes
  histos->SetReservedWords("Pair");
  //Pair classes inclusive mixed events
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(0)));
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(1)));
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(2)));
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));

  ///// add histograms to Pair classes /////
  histos->UserHistogram("Pair","","",  300,.0,300*0.04, 80,0.,80.,
			AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality); // 40MeV bins, 12GeV/c2
  histos->UserHistogram("Pair","","",  300,.0,300*0.04, 20,0.,20.,
			AliDielectronVarManager::kM, AliDielectronVarManager::kPt); // 40MeV bins, 12GeV/c2

  histos->UserHistogram("Pair","","",  100,-1.,1.,      AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","","",  400,0,20.,       AliDielectronVarManager::kPt);


  //legs from pair (fill SE)
  //  for (Int_t i=0; i<3; ++i){
  //  histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
  // }
  //Track classes
  //to fill also track info from 2nd event loop until 2
  //for (Int_t i=0; i<2; ++i){
  //  histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  //}


  die->SetHistogramManager(histos);
}
