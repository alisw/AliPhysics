void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);
//void AddHistsEleEff(AliDielectron *die);
void SetupEventCuts(AliDielectron *die, Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die,  Int_t cutDefinition);
void SetupV0cuts(AliDielectron *die, Int_t cutDefinition);
TVectorD *GetRunNumbers();
void SetupMCsignals(AliDielectron *die);
//Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
Bool_t kRot=0;
Bool_t kMix=0;
TVectorD *GetRunNumbers() {
  Double_t first=0;
  Double_t last =1;
  switch(iPeriod) {
  case k10h: first=136831; last=139517; break;
  case k11h: first=165772; last=170718; break;
  case k15o: first=244916; last=246433; break;
  }
  return (AliDielectronHelper::MakeLinBinning(last-first, first, last));
}
// Here the configuration part starts
void Config(AliAnalysisTaskMultiDielectron *task){
  vector<unsigned long> cutsToUse;
  cutsToUse.push_back(0);
  cutsToUse.push_back(1);
  cutsToUse.push_back(2);
  cutsToUse.push_back(3);
  TString names="JPsi_Any;"
    "JPsi_First;"
    "JPsi_Conversions;"
    "Jpsi_ActLength;";
  TObjArray *arrNames=names.Tokenize(";");
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  //  eventCuts->SetCentralityRange(0.0,90.0,kTRUE);
    //  task->SetTriggerOnV0AND();
  task->SetEventFilter(eventCuts);

 for(Int_t cutNumber=0; cutNumber < cutsToUse.size(); cutNumber++){
    Int_t cutDefinition = cutsToUse[cutNumber];
    TString name=Form("%02d",cutDefinition);
    if (cutNumber < arrNames->GetEntriesFast()){
      name=arrNames->At(cutNumber)->GetName();
    }
  AliDielectron *die =new AliDielectron(Form("%s",name.Data()),Form("Track cuts: %s",name.Data()));

  //
  // Manager *mgr = AliAnalysisManager::GetAnalysisManager();
  //Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  //  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  if (hasMC)SetupMCsignals(die);
  SetupTrackCuts(die,cutDefinition);
  if(cutDefinition==0 || cutDefinition==1 || cutDefinition==3){
  SetupPairCuts(die,cutDefinition);
  }
  InitHistograms(die,cutDefinition);

  if(kMix){
    //  #### Event Mixing Handler #####

    if(cutDefinition==0 || cutDefinition==1  ||  cutDefinition==3 ){
    AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
    mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
    //   mix->AddVariable(AliDielectronVarManager::kPhi ,"0,3.14,6.3");
    mix->SetDepth(10);
    die->SetMixingHandler(mix);
    }
  } //track rotator
 if(kRot){
   if(cutDefinition==0|| cutDefinition==1 ||  cutDefinition==3  ){
     AliDielectronTrackRotator *rot = new AliDielectronTrackRotator;
     rot->SetConeAnglePhi(TMath::Pi()/180.*135.);
     rot->SetIterations(10);
     die->SetTrackRotator(rot);
   }
 }
 if(cutDefinition==2){
   die->SetNoPairing();
   }
 // die->SetNoPairing();
  die->SetUseKF(kFALSE);
 task->AddDielectron(die);
 }}

//______________________________________________________________________________________
void SetupEventCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the event cuts
  //
  // trigger specific centrality cuts (reject trigger inefficiencies)
  /* Double_t minCent=0.0, maxCent=100.;
  if(!die->GetHasMC()) {
    switch(triggers) {
    case AliVEvent::kCentral:     minCent= 0.; maxCent=10.; break; //0-9
    case AliVEvent::kSemiCentral: minCent=10.; maxCent=50.; break; //12-53
    case AliVEvent::kMB:          minCent=10.; maxCent=90.; break;
    default:                      minCent= 0.; maxCent=90.; break;
    }
  }
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","eventCuts");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetCentralityRange(0.0,50.0);
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
  //eventCuts->Print();
  die->GetEventFilter().AddCuts(eventCuts);*/
}

void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the V0 cuts
  //
  AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");
  gammaV0Cuts->SetPdgCodes(22,11,11);
  gammaV0Cuts->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
  gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,TMath::Cos(0.05),1.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF, 0.0,  10.0, kFALSE);//to be checked, if properly filled
  gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist, 0.0,   0.25, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kR, 3.0,  90.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kPhivPair, 0.0,   0.02, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kM, 0.0,   0.1, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.1, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt, 0.0,   0.05, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha, -0.35,0.35, kFALSE); // not sure if it works as expected
  if(cutDefinition==2){
  gammaV0Cuts->SetExcludeTracks(kFALSE);//ktrue excludes tracks v0s,
  }else{
    gammaV0Cuts->SetExcludeTracks(kTRUE);//ktrue excludes tracks v0s,
  }
  //kfalse
  // gammaV0Cuts->Print();
  //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
  //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; && const Double_t cutQTG2 < 0.04;
  die->GetTrackFilter().AddCuts(gammaV0Cuts);
  //gammaV0Cuts->Print();	//
}


//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);

  //default quality cuts
  AliDielectronTrackCuts *refit=new AliDielectronTrackCuts("refit","refit");

  if (cutDefinition==0){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }else if (cutDefinition==1){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  }else if (cutDefinition==2){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  } else if (cutDefinition==3){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }




  cuts->AddCut(refit);

  //pt and kink mother
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");

  if (cutDefinition==0){
    pt->AddCut(AliDielectronVarManager::kPt,0.85,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    //impact parameter
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  } else if (cutDefinition==1){
    pt->AddCut(AliDielectronVarManager::kPt,0.85,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    //impact parameter
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }else if (cutDefinition==2){
    pt->AddCut(AliDielectronVarManager::kPt,0.85,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    //impact parameter
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }else if (cutDefinition==3){
    pt->AddCut(AliDielectronVarManager::kPt,0.85,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    pt->AddCut(AliDielectronVarManager::kTPCGeomLength,1.02,2.);
    //impact parameter
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }

  cuts->AddCut(pt);
}


//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //

  // add conversion rejection
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,.05);
  die->GetPairPreFilter().AddCuts(gammaCut);

  AliDielectronVarCuts *rapidityCut=new AliDielectronVarCuts("rapidityCut","rapidityCut");
  rapidityCut->AddCut(AliDielectronVarManager::kY,-.9, .9 );
  die->GetPairFilter().AddCuts(rapidityCut);
  die->SetPreFilterUnlikeOnly(kTRUE);
}


//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //

  //Setup histogram Manager
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

  Int_t pairClasses [5] = { AliDielectron::kEv1PP,  AliDielectron::kEv1PM, AliDielectron::kEv1MM, AliDielectron::kEv1MEv2P, AliDielectron::kEv1PMRot};
if(cutDefinition<2){
  // to fill also mixed event histograms loop until 10
  for (Int_t i=0; i<5; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(pairClasses[i])));
  }

  //add MC signal histograms to track and pair class
  if(die->GetMCSignals()) {
    for(Int_t isig=0; isig<die->GetMCSignals()->GetEntriesFast(); isig++) {
      TString sigMCname = die->GetMCSignals()->At(isig)->GetName();

      // mc truth
      histos->AddClass(Form("Pair_%s_MCtruth",       sigMCname.Data()));
      histos->AddClass(Form("Track_Legs_%s_MCtruth", sigMCname.Data()));
      // mc reconstructed
      histos->AddClass(Form("Pair_%s",               sigMCname.Data()));
      histos->AddClass(Form("Track_Legs_%s",         sigMCname.Data()));
    }
  }
}
  //add histograms to event class
  histos->AddClass("Event");
  histos->UserHistogram("Event","","",100,0.,100.,AliDielectronVarManager::kCentralityNew);
  histos->UserHistogram("Event","","",50,0.,50.,AliDielectronVarManager::kPairs);
  histos->UserHistogram("Event","","",3000,0.,3000.,AliDielectronVarManager::kNTrk);
  histos->UserHistogram("Event","","",4,0.,4.,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","","",300,-15.,15.,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","","",100,-2.,2.,100,-2.,2.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);
/*  histos->UserHistogram("Event","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(700,0.,700),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNTrk);
  histos->UserHistogram("Event","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(10000,0.,10),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kPairs);

  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTOFnSigmaEle);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(200,0.,200.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(150,0.,150.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(200,0.,200.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kLeg1DCAabsXY);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(20,-1.0,1.0),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(63,0.,6.32),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(200,0.,200.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNch10);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(200,0.,200.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNaccItsPureEsd10);
*/
  histos->UserHistogram("Track","","",200,0.2,20.,AliDielectronVarManager::kPt);

  histos->UserHistogram("Track","","",200,0.,20.,AliDielectronVarManager::kP);

  histos->UserHistogram("Track","","",200,-1.0,1.0,AliDielectronVarManager::kEta);

  histos->UserHistogram("Track","","",400,0.0,4.0,AliDielectronVarManager::kPhi);



    histos->UserHistogram("Track","","",
  		200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
    histos->UserHistogram("Track","","",
  		200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
    histos->UserHistogram("Track","","",
  		200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
    histos->UserHistogram("Track","","",
  			200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
    histos->UserHistogram("Track","","",
      			200,-1.0,1.0,100,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
    histos->UserHistogram("Track","","",
      			200,-1.0,1.0,100,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
    histos->UserHistogram("Track","","",
      			200,-1.0,1.0.,100,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
    histos->UserHistogram("Track","","",
      			200,-1.0,1.0,100,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
    histos->UserHistogram("Track","","",
      			100,0.0,100,100,-10.,10.,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
    histos->UserHistogram("Track","","",
      			100,0.0,100,100,-10.,10.,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
    histos->UserHistogram("Track","","",
      			100,0.0,100.0.,100,-10.,10.,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
    histos->UserHistogram("Track","","",
      			100,0.0,100.0,100,-10.,10.,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
    histos->UserHistogram("Track","","",
        		200,0.1,20.,100,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,kTRUE);
    histos->UserHistogram("Track","","",
			  200,0.1,20.,100,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaEle,kTRUE); // add by xiaozhi
    histos->UserHistogram("Track","","",
			  200,0.,4.,100,-10.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle,kTRUE); // add by xiaozhi


    histos->UserHistogram("Track","","",
            200,0.,200,AliDielectronVarManager::kNch10,kTRUE);
    histos->UserHistogram("Track","","",
  			200,0.2,20.,100,0.,1.2,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFbeta,kTRUE);
    histos->UserHistogram("Track","","",
  			200,0.2,20.,200,-20.,20.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
    histos->UserHistogram("Track","","",
  			100,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","","",
  			144,0.0,6.285,200,0.0,200,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignal);
     histos->UserHistogram("Track","","",
  			144,0.0,6.285,40,-1.0,1.0,AliDielectronVarManager::kPhi,AliDielectronVarManager::kEta);
    histos->UserHistogram("Track","","",
  			40,-1.0,1.0,100,0.0,200,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);

    histos->UserHistogram("Track","","",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","","",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
    histos->UserHistogram("Track","","",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
    histos->UserHistogram("Track","","",100,0.,3.0,AliDielectronVarManager::kTPCGeomLength);
    histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(150,0.,150.),AliDielectronVarManager::kNaccTrckltsEsd10Corr);

  histos->UserHistogram("Track","","",200,0.,20.,161,-0.5,161.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","","",200,0.0,20.0,160,0,1.1,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCrFrac);

  histos->UserHistogram("Track","","",10000,0.0,1.0,AliDielectronVarManager::kM);

  histos->UserHistogram("Track","","",500,-1.,1.,AliDielectronVarManager::kImpactParXY);

  histos->UserHistogram("Track","","",600,-3.,3.,AliDielectronVarManager::kImpactParZ);

  histos->UserHistogram("Track","","",7,0.,7.,AliDielectronVarManager::kNclsITS);

  histos->UserHistogram("Track","","",10000,0.0,1.0,AliDielectronVarManager::kM);

  histos->UserHistogram("Track","","",500,-1.,1.,AliDielectronVarManager::kImpactParXY);

  histos->UserHistogram("Track","","",600,-3.,3.,AliDielectronVarManager::kImpactParZ);

  histos->UserHistogram("Track","","",7,0.,7.,AliDielectronVarManager::kNclsITS);

//TRD part new

histos->UserHistogram("Track","","",300,0.5,300.5,AliDielectronVarManager::kNclsTRD);

histos->UserHistogram("Track","","",200,0.0,20.0,300,0.5,300.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kNclsTRD);

histos->UserHistogram("Track","","",10,0.,10.,AliDielectronVarManager::kTRDntracklets);

//histos->UserHistogram("Track","","",10,0.,10.,AliDielectronVarManager::kTRDpidQuality);

histos->UserHistogram("Track","","",300,0.,1000.,AliDielectronVarManager::kTRDchi2);

histos->UserHistogram("Track","","",20,0.,20.,AliDielectronVarManager::kTRDchi2Trklt);

histos->UserHistogram("Track","","",100,0.0,1.0,AliDielectronVarManager::kTRDpidEffLeg);

histos->UserHistogram("Track","","",
    20,-1.0,1.0,70,-3.5,3.5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTRDphi);

histos->UserHistogram("Track","","",100,0,1000,AliDielectronVarManager::kTRDsignal);

histos->UserHistogram("Track","","",20,-1.0,1.0,AliDielectronVarManager::kInTRDacceptance);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprobEle);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprobPio);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprob2DEle);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprob2DPio);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprob2DPro);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprob3DEle);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprob3DPio);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprob3DPio);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprob7DEle);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprob7DPio);

histos->UserHistogram("Track","","",20,0.0,1.0,AliDielectronVarManager::kTRDprob7DPro);



  //add histograms to Pair classes
    histos->UserHistogram("Pair","","",
			  301,-.01,6.01,AliDielectronVarManager::kM);
    histos->UserHistogram("Pair","","",
			  100,-1.,1.,AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","","",
			  100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","","",
			  301,-.01,6.01,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","","",
			  301,-.01,6.01,500,0.,5.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);
     histos->UserHistogram("Pair","","",
		       100,-3.,3.,AliDielectronVarManager::kLegDist);

     histos->UserHistogram("Pair","","",
			   100,-3.,3.,AliDielectronVarManager::kLegDistXY);

     histos->UserHistogram("Pair","","",
			   100,-3.,3.,AliDielectronVarManager::kLeg1DCAsigXY);
    /*
    //add MC signal histograms to track class
    if(die->GetMCSignals()) {

      TString className="";
      for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) {
	TString sigMCname = die->GetMCSignals()->At(i)->GetName();

	// mc truth
	histos->AddClass(Form("Pair_%s_MCtruth",      sigMCname.Data()));
	histos->AddClass(Form("Track_Legs_%s_MCtruth",sigMCname.Data()));
	// mc reconstructed
	histos->AddClass(Form("Pair_%s",              sigMCname.Data()));
	histos->AddClass(Form("Track_Legs_%s",        sigMCname.Data()));
	// tracks
	histos->AddClass(Form("Track_%s_%s",AliDielectron::PairClassName(AliDielectron::kEv1PM),sigMCname.Data()));
	histos->AddClass(Form("Track_%s_%s_MCtruth",AliDielectron::PairClassName(AliDielectron::kEv1PM),sigMCname.Data()));
      } //end: loop signals
    } //end: has signals
    */
    // add single electron histograms
    //AddHistsEleEff(die);


    //  }
  die->SetHistogramManager(histos);
}

//______________________________________________________________________________________
void SetupMCsignals(AliDielectron *die){
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
  /*
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
  die->AddSignalMC(promptJpsiNonRad);

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
  die->AddSignalMC(promptJpsiRad);
  */
}
