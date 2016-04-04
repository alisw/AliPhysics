void InitHistograms(AliDielectron *die, Int_t cutDefinition);
//void InitCF(AliDielectron* die, Int_t cutDefinition);
//void AddHistsEleEff(AliDielectron *die);
void SetupEventCuts(AliDielectron *die, Int_t cutDefinition);
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);
AliESDtrackCuts *SetupPreFilterESDtrackCuts(Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die,  Int_t cutDefinition);
AliDielectronPID *SetPIDcuts(Int_t cutDefinition);
AliDielectronPID *SetupPreFilterPIDcuts(Int_t cutDefinition);

void SetupV0cuts(AliDielectron *die, Int_t cutDefinition);
TVectorD *GetRunNumbers();

void SetupMCsignals(AliDielectron *die);

//Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
Bool_t kRot=1;
Bool_t kMix=1;

TVectorD *GetRunNumbers() {
  Double_t first=0;
  Double_t last =1;
		
  switch(iPeriod) {   
  case k10h: first=136831; last=139517; break;
  case k11h: first=165772; last=170718; break;
  case k15o: first=224917; last=226433; break;

  }
  return (AliDielectronHelper::MakeLinBinning(last-first, first, last));
}

//
// Here the configuration part starts
//


void Config(AliAnalysisTaskMultiDielectron *task){
  vector<unsigned long> cutsToUse;
  cutsToUse.push_back(0);
  cutsToUse.push_back(1);
  cutsToUse.push_back(2);
  TString names=
    "JPsi_Default;"
    "JPsi_OneSoftOneStrong;"
    "JPsi_TwoStrong;";
  TObjArray *arrNames=names.Tokenize(";");
  
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetCentralityRange(0.0,50.0);
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
  InitHistograms(die,cutDefinition);

  if(kMix){
    //  #### Event Mixing Handler #####
    AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
    mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-8,-5,0,5,8,10");
    //   mix->AddVariable(AliDielectronVarManager::kPhi ,"0,3.14,6.3");
    mix->SetDepth(10);
    die->SetMixingHandler(mix);
  } //track rotator
 if(kRot){
   AliDielectronTrackRotator *rot = new AliDielectronTrackRotator;
   rot->SetConeAnglePhi(TMath::Pi()/180.*135.);
   rot->SetIterations(10);
   die->SetTrackRotator(rot);
 }   
 task->AddDielectron(die);
  }
}

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
  gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,TMath::Cos(0.02),1.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF, 0.0,  10.0, kFALSE);//to be checked, if properly filled
  gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist, 0.0,   0.25, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kR, 3.0,  90.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair, 0.0,   0.05, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kM, 0.0,   0.05, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.1, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt, 0.0,   0.05, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha, -0.35,0.35, kFALSE); // not sure if it works as expected
  gammaV0Cuts->SetExcludeTracks(kTRUE);//ktrue excludes tracks v0s, 
  //kfalse 
  // gammaV0Cuts->Print();
  //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
  //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; && const Double_t cutQTG2 < 0.04;
  if(cutDefinition==0){
  die->GetTrackFilter().AddCuts(gammaV0Cuts);
  }
  //gammaV0Cuts->Print();	// 
}


//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  // selection or rejection of V0 tracks                                                
  //  noconv->SetExcludeTracks(kTRUE);
  //die->SetUseKF(kFALSE);
  
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);
  
  //_____________________________________________PREFILTER_________________________________________________
  if(cutDefinition==1 || cutDefinition==2){
  
    //apply soft cuts all tracks
    die->GetTrackFilter().AddCuts(SetupPreFilterESDtrackCuts(cutDefinition));
    die->GetTrackFilter().AddCuts(SetupPreFilterPIDcuts(cutDefinition));
    // die->GetPairPreFilterLegs().AddCuts(noconv);

    //apply strong cut on one leg and soft cut on the other leg 
    AliDielectronPairLegCuts *strongCuts=new AliDielectronPairLegCuts("StrongTrackCuts","StrongTrackCuts");  
    strongCuts->GetLeg1Filter().AddCuts(SetPIDcuts(cutDefinition));
    strongCuts->GetLeg1Filter().AddCuts(SetupESDtrackCuts(cutDefinition));
    strongCuts->GetLeg2Filter().AddCuts(SetPIDcuts(cutDefinition));
    strongCuts->GetLeg2Filter().AddCuts(SetupESDtrackCuts(cutDefinition));

    if(cutDefinition==1){
      strongCuts->SetCutType(AliDielectronPairLegCuts::kAnyLeg);
    }   else{
      strongCuts->SetCutType(AliDielectronPairLegCuts::kBothLegs);
    }
    
    //pair Prefilter cuts
    AliDielectronVarCuts* pairCutsInvM  = new AliDielectronVarCuts("pairCutsInvM","pairCutsInvM");
    AliDielectronVarCuts* pairCutsOpAng = new AliDielectronVarCuts("pairCutsOpAng","pairCutsOpAng");
    AliDielectronVarCuts *rapidityCut=new AliDielectronVarCuts("rapidityCut","rapidityCut");


    pairCutsInvM ->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
    pairCutsOpAng->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.050);
    rapidityCut->AddCut(AliDielectronVarManager::kY,-.9, .9 );

    
    AliDielectronCutGroup* pairCutsPF =new AliDielectronCutGroup("pairCutsPF","pairCutsPF",AliDielectronCutGroup::kCompAND);
    pairCutsPF->AddCut(pairCutsInvM);
    pairCutsPF->AddCut(pairCutsOpAng);
    die->GetPairPreFilter().AddCuts(pairCutsPF);
    die->GetPairPreFilter().AddCuts(strongCuts);
    
    die->SetPreFilterUnlikeOnly(kTRUE);
   //FinalTrackCuts after prefiltering
    //AliDielectronPairLegCuts *strongFCuts=new AliDielectronPairLegCuts("StrongFLegCuts","StrongfLegTrackCuts");  
    AliDielectronCutGroup* strongFcuts = new AliDielectronCutGroup("strongFcuts","strongFcuts",AliDielectronCutGroup::kCompAND);  
    strongFcuts->AddCut(SetPIDcuts(cutDefinition));
    strongFcuts->AddCut(SetupESDtrackCuts(cutDefinition));
    die->GetPairPreFilterLegs().AddCuts(strongFcuts);
    die->GetPairFilter().AddCuts(rapidityCut); 

    
  }//std way to cut adding the v0 rejection
  else if(cutDefinition==0) {
    die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));
    die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));
    
    // add conversion rejection
    AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
    gammaCut->AddCut(AliDielectronVarManager::kM,0.,.05);
    die->GetPairPreFilter().AddCuts(gammaCut);
    
  } 
}
//-------------------------------prefilter track cuts-----------------------------------
AliESDtrackCuts *SetupPreFilterESDtrackCuts(Int_t cutDefinition){
  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts;
  //global
  fesdTrackCuts->SetPtRange(0.3, 100. );
  fesdTrackCuts->SetEtaRange( -1.1 , 1.1 );
 // fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fesdTrackCuts->SetDCAToVertex2D(kFALSE);
 // fesdTrackCuts->SetMaxDCAToVertexZ(3.);
 // fesdTrackCuts->SetMaxDCAToVertexXY(1.);
  //ITS cuts
 // fesdTrackCuts->SetRequireITSRefit(kTRUE);
  fesdTrackCuts->SetRequireTPCRefit(kTRUE);
  fesdTrackCuts->SetMinNClustersTPC(50);
  return fesdTrackCuts;
}

//pidcuts prefilter
AliDielectronPID *SetupPreFilterPIDcuts(Int_t cutDefinition){
  AliDielectronPID *pidits = new AliDielectronPID();
  pidits->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.0,5.0,0.0, 100., kFALSE ,AliDielectronPID::kRequire,-1);
  //pidits->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,2.0 ,0.0, 100., kTRUE ,AliDielectronPID::kRequire,-1);
  //pidits->AddCut(AliDielectronPID::kTPC,AliPID::kProton,  -100. ,2.0 ,0.0, 100., kTRUE ,AliDielectronPID::kRequire,-1); 
  return pidits;
}
//-----------------------------------track cuts-----------------------------------------                                                                                             
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition){
  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts;
  //global                                                                                                                                            
  fesdTrackCuts->SetPtRange(0.85, 100. );
  fesdTrackCuts->SetEtaRange( -0.9 , 0.9 );
  fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCuts->SetDCAToVertex2D(kFALSE);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.);
  fesdTrackCuts->SetMaxDCAToVertexXY(1.);
  fesdTrackCuts->SetRequireTPCRefit(kTRUE);
  fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  fesdTrackCuts->SetMinNClustersTPC(70);
  // fesdTrackCuts->SetMinNCrossedRowsTPC(100);
  //  fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
  fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  return fesdTrackCuts;
}


//____________________PID cuts___________________________________________________

AliDielectronPID *SetPIDcuts(Int_t cutDefinition){
  AliDielectronPID *pid = new AliDielectronPID();
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,      -100. ,3.5 ,0.0, 100., kTRUE ,AliDielectronPID::kRequire,-1);
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -100. ,3.5 ,0.0, 100., kTRUE ,AliDielectronPID::kRequire,-1);  
//  pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
  // pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,    -2.,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,-1);
  return pid;
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
    // to fill also mixed event histograms loop until 10
  for (Int_t i=0; i<12; ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
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
  
      if(cutDefinition==0){
    
    //add histograms to event class
    histos->AddClass("Event");
    histos->UserHistogram("Event","","",100,0.,100.,AliDielectronVarManager::kCentrality);
    histos->UserHistogram("Event","","",100,0.,10.,AliDielectronVarManager::kPairs);
    
    histos->UserHistogram("Event","","", AliDielectronHelper::MakeLinBinning(90, 0., 90.), AliDielectronHelper::MakeLinBinning(3000, 0., 3000.),
			  AliDielectronVarManager::kCentrality, AliDielectronVarManager::kRefMultTPConly);
    histos->UserHistogram("Event","","",700,0.,700.,AliDielectronVarManager::kNTrk);
    histos->UserHistogram("Event","","",100,0.,10.,AliDielectronVarManager::kPairs);
    histos->UserHistogram("Event","","",4,0.,4.,AliDielectronVarManager::kNevents);    
    histos->UserHistogram("Event","","",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","","",100,-2.,2.,100,-2.,2.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);
    histos->UserProfile("Event","","",AliDielectronVarManager::kNTrk, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","","",AliDielectronVarManager::kPairs, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Event","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(700,0.,700),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNTrk);
    histos->UserHistogram("Event","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(10000,0.,10),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kPairs);
  }
  
  histos->UserHistogram("Track","","",200,0.2,20.,AliDielectronVarManager::kPt); 

  histos->UserHistogram("Track","","",200,0.,20.,AliDielectronVarManager::kP);

  histos->UserHistogram("Track","","",80,-1.0,1.0,AliDielectronVarManager::kEta); 

  histos->UserHistogram("Track","","",400,0.0,4.0,AliDielectronVarManager::kPhi);  

  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE); 
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

  histos->UserHistogram("Track","","",
  			160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
 
  histos->UserHistogram("Track","","",
                        160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  histos->UserHistogram("Track","","",
			100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(150,0.,150.),AliDielectronVarManager::kNaccTrckltsEsd10Corr);
  histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(5,-0.5,2.0),AliDielectronVarManager::kTOFPIDBit);

  histos->UserHistogram("Track","","",200,0.,20.,161,-0.5,161.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);

  histos->UserHistogram("Track","","",			
			200,0.0,20.0,160,0,1.1,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCrFrac);
  
  histos->UserHistogram("Track","","",10000,0.0,1.0,AliDielectronVarManager::kM);

  histos->UserHistogram("Track","","",500,-1.,1.,AliDielectronVarManager::kImpactParXY);

  histos->UserHistogram("Track","","",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
   
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

