void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);

AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition);

TString namesDieleData=("basicQ+SPDAny+pt>1+PID;basicQ+EMCal+SPDAny+pt6_20;basicQ+EMCal+SPDAny+pt8_20;EMCal+SPDAny+pt6_20;EMCal+SPDAny+pt8_20;EMCal+SPDFirst+pt6_20;EMCal+SPDFirst+pt8_20;EMCal+SPDAny");
TObjArray *arrNamesDieleData=namesDieleData.Tokenize("; ");

const Int_t nDie=arrNamesDieleData->GetEntries();

AliDielectron* ConfigJpsi_mf_pp(Int_t cutDefinition, Bool_t isAOD=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNamesDieleData->GetEntriesFast()){
    name=arrNamesDieleData->At(cutDefinition)->GetName();
  }
  AliDielectron *diele = new AliDielectron(Form("%s",name.Data()),
                             Form("Track cuts: %s",name.Data()));

  // cut setup
  SetupEventCutsDieleFilter(diele, cutDefinition, isAOD);
  SetupTrackCutsDieleData(diele, cutDefinition, isAOD);
  SetupPairCutsDieleData(diele, cutDefinition, isAOD);
  
  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielelectron framework histograms will be filled
  //
  InitHistogramsDieleData(diele, cutDefinition, isAOD);

  // the last definition uses no cuts and only the QA histograms should be filled!
  if (cutDefinition==nDie-1) InitCFDieleData(diele, cutDefinition, isAOD);

  AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
  rot->SetConeAnglePhi(TMath::Pi());
  rot->SetIterations(10);
  diele->SetTrackRotator(rot);
  return diele;
}

//______________________________________________________________________________________
void SetupEventCutsDieleFilter(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD=kFALSE)
{
  //
  // Setup the event cuts
  //
  AliDielectronVarCuts *vtxZ = new AliDielectronVarCuts("vtxZ","Vertex z cut");
  vtxZ->AddCut(AliDielectronVarManager::kZvPrim,-10.,10.);
  
  AliDielectronVarCuts *Centrality = new AliDielectronVarCuts("Centrality","Centrality Percentile");
  Centrality->AddCut(AliDielectronVarManager::kCentrality,0.,40.);
  
  diele->GetEventFilter().AddCuts(vtxZ);
  if(cutDefinition<nDie-1) diele->GetEventFilter().AddCuts(Centrality);
  
  
}


//______________________________________________________________________________________
void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setup the track cuts
  //
  
  //ESD quality cuts DielectronTrackCuts
  if (!isAOD) {
    diele->GetTrackFilter().AddCuts(SetupESDtrackCutsDieleData(cutDefinition));
  } else {
    AliDielectronTrackCuts *trackCuts=new AliDielectronTrackCuts("trackCuts","trackCuts");
    if (cutDefinition==5||cutDefinition==6)
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    else 
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetRequireITSRefit(kTRUE);
    diele->GetTrackFilter().AddCuts(trackCuts);
  }

  //Pt cut ----------------------------------------------------------
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");
  pt->AddCut(AliDielectronVarManager::kPt,1.,1e30);
  pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
  //AOD additions since there are no AliESDtrackCuts -----------------
  pt->AddCut(AliDielectronVarManager::kNclsTPC,60.,160.);
  pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
  if(cutDefinition>2&&cutDefinition!=nDie-1) pt->AddCut(AliDielectronVarManager::kTPCsignal,70.,120.);
  if(cutDefinition==nDie-1) pt->AddCut(AliDielectronVarManager::kTPCsignal,50.,150.); 

  diele->GetTrackFilter().AddCuts(pt);

  // PID cuts --------------------------------------------------------
  AliDielectronPID *pid = new AliDielectronPID("PID10","TPC nSigma |e|<3 + |Pi|>3 + |P|>3 + TOF nSigma |e|<3");
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.,3.,0.,0.,kTRUE);
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-20.,3.,0.,0.,kTRUE);
  if(cutDefinition==0) diele->GetTrackFilter().AddCuts(pid);
  // PID cuts 2 ---- For EMCal purposes (all electron band) ----------------------------------------------------
  AliDielectronPID *pid2 = new AliDielectronPID("PIDforemcal","TPC nSigma |e|<3");
  pid2->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
  if(cutDefinition==1||cutDefinition==2) diele->GetTrackFilter().AddCuts(pid2);
  
}

//______________________________________________________________________________________
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setup the pair cuts
  //
  //Invariant mass and rapidity selection
  AliDielectronVarCuts *pairCut=new AliDielectronVarCuts("0<M<5+|Y|<.9+PtEMCalleg","0<M<5 + |Y|<.9+PtEMCalleg");
  pairCut->AddCut(AliDielectronVarManager::kM,0.,5.);
  pairCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  if(cutDefinition==1&&cutDefinition==3&&cutDefinition==5) 	pairCut->AddCut(AliDielectronVarManager::kPt,6.,20.);
  if(cutDefinition==2&&cutDefinition==4&&cutDefinition==6) 	pairCut->AddCut(AliDielectronVarManager::kPt,8.,20.);
  if(cutDefinition==nDie-1)		 			pairCut->AddCut(AliDielectronVarManager::kPt,2.,100.);
  diele->GetPairFilter().AddCuts(pairCut);
  
  if(cutDefinition>1){
	  AliDielectronVarCuts *mycut = new AliDielectronVarCuts("ptCutEMCAL","cut for EMCal");
	  mycut->AddCut(AliDielectronVarManager::kEMCALnSigmaEle,-3.,3.);
	  mycut->AddCut(AliDielectronVarManager::kP,5.,1e30);
	  AliDielectronPairLegCuts *varpair=new AliDielectronPairLegCuts();
	  varpair->GetLeg1Filter().AddCuts(mycut);
	  diele->GetPairFilter().AddCuts(varpair);
  }
  
}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition)
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;

  // basic track quality cuts  (basicQ)
  esdTrackCuts->SetMaxDCAToVertexZ(3.0);
  esdTrackCuts->SetMaxDCAToVertexXY(1.0);
  esdTrackCuts->SetEtaRange( -0.9 , 0.9 );
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetPtRange(.8,1e30);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);

  if (cutDefinition==5||cutDefinition==6)
   	esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
	else esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);  

  return esdTrackCuts;
}


//______________________________________________________________________________________
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(diele->GetName(),diele->GetTitle());
  
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
  //track rotation
  histos->AddClass(Form("Pair_%s",PairClassName(AliDielectron::kEv1PMRot)));
  histos->AddClass(Form("Track_Legs_%s",PairClassName(AliDielectron::kEv1PMRot)));
  
  
  //add histograms to event class
//   if (cutDefinition==0) {
    histos->AddClass("Event");
    histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","Centrality","Centrality;Cent(%)",100,0.,100.,AliDielectronVarManager::kCentrality);
//   }
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",400,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                        100,-1,1,144,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        200,0.2,20.,100,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_Phi","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        200,0.,2*TMath::Pi(),100,0.,200.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_nSigmaEMCal","dEdx;NsigmaEmcal;TPC signal (arb units);#tracks",
                        200,-10.,10.,100,0.,200.,AliDielectronVarManager::kEMCALnSigmaEle,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
      
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        251,-.01,5.01,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  
  diele->SetHistogramManager(histos);
}


void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(diele->GetName(),diele->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"2.0, 4.0, 6., 8.0, 10.0, 100.0");
  cf->AddVariable(AliDielectronVarManager::kY,"-5,-1,-0.9,-0.8,-0.7,-0.5,-0.3,0.3,0.5,0.7,0.8,0.9,1.0,5");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);
  cf->AddVariable(AliDielectronVarManager::kOpeningAngle,"0.,0.2,0.5,1.,2.,3.14");
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 70, 75, 80, 85, 90, 100, 120, 160",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-5,-1,-0.9,-0.85,-0.8,-0.75,0.75,0.8,0.85,0.9,1.0,5",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-2.5,-2,-1.5,-1,-0.5,4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.,3.5,4.,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.,3.5,4.,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,"60.,65.,68.,70.,80.,90.,100.,120.,150.",kTRUE);

  //event variables
  cf->AddVariable(AliDielectronVarManager::kNaccTrcklts,"0.0, 9.0, 17.0, 25.0, 36.0, 55.0, 500.0");
  cf->AddVariable(AliDielectronVarManager::kCentrality,"0.0,20.0,40.0,80.0");
  if (!isAOD){
    Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
    if (hasMC){
      cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
      cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
      cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,10000,-5000.5,4999.5,kTRUE);
    }
  }
    //only in this case write MC truth info
  if (cutDefinition==0){
    cf->SetStepForMCtruth();
  }

  diele->SetCFManagerPair(cf);
  
}

