void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);

AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition);

TString namesDieleData=("EMCal+SPDAny;EMCal+SPDFirst");
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
  InitCFDieleData(diele, cutDefinition, isAOD);

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
  Centrality->AddCut(AliDielectronVarManager::kCentrality,0.,80.);
  
  diele->GetEventFilter().AddCuts(vtxZ);
  diele->GetEventFilter().AddCuts(Centrality);
  
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
    if (cutDefinition==1)
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
  pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
  pt->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  pt->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  pt->AddCut(AliDielectronVarManager::kTPCsignal,68.,150.); 
  pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
  diele->GetTrackFilter().AddCuts(pt);
}

//______________________________________________________________________________________
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setup the pair cuts
  //
  //Invariant mass and rapidity selection
  AliDielectronVarCuts *pairCut=new AliDielectronVarCuts("0<M<5+|Y|<.9+PtEMCalleg","0<M<5 + |Y|<.9+PtEMCalleg");
  pairCut->AddCut(AliDielectronVarManager::kM,1.,5.);
  pairCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  pairCut->AddCut(AliDielectronVarManager::kPt,2.,100.);
  diele->GetPairFilter().AddCuts(pairCut);
  
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
  esdTrackCuts->SetPtRange(1.,1e30);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);

  if (cutDefinition==1)
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
    histos->UserHistogram("Event","Multiplicity","Multiplicity V0;Multiplicity V0",
                        500,0.,25000.);
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

  histos->UserHistogram("Track","dEdx_nSigmaEMCal","dEdx;NsigmaEmcal;TPC signal (arb units);NSigmaEMCAL",
                        200,-10.,10.,100,0.,200.,AliDielectronVarManager::kEMCALnSigmaEle,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  
  histos->UserHistogram("Track","dEdx_EoverP","dEdx;EoverP;TPC signal (arbunits);E/P",100,0.,5.,100,0.,200.,AliDielectronVarManager::kEMCALnSigmaEle,AliDielectronVarManager::kTPCsignal,kTRUE);
  

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
  cf->AddVariable(AliDielectronVarManager::kPt,"2.0, 4.0, 6.0, 8.0, 10.0, 100.0");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);
  cf->AddVariable(AliDielectronVarManager::kOpeningAngle,"0.,0.2,0.5,0.8,1.,2.,3.14");
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"1.0, 2.0, 3.0, 4.0, 5.0, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kP,"0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 10.0, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 70, 80, 90, 120, 160",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-5,-1,-0.9,-0.7,0.7,0.9,1,5",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPhi,"0.,1.,2.,3.,4.,5.,7.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEMCALnSigmaEle,"-3.,-2,-1.,1.,2.,3.,4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEMCALM02,"0.,0.1,0.2,0.4,2.,4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEMCALNCells,"0,1,2,3,4,20",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-3.,-2,-1.,1.,2.,3.,4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,"68.,70.,72.,80.,100.,120.,150.",kTRUE);
  //event variables
  cf->AddVariable(AliDielectronVarManager::kCentrality,"0.0,10.0,20.0,40.0,60.0,80.0");

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

