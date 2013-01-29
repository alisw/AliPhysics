void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);

//AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition);

TString namesDieleData=("basicQ+SPDfirst+pt>1+PID; basicQ+SPDany+pt>1+PID; basicQ+ITS012+pt>1+PID; basiQ+SPDany+pt>1+p>1.2+TPCsignalfrom70to110");
//TString namesDieleData=("basicQ+SPDfirst+pt>1+PID");

TObjArray *arrNamesDieleData=namesDieleData.Tokenize("; ");

const Int_t nDie=arrNamesDieleData->GetEntries();

AliDielectron* ConfigJpsi_mw_pPb(Int_t cutDefinition, Bool_t isAOD=kFALSE/*must be kTRUE for AODs old setting*/) 
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

  
  // estimators filename
  //NOTE: what does this mean?: estimator for pp multiplicity, not needed for instance for my pA-purpose(mwinn 16.1.2012)..
  diele->SetEstimatorFilename("$ALICE_ROOT/PWGDQ/dielectron/files/estimators.root");
  //diele->SetEstimatorFilename("estimators.root");
  // cut setup
  SetupTrackCutsDieleData(diele, cutDefinition, isAOD);
  SetupPairCutsDieleData(diele, cutDefinition, isAOD);
  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielelectron framework histograms will be filled
  //
  InitHistogramsDieleData(diele, cutDefinition, isAOD);
  
  // the last definition uses no cuts and only the QA histograms should be filled!, now for all cuts
  if(cutDefinition == 0 || cutDefinition == 1 || cutDefinition == 2 || cutDefinition == 3) InitCFDieleData(diele, cutDefinition, isAOD);

  AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
  if(cutDefinition == 0 || cutDefinition == 1 || cutDefinition == 2 || cutDefinition == 3) rot->SetConeAnglePhi(TMath::Pi());
  //else if(cutDefinition == 2) rot->SetConeAnglePhi(TMath::Pi()/180.*135.);
  //rot->SetIterations(10);
  rot->SetIterations(20);
  diele->SetTrackRotator(rot);
 
  
  return diele;
}

//______________________________________________________________________________________
void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setup the track cuts
  //
  
  if (!isAOD) {
    //exclude conversion electrons selected by the tender
    AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
    noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
    diele->GetTrackFilter().AddCuts(noconv);
    }/* else {
    //this is only for AODs
    AliDielectronTrackCuts *trackCuts=new AliDielectronTrackCuts("trackCuts","trackCuts");
    if (cutDefinition==0)
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    else if ((cutDefinition==1)||(cutDefinition ==3))
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
      trackCuts->SetRequireTPCRefit(kTRUE);
      trackCuts->SetRequireITSRefit(kTRUE);
      diele->GetTrackFilter().AddCuts(trackCuts);
      }*/
    AliDielectronTrackCuts *trackCuts=new AliDielectronTrackCuts("ITSandgeneral_trackCuts","ITSandgeneral_trackCuts");
    //ITS related cuts
    if (cutDefinition==0)
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    else if ((cutDefinition==1)||(cutDefinition ==3))
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
      trackCuts->SetRequireTPCRefit(kTRUE);
      trackCuts->SetRequireITSRefit(kTRUE);
      diele->GetTrackFilter().AddCuts(trackCuts);
  
  //Pt cut ----------------------------------------------------------
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("trackkineandTPCQ","trackkine_and_TPCQ");
  pt->AddCut(AliDielectronVarManager::kPt,1.,1e30);
  pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
  // ITS 0 1 2 : only for third variant
  if (cutDefinition==2) pt->AddCut(AliDielectronVarManager::kITSLayerFirstCls,0.,2.5);
  //AOD additions since there are no AliESDtrackCuts -----------------
  //
  // TPC #clusteres cut
  pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
  //chi square cut TOC missing...
  pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
  pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
  //TODO: DCA cuts to be investigated!!! NOTE: why?? (mwinn, 15.01.2013)
  pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
  pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  
  diele->GetTrackFilter().AddCuts(pt);
    
  // PID cuts --------------------------------------------------------
  if(cutDefinition ==0 || cutDefinition ==1 ||cutDefinition ==2){
    AliDielectronPID *pid = new AliDielectronPID("PID10","TPC nSigma |e|<3 + |Pi|>3.5 + P>3");
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.5,3.5,0.,0.,kTRUE);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-20.,3.,0.,0.,kTRUE);
    diele->GetTrackFilter().AddCuts(pid);
  }
  if(cutDefinition ==3){
    AliDielectronVarCuts *pidsubs = new AliDielectronVarCuts("pidSubs","pidsubs cut");
     pidsubs->AddCut(AliDielectronVarManager::kP,1.2,1e30);
     pidsubs->AddCut(AliDielectronVarManager::kTPCsignal,70.,110.);
     diele->GetTrackFilter().AddCuts(pidsubs);
  }
}

//______________________________________________________________________________________
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setup the pair cuts
  //
  // conversion rejection
  //Double_t gCut = 0.05;             // default
   Double_t gCut = 0.100;             // default

  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,gCut);
  diele->GetPairPreFilter().AddCuts(gammaCut);
  diele->SetPreFilterUnlikeOnly();
 
 
  //Invariant mass and rapidity selection
  AliDielectronVarCuts *pairCut=new AliDielectronVarCuts("2<M<4+|Y|<.9","2<M<4 + |Y|<.9");
  // pairCut->AddCut(AliDielectronVarManager::kM,2.,4.);
  pairCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  diele->GetPairFilter().AddCuts(pairCut);
}
/*
//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition)
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;

  // basic track quality cuts  (basicQ)
 //done esdTrackCuts->SetMaxDCAToVertexZ(3.0);
 //done esdTrackCuts->SetMaxDCAToVertexXY(1.0);

//done  esdTrackCuts->SetEtaRange( -0.9 , 0.9 );

//done esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
//done  esdTrackCuts->SetRequireITSRefit(kTRUE);
//done esdTrackCuts->SetRequireTPCRefit(kTRUE);

//done  esdTrackCuts->SetPtRange(1.,1e30);


//done  esdTrackCuts->SetMinNClustersTPC(70);

//done  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   // default SPD any
  if ((cutDefinition==1)||(cutDefinition==3)) esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);

  if (cutDefinition==0)
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

  return esdTrackCuts;
}
*/

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
  if (cutDefinition==0) {
    histos->AddClass("Event");
    histos->AddClass("Event_noCuts"); 
   histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    // nAcc
    histos->UserHistogram("Event","NAccRaw","Accepted raw SPD tracklets, |y|<1; nTrackl; #Entries",101,-0.5,100.5, AliDielectronVarManager::kNaccTrckltsEsd10); 
    histos->UserHistogram("Event","NAccCorr","Accepted corr SPD tracklets, |y|<1; nTrackl; #Entries",101,-0.5,100.5, AliDielectronVarManager::kNaccTrckltsEsd10Corr); 
   // nAcc vs Zvtx
   histos->UserHistogram("Event","NAccRaw_vs_Zvtx","Accepted raw SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,101,-0.5,100.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10);
    histos->UserHistogram("Event","NAccCorr_vs_Zvtx","Accepted corr SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,101,-0.5,100.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10Corr);

 
   // no event cuts 
   histos->UserHistogram("Event_noCuts","VtxZ","Vertex Z;Z[cm]",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    // nAcc
    histos->UserHistogram("Event_noCuts","NAccRaw","Accepted raw SPD tracklets, |y|<1; nTrackl; #Entries",101,-0.5,100.5, AliDielectronVarManager::kNaccTrckltsEsd10);
    histos->UserHistogram("Event_noCuts","NAccCorr","Accepted corr SPD tracklets, |y|<1; nTrackl; #Entries",101,-0.5,100.5, AliDielectronVarManager::kNaccTrckltsEsd10Corr);
   // nAcc vs Zvtx
   histos->UserHistogram("Event_noCuts","NAccRaw_vs_Zvtx","Accepted raw SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,101,-0.5,100.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10);
    histos->UserHistogram("Event_noCuts","NAccCorr_vs_Zvtx","Accepted corr SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,101,-0.5,100.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10Corr);

   }
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",400,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                        100,-1,1,144,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        200,0.2,20.,100,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
      
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        125,0.,125*.04,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
 

 // 3D histos: invMass - Multiplicity - ptJpsi
  histos->UserHistogram("Pair","InvMass_NaccRaw_PtJpsi","Inv.Mass - NaccRaw - PtJpsi;Inv. Mass [GeV];NaccRaw; pTJpsi[GeV/c]", 125,0.,125*.04,101,-0.5,100.5, 100, 0.,10., AliDielectronVarManager::kM,AliDielectronVarManager::kNaccTrckltsEsd10, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","InvMass_NaccCorr_PtJpsi","Inv.Mass - NaccCorr - PtJpsi;Inv. Mass [GeV];NaccCor; pTJpsi[GeV/c]", 125,0.,125*.04,101,-0.5,100.5,100,0.,10., AliDielectronVarManager::kM,AliDielectronVarManager::kNaccTrckltsEsd10Corr, AliDielectronVarManager::kPt);
 
  diele->SetHistogramManager(histos);
}


void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(diele->GetName(),diele->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 1.0, 1.3, 2.0, 3.0, 5., 7.0, 10.0, 100.0");
  cf->AddVariable(AliDielectronVarManager::kY,"-5,-1,-0.9,-0.8,-0.7,-0.5,-0.3,0.3,0.5,0.7,0.8,0.9,1.0,5");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPseudoProperTime,150,-0.3,0.3);
  cf->AddVariable(AliDielectronVarManager::kPseudoProperTimeErr,600,0.,0.3);
  cf->AddVariable(AliDielectronVarManager::kPairType,"-0.5,0.5,1.5,2.5,9.5,10.5");
  

  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 70, 75, 80, 85, 90, 100, 120, 160",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCchi2Cl,"0, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.",kTRUE);  
  cf->AddVariable(AliDielectronVarManager::kEta,"-5,-1,-0.9,-0.85,-0.8,-0.75,0.75,0.8,0.85,0.9,1.0,5",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-2.5,-2,-1.5,-1,-0.5,4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.,3.5,4.,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.,3.5,4.,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,-0.5,5.5,kTRUE);

  //event variables
  //cf->AddVariable(AliDielectronVarManager::kNaccTrcklts,"0.0, 9.0, 17.0, 25.0, 36.0, 55.0, 500.0");
  cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd10,101,-0.5,100.5);
  cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd10Corr,101,-0.5,100.5);
  cf->AddVariable(AliDielectronVarManager::kZvPrim,"-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.");
  if (!isAOD){
    Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
    if (hasMC){
      cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
      cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
      cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,10000,-5000.5,4999.5,kTRUE);
    }
  }
  
  //only in this case write MC truth info
  //if (cutDefinition==0){
  //  cf->SetStepForMCtruth();
  //}

  diele->SetCFManagerPair(cf);
  
}

