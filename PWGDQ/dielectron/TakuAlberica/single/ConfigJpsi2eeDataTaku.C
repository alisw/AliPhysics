void SetupTrackCutsDieleData(AliDielectronTaku *diele, Int_t cutDefinition, Bool_t isAOD);
void SetupPairCutsDieleData(AliDielectronTaku *diele, Int_t cutDefinition, Bool_t isAOD);
void InitHistogramsDieleData(AliDielectronTaku *diele, Int_t cutDefinition, Bool_t isAOD);
void InitTreesDieleData(AliDielectronTaku *diele, Int_t cutDefinition, Bool_t isAOD);
void InitCFDieleData(AliDielectronTaku *diele, Int_t cutDefinition, Bool_t isAOD);

AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition);

//TString namesDieleData=("basicQ+SPDfirst+pt>1+PID; basicQ+SPDany+pt>1+PID");
TString namesDieleData=("basicQ+SPDfirst+pt>1+PID");

TObjArray *arrNamesDieleData=namesDieleData.Tokenize("; ");

const Int_t nDie=arrNamesDieleData->GetEntries();

AliDielectronTaku* ConfigJpsi2ee(Int_t cutDefinition, Bool_t isAOD=kFALSE)
{
  //
  //cout<<" Setup the instance of AliDielectronTaku "<<endl;
  //
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNamesDieleData->GetEntriesFast()){
    name=arrNamesDieleData->At(cutDefinition)->GetName();
  }
  AliDielectronTaku *diele = new AliDielectronTaku(Form("%s",name.Data()),
						   Form("Track cuts: %s",name.Data()));

  // cut setup
  SetupTrackCutsDieleData(diele, cutDefinition, isAOD);
  SetupPairCutsDieleData(diele, cutDefinition, isAOD);
  
  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielelectron framework histograms will be filled
  //
  InitHistogramsDieleData(diele, cutDefinition, isAOD);
  //InitTreesDieleData(diele, cutDefinition, isAOD);

  // the last definition uses no cuts and only the QA histograms should be filled!
//   if (cutDefinition<nDie-1)
  InitCFDieleData(diele, cutDefinition, isAOD);

  AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
  rot->SetConeAnglePhi(TMath::Pi());
  rot->SetIterations(10);
  diele->SetTrackRotator(rot);
  return diele;
}

//______________________________________________________________________________________
void SetupTrackCutsDieleData(AliDielectronTaku *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setup the track cuts
  //
  cout<<" SetupTrackCutsDieleData "<<endl;
  //ESD quality cuts DielectronTrackCuts
  if (!isAOD) {
    diele->GetTrackFilter().AddCuts(SetupESDtrackCutsDieleData(cutDefinition));
  } else {
    AliDielectronTrackCuts *trackCuts=new AliDielectronTrackCuts("trackCuts","trackCuts");
    if (cutDefinition==0)
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    else if (cutDefinition==1)
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetRequireITSRefit(kTRUE);
    diele->GetTrackFilter().AddCuts(trackCuts);
  }

   //Pt cut ----------------------------------------------------------

  AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");
  pt->AddCut(AliDielectronVarManager::kPt,0.2,10.);
  pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
  
  //AOD additions since there are no AliESDtrackCuts -----------------
  //
  if (isAOD){
    // TPC #clusteres cut
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    //TODO: DCA cuts to be investigated!!!
//       pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
//       pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }
  diele->GetTrackFilter().AddCuts(pt);

  // PID cuts --------------------------------------------------------
  
  AliDielectronPID *pid = new AliDielectronPID("PID10","TPC nSigma |e|<3 + |Pi|>3 + |P|>3 + TOF nSigma |e|<3");
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-2.5,3.);
  pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.);

  //  pid->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-1.8.,1.8.,0.,0.,kTRUE);
  //  pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-2.,2.0.,0.,0.,kTRUE);

  //pid->AddCut(AliDielectronPID::kTOF,AliPID::kKaon,-4.,4.,0.,0.,kTRUE);
  //pid->AddCut(AliDielectronPID::kTOF,AliPID::kProton,-4.,4.,0.,0.,kTRUE);
  
  //pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.,3.,0.,0.,kTRUE);
  //  pid->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-1.75.,1.8.,0.,0.,kTRUE);
  //pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-1.3.,2.0.,0.,0.,kTRUE);



  //pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.,3.);
  //pid->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3.,3.);
  //pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.,3.);

  //pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.);

  //pid->AddCut(AliDielectronPID::kTOF,AliPID::kPion,-3.,3.);
  //pid->AddCut(AliDielectronPID::kTOF,AliPID::kKaon,-3.,3.);
  //pid->AddCut(AliDielectronPID::kTOF,AliPID::kProton,-3.,3.);

  diele->GetTrackFilter().AddCuts(pid);

}

//______________________________________________________________________________________
void SetupPairCutsDieleData(AliDielectronTaku *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setup the pair cuts
  //
  
  //Invariant mass and rapidity selection
  AliDielectronVarCuts *pairCut=new AliDielectronVarCuts("2<M<4+|Y|<.9","2<M<4 + |Y|<.9");
  //   pairCut->AddCut(AliDielectronVarManager::kM,2.,4.);
  pairCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  //pairCut->AddCut(AliDielectronVarManager::kOpeningAngle, 0, 0.1);
  //diele->GetPairFilter().AddCuts(pairCut);
  diele->GetPairPreFilter().AddCuts(pairCut);
}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition)
{
  cout<<" SetupESDtrackCutsDieleData "<<endl;
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

  //esdTrackCuts->SetPtRange(.8,1e30);
  esdTrackCuts->SetPtRange(0.2,10);
  esdTrackCuts->SetPtRange(0.2,10);

  //esdTrackCuts->SetMinNClustersTPC(120);
  esdTrackCuts->SetMinNCrossedRowsTPC(120);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);

  // default SPD any
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);

  if (cutDefinition==0)
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

  return esdTrackCuts;
}


//______________________________________________________________________________________
void InitHistogramsDieleData(AliDielectronTaku *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistosTaku *histos=new AliDielectronHistosTaku(diele->GetName(),diele->GetTitle());
  
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
//   histos->AddClass(Form("Pair_%s",PairClassName(AliDielectron::kEv1PMRot)));
//   histos->AddClass(Form("Track_Legs_%s",PairClassName(AliDielectron::kEv1PMRot)));
  
  
  //add histograms to event class
  if (cutDefinition==0) {
    histos->AddClass("Event");
    histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",300,-15.,15.,AliDielectronVarManager::kZvPrim);
  }
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",400,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                        100,-1,1,144,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        200,0.2,20.,100,0.,200.,AliDielectronVarManager::kPIn,
			AliDielectronVarManager::kTPCsignal,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);


  histos->UserHistogram("Track","TOF_Beta_P","TOF beta;P [GeV];TOF beta;#tracks",
                        200,0.2,20.,100,0.,1.2.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta,kTRUE);


  histos->UserHistogram("Track","TOFnSigmaPi_P","TOF nSigmaPi;P [GeV];TOF beta;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaK_P","TOF nSigmaK;P [GeV];TOF beta;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaPr_P","TOF nSigmaPr;P [GeV];TOF beta;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro,kTRUE);


      
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        201,-.01,4.01,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  
  diele->SetHistogramManager(histos);
}


void InitCFDieleData(AliDielectronTaku *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(diele->GetName(),diele->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 1.0, 2.0, 3.0, 5., 7.0, 10.0, 100.0");
  
  cf->AddVariable(AliDielectronVarManager::kY,"-5,-1,-0.9,-0.8,-0.7,-0.5,-0.3,0.3,0.5,0.7,0.8,0.9,1.0,5");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,"-0.5,0.5,1.5,2.5");
  cf->AddVariable(AliDielectronVarManager::kThetaHE, "-2.0, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 2.0");
  cf->AddVariable(AliDielectronVarManager::kThetaCS, "-2.0, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 2.0");
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 70, 75, 80, 85, 90, 100, 120, 160",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-5,-1,-0.9,-0.85,-0.8,-0.75,0.75,0.8,0.85,0.9,1.0,5",kTRUE);

//   cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-2.5,-2,-1.5,-1,-0.5,4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.,3.5,4.,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.,3.5,4.,100",kTRUE);

  //event variables
  cf->AddVariable(AliDielectronVarManager::kNaccTrcklts,"0.0, 9.0, 17.0, 25.0, 36.0, 55.0, 500.0");

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

//______________________________________________________________________________________
void InitTreesDieleData(AliDielectronTaku *diele, Int_t cutDefinition, Bool_t isAOD)
{

  //Setup histogram Manager
  AliDielectronHistosTaku *histos=new AliDielectronHistosTaku(diele->GetName(),diele->GetTitle());
  histos->UserTree("tree","sinlge tree");
  
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPx);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPy);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPz);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPt);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kP);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kXv);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kYv);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kZv);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kOneOverPt);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPhi);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTheta);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kEta);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kY);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kE);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kM);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kCharge);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNclsITS);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNclsTPC);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNclsTPCiter1);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNFclsTPC);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNFclsTPCr);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNFclsTPCrFrac);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTPCsignalN);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTPCsignalNfrac);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTPCchi2Cl);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTrackStatus);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNclsTRD);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTRDntracklets);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTRDpidQuality);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTRDprobEle);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTRDprobPio);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kImpactParXY);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kImpactParZ);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTrackLength);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPdgCode);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPdgCodeMother);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPdgCodeGrandMother);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNumberOfDaughters);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kHaveSameMother);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kIsJpsiPrimary);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSsignal);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSsignalSSD1);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSsignalSSD2);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSsignalSDD1);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSsignalSDD2);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSclusterMap);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSnSigmaEle);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSnSigmaPio);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSnSigmaMuo);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSnSigmaKao);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kITSnSigmaPro);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPIn);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTPCsignal);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTOFsignal);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTOFbeta);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTPCnSigmaEle);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTPCnSigmaPio);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTPCnSigmaMuo);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTPCnSigmaKao);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTPCnSigmaPro);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTOFnSigmaEle);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTOFnSigmaPio);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTOFnSigmaMuo);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTOFnSigmaKao);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTOFnSigmaPro);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kKinkIndex0);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kChi2NDF);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kDecayLength);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kR);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kOpeningAngle);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kThetaHE);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPhiHE);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kThetaCS);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPhiCS);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kLegDist);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kLegDistXY);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kDeltaEta);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kDeltaPhi);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kMerr);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kDCA);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPairType);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kPseudoProperTime);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kXvPrim);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kYvPrim);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kZvPrim);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kXRes);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kYRes);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kZRes);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNTrk);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kTracks);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNacc);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNaccTrcklts);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNch);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kCentrality);
  histos->SetReserveVariableInTree(AliDielectronVarManager::kNevents);

  diele->SetTreeManager(histos);
  
}
