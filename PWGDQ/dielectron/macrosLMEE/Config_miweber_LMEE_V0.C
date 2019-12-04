void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void SetupCuts(AliDielectron *die, Int_t cutDefinition);
AliDielectronPID *SetPIDcuts(Int_t cutDefinition);
AliDielectronPID *SetPreFilterPIDcuts(Int_t cutDefinition);
const AliDielectronEventCuts *GetEventCuts();

Bool_t isRandomRejTask=kFALSE;//needed for InitHistograms() //dont change!!!
Bool_t kRot = kFALSE;
Bool_t kMix = kFALSE;
Bool_t kNoPairing   = kFALSE;
Bool_t randomizeDau = kTRUE;
     
TString names("V0_noPID;V0_TPC;V0_TPC_TOF;V0_TPC_TOFpid;V0_TPC_ITS;V0_TPC_ITSpid;V0_ITSpid_TPCpid_TOF;V0_ITS_TPCpid_TOFpid;V0_ITSpid_TPCpid_TOFpid");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();

AliDielectron* Config_miweber_LMEE_V0(Int_t cutDefinition=1, Bool_t bESDANA = kFALSE, Bool_t bCutQA = kFALSE, Bool_t bCutV0 = kTRUE, Bool_t bCutV0Tight = kFALSE, Bool_t isRandomRej=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  // (inspired from simpler task of Theo Broeker w/o the usage of a CutLib file)
  //

  isRandomRejTask=isRandomRej;

  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  
  AliDielectron *die =
    new AliDielectron(Form("%s",name.Data()),
                      Form("Track cuts: %s",name.Data()));

  if(bCutQA){
    die->SetCutQA(bCutQA);
  }
  
  if(kRot){
    AliDielectronTrackRotator *rot = new AliDielectronTrackRotator;
    rot->SetConeAnglePhi(TMath::Pi());
    rot->SetIterations(10);
    die->SetTrackRotator(rot);
  }//kRot
  
  
  if(kMix && !(die->GetHasMC()) ){ // need second since there is a problem when mixing MC events (TRef?)
    AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;

    mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
    mix->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,30,50,80");
    mix->SetDepth(15);
    mix->SetMixType(AliDielectronMixingHandler::kAll);
    
    // using TPC event plane, uncorrected. (also, the old phi range was wrong, now same effective binning.)
    // mix->AddVariable(AliDielectronVarManager::kTPCrpH2uc, 6, TMath::Pi()/-2., TMath::Pi()/2.);
    
    die->SetMixingHandler(mix);
  }//kMix


  // set track cuts
  SetupCuts(die,cutDefinition,bESDANA,bCutV0,bCutV0Tight);

 //
 // histogram setup
 // only if an AliDielectronHistos object is attached to the
 // dielectron framework histograms will be filled
 //
 
 InitHistograms(die,cutDefinition);
 //  InitCF(die,cutDefinition);
 
 die->SetNoPairing(kNoPairing);
 
 return die;

}

//______________________________________________________________________________________
void SetupCuts(AliDielectron *die, Int_t cutDefinition, Bool_t bESDANA = kFALSE, Bool_t bCutV0 = kTRUE, Bool_t bCutV0Tight = kFALSE)
{
  // Setup the track cuts

  //pairing with TLorentzVector
  die->SetUseKF(kFALSE);
      
  // V0 cuts
  if(bCutV0)
    die->GetTrackFilter().AddCuts(SetupV0Cuts(bCutV0Tight));

  // track cuts
  die->GetTrackFilter().AddCuts(SetupTrackCuts());
  
  // PID cuts
  if(cutDefinition > 0)
    die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));
  
}


//______________________________________________________________________________________
AliDielectronV0Cuts *SetupV0Cuts(Bool_t bCutV0Tight = kFALSE){

  AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");

  gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kAll);  // kAll by default
  gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
    
  if(!bCutV0Tight){
    gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);// ( 0.02 -- 0.05 )
    gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE); 
    gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.2, kFALSE);// ( 0.05 -- 0.2 )
    gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.1, kFALSE);// ( 0.05 -- 0.1 )
    gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  }
  else{
    gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.05),   1.0,  kFALSE); // ( 0.02 -- 0.05 )
    gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE); // ( 0.05 -- 0.2 )
    gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE); // ( 0.05 -- 0.1 )
    gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  }

  gammaV0Cuts->SetExcludeTracks(kFALSE);

  return gammaV0Cuts;
}

//______________________________________________________________________________________
AliESDtrackCuts *SetupTrackCuts(){

  AliESDtrackCuts* trackCuts=new AliESDtrackCuts();
  trackCuts->SetTitle("track cuts");
  
  // pT and eta
  trackCuts->SetPtRange(0.05, 1e30);
  trackCuts->SetEtaRange(-0.8, 0.8);
  //TPC
  trackCuts->SetRequireTPCRefit(kTRUE);
  trackCuts->SetMinNCrossedRowsTPC(100);
  trackCuts->SetMinNClustersTPC(80);
  trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  trackCuts->SetMaxChi2PerClusterTPC(4.0);
  trackCuts->SetMaxFractionSharedTPCClusters(0.4);
  
  return trackCuts;
}

//______________________________________________________________________________________
AliDielectronPID *SetPIDcuts(Int_t cutDefinition){
  
  AliDielectronPID *pid = new AliDielectronPID();
  pid->SetTitle("PID cuts");

  if(cutDefinition == 1){
    
    // TPC PID
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);

  }

  // TPC PID + TOF
  else if(cutDefinition == 2){
    
    //TPC
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);

    // TOF matching required
    pid->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -999. ,  1e30,   0.0, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
  }
  
  // TPC PID + TOF PID
  else if(cutDefinition == 3){
    
    //TPC
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);

    // TOF required
    pid->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.0, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
  }

  // TPC PID + ITS
  else if(cutDefinition == 4){
    
    //TPC
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);

    // ITS matching required
    pid->AddCut(AliDielectronPID::kITS, AliPID::kElectron,     -999. ,  1e30,   0.0, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
  }

   // TPC PID + ITS PID
  else if(cutDefinition == 5){
    
    //TPC
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);

    // ITS required
    pid->AddCut(AliDielectronPID::kITS, AliPID::kElectron,     -3. ,  3.,   0.0, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
  }

  // ITS PID + TPC PID + TOF
  else if(cutDefinition == 6){
    
    //TPC
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);

    // ITS required
    pid->AddCut(AliDielectronPID::kITS, AliPID::kElectron,     -3. ,  3.,   0.0, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);

    // TOF matching required
    pid->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -999. ,  1e30,   0.0, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
  }

    // ITS  + TPC PID + TOF PID
  else if(cutDefinition == 7){
    
    //TPC
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);

    // ITS matching required
    pid->AddCut(AliDielectronPID::kITS, AliPID::kElectron,     -999. ,  1e30,   0.0, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);

     // TOF required
    pid->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.0, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
  }

  // ITS PID  + TPC PID + TOF PID
  else if(cutDefinition == 8){
    
    //TPC
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);

    // ITS required
    pid->AddCut(AliDielectronPID::kITS, AliPID::kElectron,     -3. ,  3.,   0.0, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
    
    // TOF required
    pid->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.0, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
  }
  
  return pid;
}


//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  
  //Setup histogram classes
  AliDielectronHistos *histos=
    new AliDielectronHistos(die->GetName(),
                           die->GetTitle());
  
  //Initialise histogram classes
  //histos->SetReservedWords("Track;Pair");
  histos->SetReservedWords("Track;Pair;Track_Legs;Pre;RejPair;RejTrack;Random");
  //histos->SetReservedWords("Track");  

  //Event class
  histos->AddClass("Event");

  if(!isRandomRejTask){
    //Track classes
    //to fill also track info from 2nd event loop until 2
    for (Int_t i=0; i<2; ++i){
      histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
    }
    //Pair classes
    for (Int_t i=0; i<3; ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
    }
    //ME and track rot
    if (die->GetMixingHandler()) {
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
    }
    if (die->GetTrackRotator()) {
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));
    }
  }
  if(isRandomRejTask){
    //
    // _____ histograms for AliAnalysisTaskMultiDielectronPR _____
    //
    //    histos->AddClass("Rand_Pair");
    //    histos->AddClass("Rand_RejPair");
    const char* cRandomPairClassNames[2] = { "Testpart", "RejTestpart" };
    for (Int_t i=0; i<2; ++i){
      histos->AddClass(Form("Random_%s",cRandomPairClassNames[i]));
    }
    histos->UserHistogram("Random","Pt","",200,0,10.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Random","Eta","",200,-2,2,AliDielectronVarManager::kEta);
    histos->UserHistogram("Random","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    histos->UserHistogram("Random","Px","",200,0,10.,AliDielectronVarManager::kPx);
    histos->UserHistogram("Random","Py","",200,0,10.,AliDielectronVarManager::kPy);
    histos->UserHistogram("Random","Pz","",200,0,10.,AliDielectronVarManager::kPz);
    histos->UserHistogram("Random","Pt_Eta_Phi","",
                          500,0.,10.,16,-0.8,0.8,30,0.,2*TMath::Pi(),
                          AliDielectronVarManager::kPt,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  }
  
  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0,1,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","ZVertex","ZVertex;ZVertex/cm",120,-12.,12.,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","Centrality","Centrality;Centrality/%",202,-1.,100.,AliDielectronVarManager::kCentralityNew);

  //add histograms to track class
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","P","P;P [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kP);
  histos->UserHistogram("Track","PIn","PIn;PIn [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track","Eta_phi","Eta vs Phi;Eta;Phi",90,-0.9,0.9,160,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Pt_phi","Pt vs Phi;Pt;Phi [GeV];#tracks",500,0.,5.,320,0.,6.4,AliDielectronVarManager::kPt,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","ImpParXY","ImpParXY; ImpParXY ;#tracks",100,-5.,5.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","ImpParZ","ImpParZ; ImpParZ ;#tracks",100,-5.,5.,AliDielectronVarManager::kImpactParZ);
  
  histos->UserHistogram("Track","NClusterTPC","NClusterTPC; NClusterTPC ;#tracks",200,-0.5,199.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","NClusterTPCShared","NClusterTPCShared; NClusterTPCShared ;#tracks",200,-0.5,199.5,AliDielectronVarManager::kNclsSTPC);
  histos->UserHistogram("Track","CrossedRows","CrossedRows; CrossedRows ;#tracks",200,-0.5,199.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","CrossedRowsOverFindable","CrRowsOverFindable; CrRows/FindableCls ;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCchi2perCls","TPCchi2perCls; TPCchi2perCls ;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  
  histos->UserHistogram("Track","NClusterITS","NClusterITS; NClusterITS ;#tracks",8,-0.5,7.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","ITSchi2perCls","ITSchi2perCls; ITSchi2perCls ;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);

  histos->UserHistogram("Track","TPCdEdx_P","dEdx;P [GeV];TPC signal (arb units) vs Momentum;Mom;TPCsignal",     200,0.,10.,150,  0.,150. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);  
  histos->UserHistogram("Track","TPCnSigma_MomEle","TPC number of sigmas Electrons vs Momentum;Mom;TPCsigmaEle", 200,0.,10.,300,-30., 30. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","ITSdEdx_P","dEdx;P [GeV];ITS signal (arb units) vs Momentum;Mom;ITSsignal",     200,0.,10.,150,  0.,150. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","ITSnSigma_MomEle","ITS number of sigmas Electrons vs Momentum;Mom;ITSsigmaEle"                           ,     200,0.,10.,300,-30., 30. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","TOFbeta_Mom","kTOFbeta vs Momentum;Mom;TOFbeta"                           ,     200,0.,10.,120,  0.,  1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFnSigma_MomEle","TOF number of sigmas Electrons vs Momentum;Mom;TOFsigmaEle"                           ,     200,0.,10.,300,-30., 30. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
  histos->UserHistogram("Track","TOFnSigma_MomEle_large","TOF number of sigmas Electrons vs Momentum;Mom;TOFsigmaEle"                           ,     200,0.,10.,300,-1500., 1500. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);

 
  //add histograms to pair classes
  histos->UserHistogram("Pair",
                        "InvMass_pPt","Inv.Mass:PairPt;Inv. Mass (GeV/c^{2});Pair Pt (GeV/c)",
                        500,0.,0.05,250,0.,5.,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt);

  histos->UserHistogram("Pair","InvMass_R","Inv.Mass:R;Inv. Mass (GeV/c^{2});R",
                        500,0.,0.05,500,0.,50.,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kR);  
  						
   histos->UserHistogram("Pair","ArmAlpha_ArmPt","ArmAlpha:ArmPt;#alpha;Armenteros Pt (GeV/c)",
                        500,-1.,1.,500,0.,0.5,
                        AliDielectronVarManager::kArmAlpha, AliDielectronVarManager::kArmPt); 

  die->SetHistogramManager(histos);

}

const AliDielectronEventCuts *GetEventCuts(){

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex SPD && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny); // AOD
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1); 

  return eventCuts;
}



