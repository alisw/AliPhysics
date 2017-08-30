void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void SetupCuts(AliDielectron *die, Int_t cutDefinition);
AliDielectronPID *SetPIDcuts(Int_t cutDefinition);
AliDielectronPID *SetPreFilterPIDcuts(Int_t cutDefinition);
const AliDielectronEventCuts *GetEventCuts();

Bool_t isRandomRejTask=kFALSE;//needed for InitHistograms() //dont change!!!
Bool_t kRot = kFALSE;
Bool_t kMix = kTRUE;
Bool_t kNoPairing   = kFALSE;
Bool_t randomizeDau = kTRUE;
     
TString names("PbPbData0;PbPbData1;PbPbData2;PbPbData3;PbPbData4");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();
const Int_t nPF = 2; // use prefiltering for cuts lower than nPF

AliDielectron* Config_miweber_LMEE_PbPb_woCutLib(Int_t cutDefinition=1, Bool_t bESDANA = kFALSE, Bool_t bCutQA = kFALSE, Bool_t isRandomRej=kFALSE)
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
  else if(cutDefinition==666)
    name="PbPbData666";//special cuts nano AOD filtering
  else if(cutDefinition==667)
    name="PbPbData667";//special cuts nano AOD filtering
  
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
  SetupCuts(die,cutDefinition,bESDANA);

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
void SetupCuts(AliDielectron *die, Int_t cutDefinition, Bool_t bESDANA = kFALSE)
{
  // Setup the track cuts

  //pairing with TLorentzVector
  die->SetUseKF(kFALSE);
      
  AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","conversion tagging");
  noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
     
  AliDielectronVarCuts* pairCutsInvM  = new AliDielectronVarCuts("pairCutsInvM","pairCutsInvM");
  AliDielectronVarCuts* pairCutsOpAng = new AliDielectronVarCuts("pairCutsOpAng","pairCutsOpAng");  
        
  if(cutDefinition < nPF){

    Printf("Use prefiltering!");

    if(bESDANA){
      die->GetTrackFilter().AddCuts(SetupPreFilterESDtrackCuts(cutDefinition));
    }
    else{
      die->GetTrackFilter().AddCuts(SetupPreFilterAODtrackCuts(cutDefinition));
    }
    die->GetTrackFilter().AddCuts(SetPreFilterPIDcuts(cutDefinition));


    //pairPrefilter
    AliAnalysisCuts* pairPreCuts=0x0;

    pairCutsInvM ->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
    pairCutsOpAng->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.050); 

    AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
    pairCutsCG->AddCut(pairCutsInvM);
    pairCutsCG->AddCut(pairCutsOpAng);
    pairPreCuts = pairCutsCG;
    
    die->GetPairPreFilter().AddCuts(pairCutsOpAng);
    
    //FinalTrackCuts after prefiltering
    die->GetPairPreFilterLegs().AddCuts(SetPIDcuts(cutDefinition));
    if(bESDANA){
      die->GetPairPreFilterLegs().AddCuts(SetupESDtrackCuts(cutDefinition));
    }
    else{
      die->GetPairPreFilterLegs().AddCuts(SetupAODtrackCuts(cutDefinition));
    }
    // die->GetPairPreFilterLegs().AddCuts(noconv);
    
    die->SetPreFilterUnlikeOnly(kTRUE);
  }
  else if( cutDefinition == 666 ){
    // loose cuts as used for nano AOD filtering (only for AODs)
    if(!bESDANA){      
      SetupAODtrackCutsLoose(die);
    }
  }
  else if( cutDefinition == 667 ){
    // tighter cuts than used for nano AOD filtering (only for AODs)
    if(!bESDANA){      
      SetupAODtrackCutsTight(die);
    }
  }
  else{
    if(bESDANA){
      die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));
    }
    else{
      die->GetTrackFilter().AddCuts(SetupAODtrackCuts(cutDefinition));
    }    
    die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));
    // die->GetTrackFilter().AddCuts(noconv);
  }

  // removing strange mixing pairs
  AliDielectronVarCuts *opAngleCut = new AliDielectronVarCuts("opAngleCut","opAngleCut");
  opAngleCut->AddCut(AliDielectronVarManager::kOpeningAngle,-1.e-10,1.e-10,kTRUE);   
  //die->GetPairFilter().AddCuts(opAngleCut); 
          
}
//______________________________________________________________________________________
//-----------------------------------pid------------------------------------------------

AliDielectronPID *SetPIDcuts(Int_t cutDefinition){
  
  AliDielectronPID *pid = new AliDielectronPID();
  pid->SetTitle("PID cuts");

  if(cutDefinition == 0){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,5. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -2., 2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition == 1){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,5. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -2., 2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition == 2){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,5. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -2., 2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition == 3){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,5. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.8,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition == 4){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  
 return pid;
}


//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition){

  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts;
  fesdTrackCuts->SetTitle("ESD track cuts");

  // //global
  // fesdTrackCuts->SetPtRange( 0.2 , 100. );
  // fesdTrackCuts->SetEtaRange( -0.8 , 0.8 );
  // fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  // fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  // fesdTrackCuts->SetDCAToVertex2D(kFALSE);
  // fesdTrackCuts->SetMaxDCAToVertexZ(3.);
  // fesdTrackCuts->SetMaxDCAToVertexXY(1.);
  
  // fesdTrackCuts->SetRequireITSRefit(kTRUE);
  // fesdTrackCuts->SetRequireTPCRefit(kTRUE);
 
  // if(cutDefinition == 0){
  //   fesdTrackCuts->SetMinNClustersITS(4);
  //   fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
  //   fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  //   fesdTrackCuts->SetMinNClustersTPC(80);
  //   fesdTrackCuts->SetMinNCrossedRowsTPC(100);
  //   fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
  //   fesdTrackCuts->SetMaxChi2PerClusterTPC(4);  
  // }
  // if(cutDefinition == 1){
  //   fesdTrackCuts->SetMinNClustersITS(5);
  //   fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
  //   fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  //   fesdTrackCuts->SetMinNClustersTPC(100);
  //   fesdTrackCuts->SetMinNCrossedRowsTPC(130);
  //   fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
  //   fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  // }
  
  return fesdTrackCuts;
}

AliDielectronPID *SetPreFilterPIDcuts(Int_t cutDefinition){
  AliDielectronPID *pid = new AliDielectronPID();

  if(cutDefinition==0){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -4., 4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  else if(cutDefinition==1){
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-5.,5.,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -4., 4. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable    ,AliDielectronVarManager::kPt);
  }   
  return pid;
}

AliESDtrackCuts *SetupPreFilterESDtrackCuts(Int_t cutDefinition){

  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts;
  fesdTrackCuts->SetTitle("ESD track cuts (pre)");

  // //global
  // fesdTrackCuts->SetPtRange( 0.08 , 100. );
  // fesdTrackCuts->SetEtaRange( -1.1 , 1.1 );
  // fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  // fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  // fesdTrackCuts->SetDCAToVertex2D(kFALSE);
  // fesdTrackCuts->SetMaxDCAToVertexZ(3.);
  // fesdTrackCuts->SetMaxDCAToVertexXY(1.);
 
  // //ITS cuts
  // fesdTrackCuts->SetRequireITSRefit(kTRUE);
  // fesdTrackCuts->SetMinNClustersITS(3);
  // fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);

  return fesdTrackCuts;
  
}

//______________________________________________________________________________________
AliAnalysisCuts *SetupAODtrackCuts(Int_t cutDefinition){

  // so far only one cut definition 

  AliAnalysisCuts* faodTrackCuts=0x0;

  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 2 with PID
  trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
  trackCutsAOD->AddCut(AliDielectronVarManager::kPt, 0.2,100.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -0.8,0.8);


  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2010(kFALSE); loose DCA, 2D cut
  trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  
  AliDielectronCutGroup* cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
  cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
  cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);

  faodTrackCuts = cgTrackCutsAnaSPDfirst;
  faodTrackCuts->SetTitle("AOD track cuts");

  return faodTrackCuts;
}


AliAnalysisCuts *SetupPreFilterAODtrackCuts(Int_t cutDefinition){

  AliAnalysisCuts* faodTrackCuts=0x0;

  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kPt, 0.08,100.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -1.1,1.1);

  AliDielectronTrackCuts* trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");

  if(cutDefinition==1){
    trackCutsDiel->SetAODFilterBit(1<<0|1<<1); // ITSsa cuts + TPC only cuts
    trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }
  else{
    trackCutsDiel->SetAODFilterBit(1<<0); // TPC only cuts 
    trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }
  AliDielectronCutGroup* cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
  cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
  cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);

  faodTrackCuts = cgTrackCutsAnaSPDfirst;
  faodTrackCuts->SetTitle("AOD track cuts (pre)");
  
  return faodTrackCuts;
}

//______________________________________________________________________________________
void SetupAODtrackCutsLoose(AliDielectron *die)
{
  //
  // Setup the track cuts
  // - these are the cuts used for nano AOD creation LHC15o (ConfigLMEE_nano_PbPb2015.C)
  //

  Bool_t hasMC=die->GetHasMC();

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany); // I think we loose the possibility to use prefilter?

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  // specific cuts
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);

  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 8.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  // TOF
  // pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-4.,4.,0.,0.,kFALSE, AliDielectronPID::kIfAvailable);
  // TPC
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
  // ITS
  pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-4.,4.);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  pidCutsMC->SetCutOnMCtruth(kTRUE);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11., -11.);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +11., +11.);

  // activate the cut sets (order might be CPU timewise important)
  if(hasMC) {
    // cuts->AddCut(pidCutsMC); // commnted out for the moment
    // die->GetTrackFilter().AddCuts(pidCutsMC);
    cuts->AddCut(trkFilter);
    die->GetTrackFilter().AddCuts(trkFilter);
    cuts->AddCut(varCuts);
    die->GetTrackFilter().AddCuts(varCuts);
    cuts->AddCut(trkCuts);
    die->GetTrackFilter().AddCuts(trkCuts);
    cuts->AddCut(pidCuts);
    die->GetTrackFilter().AddCuts(pidCuts);
  }
  else {
    cuts->AddCut(trkFilter);
    die->GetTrackFilter().AddCuts(trkFilter);
    cuts->AddCut(varCuts);
    die->GetTrackFilter().AddCuts(varCuts);
    cuts->AddCut(trkCuts);
    die->GetTrackFilter().AddCuts(trkCuts);
    cuts->AddCut(pidCuts);
    die->GetTrackFilter().AddCuts(pidCuts);
  }
  cuts->Print();
  //die->GetTrackFilter().AddCuts(cuts); // done now separately

}

//______________________________________________________________________________________
void SetupAODtrackCutsTight(AliDielectron *die)
{
  //
  // Setup the track cuts
  // - these are similar cuts used for nano AOD creation LHC15o (ConfigLMEE_nano_PbPb2015.C), but tighter
  //

  Bool_t hasMC=die->GetHasMC();

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany); // I think we loose the possibility to use prefilter?

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  // specific cuts
  trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
  varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0, 6.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster
  varCuts->AddCut(AliDielectronVarManager::kNclsSTPC,     0.0, 1.0);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   6.0);
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);

  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 8.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  // TPC
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,5,0.,0.,kTRUE);
  // ITS
  pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-3.,3.);
  // TOF
  pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.);
  // pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-4.,4.,0.,0.,kFALSE, AliDielectronPID::kIfAvailable);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  pidCutsMC->SetCutOnMCtruth(kTRUE);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11., -11.);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +11., +11.);

  // activate the cut sets (order might be CPU timewise important)
  if(hasMC) {
    // cuts->AddCut(pidCutsMC); // commnted out for the moment
    // die->GetTrackFilter().AddCuts(pidCutsMC);
    cuts->AddCut(trkFilter);
    die->GetTrackFilter().AddCuts(trkFilter);
    cuts->AddCut(varCuts);
    die->GetTrackFilter().AddCuts(varCuts);
    cuts->AddCut(trkCuts);
    die->GetTrackFilter().AddCuts(trkCuts);
    cuts->AddCut(pidCuts);
    die->GetTrackFilter().AddCuts(pidCuts);
  }
  else {
    cuts->AddCut(trkFilter);
    die->GetTrackFilter().AddCuts(trkFilter);
    cuts->AddCut(varCuts);
    die->GetTrackFilter().AddCuts(varCuts);
    cuts->AddCut(trkCuts);
    die->GetTrackFilter().AddCuts(trkCuts);
    cuts->AddCut(pidCuts);
    die->GetTrackFilter().AddCuts(pidCuts);
  }
  cuts->Print();
  //die->GetTrackFilter().AddCuts(cuts); // done now separately

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
      // Legs of final Pairs. Both charges together. No duplicate entries.
      if(cutDefinition < nPF) 
	histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
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
    if(cutDefinition < nPF){  
      //PreFilter Classes
      //to fill also track info from 2nd event loop until 2
      for (Int_t i=0; i<2; ++i){
        histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
      }
      
      //Create Classes for Rejected Tracks/Pairs:
      for (Int_t i=0; i<3; ++i){
        histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
        // Legs of rejected Pairs. Both charges together. One track can and will make multiple entries.
        histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistogramsPair(...)'
      }
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

  //add histograms to pair classes
  histos->UserHistogram("Pair",
                        "InvMass_pPt","Inv.Mass:PairPt;Inv. Mass (GeV/c^{2});Pair Pt (GeV/c)",
                        500,0.,5.,250,0.,5.,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt);

  histos->UserHistogram("Pair","InvMass_PairPt_PhivPair","InvMass:PairPt:PhivPair;Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                        600,0.,6., 600,0.,6., 20,0.,TMath::Pi(),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);  
  						
  histos->UserHistogram("Pair",
                        "Eta_phi_pair","Eta vs Phi (pair);Eta;Phi",
                        50,-1.,1.,80,0.,6.4,
                        AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
			
  histos->UserHistogram("Pair",
                        "InvMass_PhivPair","InvMass_PhivPair;InvMass;PhivPair",
                         50, 0. , 0.5, 160 , 0., 3.2,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPhivPair );

  histos->UserHistogram("Pair",
		            	"InvMass_OpAngle","InvMass_OpAngle;Invariant Mass;Opening angle",
		            	100, 0., 0.5 ,160, 0., 3.2,
		            	AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);

  histos->UserHistogram("Pair",
                        "Y","Y;counts;Y",
                        60, -1.2 , 1.2, 
                        AliDielectronVarManager::kY);

  if(cutDefinition < nPF){ 
    histos->UserHistogram("Pre","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);

    histos->UserHistogram("RejPair","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);

    histos->UserHistogram("RejPair",
                          "OpAngle_InvMass","InvMass_openingAngle;Invariant Mass;opening angle",
                          100, 0., 0.2, 100, 0. ,0.2,
                          AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
                        
    histos->UserHistogram("RejTrack","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt); 
  
    histos->UserHistogram("Track_Legs","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt); 
  }
  die->SetHistogramManager(histos);

}

const AliDielectronEventCuts *GetEventCuts(){

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex SPD && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny); // AOD
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1); 

  //no centrality cuts for the moment
  //Bool_t isRun2 = kTRUE;
  //eventCuts->SetCentralityRange(0,80,isRun2);

  return eventCuts;
}

