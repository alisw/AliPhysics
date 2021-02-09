void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void SetTPCCorr(AliDielectron *die);
void SetupMCsignals(AliDielectron* die);
void SetupCuts(AliDielectron *die, Int_t cutDefinition, Bool_t bESDANA);
AliDielectronEventCuts *GetEventCuts(Double_t centMin, Double_t centMax, Bool_t reqAliEventCuts, Bool_t reqAliEventCutsCorrelated);
AliDielectronPID *SetPIDcuts(Int_t cutDefinition);
AliDielectronPID *SetPreFilterPIDcuts(Int_t cutDefinition);
AliESDtrackCuts *SetupPreFilterESDtrackCuts(Int_t cutDefinition);
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);
AliAnalysisCuts *SetupPreFilterAODtrackCuts(Int_t cutDefinition);
AliAnalysisCuts *SetupAODtrackCuts(Int_t cutDefinition);
void SetupAODtrackCutsLoose(AliDielectron *die);
void SetupAODtrackCutsTight(AliDielectron *die);
void SetupAODtrackCutsTight2(AliDielectron *die);
void SetupAODtrackCutsTight3(AliDielectron *die);
void SetupAODtrackCutsCarsten1(AliDielectron *die);
void SetupAODtrackCutsCarsten2(AliDielectron *die);
void SetupAODtrackCutsRun3(AliDielectron *die);
void SetupAODtrackCutsRun3_FT2(AliDielectron *die);
void SetupAODtrackCutsRun3_FT2_noPID(AliDielectron *die);
void SetupAODtrackCutsRun3_FT2_ESD(AliDielectron *die);
void SetupAODtrackCutsRun3_FT2_noPID_ESD(AliDielectron *die);


Bool_t isRandomRejTask=kFALSE;//needed for InitHistograms() //dont change!!!

// available cut defintions
const Int_t nPF       = 2; // use prefiltering for cuts < nPF
const Int_t nNoPF     = 5; // use no prefiltering for nPF <= cuts < nNoPF 
const Int_t nExtraMin = 666; // use extra cuts for nExtraMin <= cuts < nExtraMax 
const Int_t nExtraMax = 677; // use extra cuts for nExtraMin <= cuts < nExtraMax 

AliDielectron* Config_miweber_LMEE_PbPb_woCutLib(Int_t cutDefinition=1,
						 Bool_t bESDANA = kFALSE,
						 Bool_t bCutQA = kFALSE,
						 Bool_t isRandomRej=kFALSE,
						 Bool_t useTPCCorr=kFALSE,
						 Bool_t useRotation=kFALSE,
						 Bool_t useMixing=kTRUE,
						 Bool_t noPairing=kFALSE,
						 Bool_t hasMC=kFALSE
						 )
{
  //
  // Setup the instance of AliDielectron
  // (inspired from simpler task of Theo Broeker w/o the usage of a CutLib file)
  //

  isRandomRejTask=isRandomRej;

  // First check if cut definition is defined
  if(!((cutDefinition >= 0 && cutDefinition < nNoPF) || (cutDefinition >= nExtraMin && cutDefinition < nExtraMax))){
    Printf("=====================================");
    Printf("Cut set = %d is not available.",cutDefinition);
    Printf("Use [0,%d[ for cuts with prefilter",nPF);
    Printf("Use [%d,%d[ for cuts w/o prefilter",nPF,nNoPF);
    Printf("Use [%d,%d[ for extra cuts",nExtraMin,nExtraMax);
    Printf("=====================================");
    return NULL;
  }
  else{
    Printf("=====================================");
    Printf("Configuring cut set %d (MC = %d)",cutDefinition,hasMC);
    Printf("=====================================");
  }
    

  // create the actual framework object
  TString name=Form("PbPbData%d",cutDefinition);
  
  AliDielectron *die =
    new AliDielectron(Form("%s",name.Data()),
                      Form("Track cuts: %s",name.Data()));

  die->SetHasMC(hasMC);

  if(bCutQA){
    die->SetCutQA(bCutQA);
  }
  
  if(useRotation){
    AliDielectronTrackRotator *rot = new AliDielectronTrackRotator;
    rot->SetConeAnglePhi(TMath::Pi());
    rot->SetIterations(10);
    die->SetTrackRotator(rot);
  }//useRotation
  
  if(useMixing){
    
    AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
    
    mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
    mix->AddVariable(AliDielectronVarManager::kCentralityNew,"0,5,10,20,30,50,80");
    mix->SetDepth(15);
    mix->SetMixType(AliDielectronMixingHandler::kAll);
    
    // using TPC event plane, uncorrected. (also, the old phi range was wrong, now same effective binning.)
    if(cutDefinition==671 || cutDefinition==672){
     mix->AddVariable(AliDielectronVarManager::kQnTPCrpH2, 6, TMath::Pi()/-2., TMath::Pi()/2.);
    }
    
    die->SetMixingHandler(mix);
  }//useMixing

  // set track cuts
  SetupCuts(die,cutDefinition,bESDANA);
  
  if(useTPCCorr)  SetTPCCorr(die);
 
 //
 // histogram setup
 // only if an AliDielectronHistos object is attached to the
 // dielectron framework histograms will be filled
 //
  if (die->GetHasMC()) SetupMCsignals(die);
 
 InitHistograms(die,cutDefinition);
 //  InitCF(die,cutDefinition);
 
 die->SetNoPairing(noPairing);
 
 return die;

}

//______________________________________________________________________________________
void SetupCuts(AliDielectron *die, Int_t cutDefinition, Bool_t bESDANA)
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
  else if( cutDefinition == 668 ){
    // another set of tight cuts (only for AODs)
    if(!bESDANA){      
      SetupAODtrackCutsTight2(die);
    }
  }
  else if( cutDefinition == 669 ){
    // another set of tight cuts (only for AODs)
    if(!bESDANA){      
      SetupAODtrackCutsTight3(die);
    }
  }
  else if( cutDefinition == 670 || cutDefinition == 671 || cutDefinition == 672 ){
    // another set of tight cuts (only for AODs)
    if(!bESDANA){      
      SetupAODtrackCutsCarsten1(die);
    }
  }
  else if( cutDefinition == 673 ){
    // set for Run 3 analysis (as 669, but w/o ITS PID)
    if(!bESDANA){      
      SetupAODtrackCutsRun3(die);
    }
  }
  else if( cutDefinition == 674 ){
    // set for Run 3 analysis (as close as possible as for FT2)
    if(!bESDANA){      
      SetupAODtrackCutsRun3_FT2(die);
    }
    else{
      SetupAODtrackCutsRun3_FT2_ESD(die);
    }
  }
  else if( cutDefinition == 675 ){
    // set to compare to Carsten's analysis
    if(!bESDANA){      
      SetupAODtrackCutsCarsten2(die);
    }
  }
  else if( cutDefinition == 676 ){
    // set for Run 3 analysis (as close as possible as for FT2) - w/o exp PID cuts
    if(!bESDANA){      
      SetupAODtrackCutsRun3_FT2_noPID(die);
    }
    else{
      SetupAODtrackCutsRun3_FT2_noPID_ESD(die);
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
  trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2010(kFALSE); loose DCA, 2D cut
  trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
  
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
    trackCutsDiel->SetAODFilterBit((AliDielectronTrackCuts::EFilterBit)(AliDielectronTrackCuts::kTPCqual | AliDielectronTrackCuts::kITSonly)); // ITSsa cuts + TPC only cuts
    trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kAny);
  }
  else{
    trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual); // TPC only cuts 
    trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kAny);
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
  pidCutsMC->SetCutOnMCtruth(kTRUE);//un-comment to use only positive labels
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11.1, -10.9);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +10.9, +11.1);

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
  trkCuts->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
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
  pidCutsMC->SetCutOnMCtruth(kTRUE);//un-comment to use only positive labels
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11.1, -10.9);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +10.9, +11.1);

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
void SetupAODtrackCutsTight2(AliDielectron *die)
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
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts    = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronVarCuts *varCuts2   = new AliDielectronVarCuts("VarCuts2","VarCuts2");
  AliDielectronTrackCuts *trkCuts  = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  
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

  // additional 
  varCuts2->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.95, 1.05);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  // TPC
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
  // ITS
  pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-3.,3.);
  // TOF
  pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  pidCutsMC->SetCutOnMCtruth(kTRUE);//un-comment to use only positive labels
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11.1, -10.9);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +10.9, +11.1);

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
    cuts->AddCut(varCuts2);
    die->GetTrackFilter().AddCuts(varCuts2);
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
    cuts->AddCut(varCuts2);
    die->GetTrackFilter().AddCuts(varCuts2);
    cuts->AddCut(pidCuts);
    die->GetTrackFilter().AddCuts(pidCuts);
  }
  cuts->Print();
  //die->GetTrackFilter().AddCuts(cuts); // done now separately

}

//______________________________________________________________________________________
void SetupAODtrackCutsTight3(AliDielectron *die)
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
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts    = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronVarCuts *varCuts2   = new AliDielectronVarCuts("VarCuts2","VarCuts2");
  AliDielectronTrackCuts *trkCuts  = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  
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

  // additional 
  varCuts2->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.95, 1.05);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  // TPC
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
  // ITS
  pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-3.,3.);
  // TOF
  pidCuts->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.4, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  pidCutsMC->SetCutOnMCtruth(kTRUE);//un-comment to use only positive labels
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11.1, -10.9);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +10.9, +11.1);

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
    cuts->AddCut(varCuts2);
    die->GetTrackFilter().AddCuts(varCuts2);
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
    cuts->AddCut(varCuts2);
    die->GetTrackFilter().AddCuts(varCuts2);
    cuts->AddCut(pidCuts);
    die->GetTrackFilter().AddCuts(pidCuts);
  }
  cuts->Print();
  //die->GetTrackFilter().AddCuts(cuts); // done now separately

}

//______________________________________________________________________________________
void SetupAODtrackCutsRun3(AliDielectron *die)
{
  //
  // Setup the track cuts
  // - these are similar cuts used for nano AOD creation LHC15o (ConfigLMEE_nano_PbPb2015.C), but tighter
  // - no ITS PID
  //

  Bool_t hasMC=die->GetHasMC();

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts    = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronVarCuts *varCuts2   = new AliDielectronVarCuts("VarCuts2","VarCuts2");
  AliDielectronTrackCuts *trkCuts  = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  
  // specific cuts
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
  //varCuts->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster (cannot be used for Run3)
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 8.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

  // additional 
  varCuts2->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.95, 1.05);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  // TPC
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
  // TOF
  pidCuts->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.2, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  //pidCutsMC->SetCutOnMCtruth(kTRUE);//un-comment to use only positive labels
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11.1, -10.9);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +10.9, +11.1);

  // activate the cut sets (order might be CPU timewise important)
  if(hasMC){
    cuts->AddCut(pidCutsMC); // commented in for the moment -> only true electrons
    die->GetTrackFilter().AddCuts(pidCutsMC);
  }
  cuts->AddCut(trkFilter);
  die->GetTrackFilter().AddCuts(trkFilter);
  cuts->AddCut(varCuts);
  die->GetTrackFilter().AddCuts(varCuts);
  cuts->AddCut(trkCuts);
  die->GetTrackFilter().AddCuts(trkCuts);
  cuts->AddCut(varCuts2);
  die->GetTrackFilter().AddCuts(varCuts2);
  cuts->AddCut(pidCuts);
  die->GetTrackFilter().AddCuts(pidCuts);
  
  cuts->Print();
  //die->GetTrackFilter().AddCuts(cuts); // done now separately
}

//______________________________________________________________________________________
void SetupAODtrackCutsRun3_FT2(AliDielectron *die)
{
  //
  // Setup the track cuts
  // - these are similar cuts used for nano AOD creation LHC15o (ConfigLMEE_nano_PbPb2015.C), but tighter
  // - no ITS PID
  //

  Bool_t hasMC=die->GetHasMC();

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts    = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronVarCuts *varCuts2   = new AliDielectronVarCuts("VarCuts2","VarCuts2");
  AliDielectronTrackCuts *trkCuts  = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  
  // specific cuts
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kAtLeast, AliDielectronTrackCuts::kSPD0); // SPD first
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,      70.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0, 3.5);
  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 200.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8, 0.8);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0, 1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0, 3.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

 
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  // TPC
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
  // TOF
  pidCuts->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.2, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  //pidCutsMC->SetCutOnMCtruth(kTRUE);//un-comment to use only positive labels
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11.1, -10.9);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +10.9, +11.1);

  // activate the cut sets (order might be CPU timewise important)
  if(hasMC){
    cuts->AddCut(pidCutsMC); // commented in for the moment -> only true electrons
    die->GetTrackFilter().AddCuts(pidCutsMC);
  }
  cuts->AddCut(trkFilter);
  die->GetTrackFilter().AddCuts(trkFilter);
  cuts->AddCut(varCuts);
  die->GetTrackFilter().AddCuts(varCuts);
  cuts->AddCut(trkCuts);
  die->GetTrackFilter().AddCuts(trkCuts);
  cuts->AddCut(varCuts2);
  die->GetTrackFilter().AddCuts(varCuts2);
  cuts->AddCut(pidCuts);
  die->GetTrackFilter().AddCuts(pidCuts);
  
  cuts->Print();
  //die->GetTrackFilter().AddCuts(cuts); // done now separately
}

//______________________________________________________________________________________
void SetupAODtrackCutsRun3_FT2_ESD(AliDielectron *die)
{
  //
  // Setup the track cuts
  // - these are similar cuts used for nano AOD creation LHC15o (ConfigLMEE_nano_PbPb2015.C), but tighter
  // - no ITS PID
  // - ESD track cuts

  Bool_t hasMC=die->GetHasMC();

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  
  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts;
  fesdTrackCuts->SetTitle("ESD track cuts");
  
  // original upgrade study:
  fesdTrackCuts->SetDCAToVertex2D(kTRUE);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.0);
  fesdTrackCuts->SetMaxDCAToVertexXY(1.0);
  fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  
  // ITS:
  fesdTrackCuts->SetRequireITSRefit(kTRUE);
  fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  fesdTrackCuts->SetMinNClustersITS(3);
  // // TPC:
  fesdTrackCuts->SetRequireTPCRefit(kTRUE);
  fesdTrackCuts->SetMinNClustersTPC(70);	
  fesdTrackCuts->SetMaxChi2PerClusterTPC(3.5);

   /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts    = new AliDielectronVarCuts("VarCuts","VarCuts");
  
  // kinematic cuts
  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 200.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8, 0.8);
  

   /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  // TPC
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
  // TOF
  pidCuts->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.2, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
  

   /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  //pidCutsMC->SetCutOnMCtruth(kTRUE);//un-comment to use only positive labels
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11.1, -10.9);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +10.9, +11.1);

  // activate the cut sets (order might be CPU timewise important)
  if(hasMC){
    cuts->AddCut(pidCutsMC); // commented in for the moment -> only true electrons
    die->GetTrackFilter().AddCuts(pidCutsMC);
  }
  cuts->AddCut(varCuts);
  die->GetTrackFilter().AddCuts(varCuts);
  cuts->AddCut(fesdTrackCuts);
  die->GetTrackFilter().AddCuts(fesdTrackCuts);
  cuts->AddCut(pidCuts);
  die->GetTrackFilter().AddCuts(pidCuts);
  
  cuts->Print();
  //die->GetTrackFilter().AddCuts(cuts); // done now separately
}

//______________________________________________________________________________________
void SetupAODtrackCutsRun3_FT2_noPID(AliDielectron *die)
{
  //
  // Setup the track cuts
  // - these are similar cuts used for nano AOD creation LHC15o (ConfigLMEE_nano_PbPb2015.C), but tighter
  // - no PID
  //

  Bool_t hasMC=die->GetHasMC();

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts    = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronVarCuts *varCuts2   = new AliDielectronVarCuts("VarCuts2","VarCuts2");
  AliDielectronTrackCuts *trkCuts  = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  
  // specific cuts
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kAtLeast, AliDielectronTrackCuts::kSPD0); // SPD first
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,      70.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0, 3.5);
  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 200.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8, 0.8);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0, 1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0, 3.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

 
   /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  //pidCutsMC->SetCutOnMCtruth(kTRUE);//un-comment to use only positive labels
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11.1, -10.9);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +10.9, +11.1);

  // activate the cut sets (order might be CPU timewise important)
  if(hasMC){
    cuts->AddCut(pidCutsMC); // commented in for the moment -> only true electrons
    die->GetTrackFilter().AddCuts(pidCutsMC);
  }
  cuts->AddCut(trkFilter);
  die->GetTrackFilter().AddCuts(trkFilter);
  cuts->AddCut(varCuts);
  die->GetTrackFilter().AddCuts(varCuts);
  cuts->AddCut(trkCuts);
  die->GetTrackFilter().AddCuts(trkCuts);
  cuts->AddCut(varCuts2);
  die->GetTrackFilter().AddCuts(varCuts2);
  
  cuts->Print();
  //die->GetTrackFilter().AddCuts(cuts); // done now separately
}

//______________________________________________________________________________________
void SetupAODtrackCutsRun3_FT2_noPID_ESD(AliDielectron *die)
{
  //
  // Setup the track cuts
  // - these are similar cuts used for nano AOD creation LHC15o (ConfigLMEE_nano_PbPb2015.C), but tighter
  // - no PID
  // - ESD track cuts
  //

  Bool_t hasMC=die->GetHasMC();

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  
  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts;
  fesdTrackCuts->SetTitle("ESD track cuts");
  
  // original upgrade study:
  fesdTrackCuts->SetDCAToVertex2D(kTRUE);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.0);
  fesdTrackCuts->SetMaxDCAToVertexXY(1.0);
  fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  
  // ITS:
  fesdTrackCuts->SetRequireITSRefit(kTRUE);
  fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  fesdTrackCuts->SetMinNClustersITS(3);
  // // TPC:
  fesdTrackCuts->SetRequireTPCRefit(kTRUE);
  fesdTrackCuts->SetMinNClustersTPC(70);	
  fesdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts    = new AliDielectronVarCuts("VarCuts","VarCuts");
  
  // kinematic cuts
  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 200.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8, 0.8);


   /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  //pidCutsMC->SetCutOnMCtruth(kTRUE);//un-comment to use only positive labels
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11.1, -10.9);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +10.9, +11.1);

  // activate the cut sets (order might be CPU timewise important)
  if(hasMC){
    cuts->AddCut(pidCutsMC); // commented in for the moment -> only true electrons
    die->GetTrackFilter().AddCuts(pidCutsMC);
  }
  cuts->AddCut(varCuts);
  die->GetTrackFilter().AddCuts(varCuts);
  cuts->AddCut(fesdTrackCuts);
  die->GetTrackFilter().AddCuts(fesdTrackCuts);
  
  cuts->Print();
  //die->GetTrackFilter().AddCuts(cuts); // done now separately
}

//______________________________________________________________________________________
void SetupAODtrackCutsCarsten1(AliDielectron *die)
{
  //
  // Setup the track cuts - does not include Carstens's conversion rejection relevant cuts
  // - this is similar cut set than one exemplary one used by Carsten, including nanoAOD-prefiltering cuts
  //
  
  //Cuts used for nano AOD filtering (taken from ConfigLMEE_nano_PbPb2015.C on 22082018) - redundant for nano AODs but necessary for MC

  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");    
//  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany); // I think we loose the possibility to use prefilter?

  AliDielectronVarCuts *varCutsFilter   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCutsFilter = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  // specific cuts
  trkCutsFilter->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkCutsFilter->SetRequireITSRefit(kTRUE);
  trkCutsFilter->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts - nanoAOD filtering
  varCutsFilter->AddCut(AliDielectronVarManager::kPt,           0.2, 8.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
//  varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster    // did not work on ESD when filtering nanoAODs
  varCutsFilter->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

  //PID Cuts as used in nanoAOD Filtering
//  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
//  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);
//  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
//  pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-4.,4.);

  //Carsten's additional PID cuts: (Physics Forum 12.04.18)
  AliDielectronPID *pidCut_3        = new AliDielectronPID("PIDCuts","PIDCuts");
  pidCut_3->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2, 3.0 , 0. ,100., kFALSE);
  pidCut_3->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
  pidCut_3->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  pidCut_3->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  AliDielectronCutGroup* trackCuts=0x0;

  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");     
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
            
  AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
  double delta = 0.00001;
  AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
  trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
  AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
  trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
  AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
  trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
  AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
  trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
  AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
  trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
  AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
  trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
  SharedClusterCut->AddCut(trackCutsSharedCluster0);
  SharedClusterCut->AddCut(trackCutsSharedCluster2);
  SharedClusterCut->AddCut(trackCutsSharedCluster4);
  SharedClusterCut->AddCut(trackCutsSharedCluster8);
  SharedClusterCut->AddCut(trackCutsSharedCluster16);
  SharedClusterCut->AddCut(trackCutsSharedCluster32);
  
  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);
  trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
    
  trackCuts = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
  //Add nanoAOD prefiltering cuts
  
  //cuts->AddCut(trkFilter);    //not used on ESDs in ConfigLMEE_nano_PbPb2015.C
  cuts->AddCut(trackCuts);
  cuts->AddCut(pidCut_3);
  cuts->Print();
  
  
  trackCuts->AddCut(varCutsFilter);
  trackCuts->AddCut(trkCutsFilter);
  
  trackCuts->AddCut(trackCutsDiel);
  trackCuts->AddCut(trackCutsAOD);

  trackCuts->AddCut(pidCut_3);
  
  
  die->GetTrackFilter().AddCuts(trackCuts);
}

void SetupAODtrackCutsCarsten2(AliDielectron *die)
{
  //
  // Setup the track cuts - includes Carstens's conversion rejection relevant cuts
  // - this is similar cut set than one exemplary one used by Carsten, including nanoAOD-prefiltering cuts
  //
  
  //Cuts used for nano AOD filtering (taken from ConfigLMEE_nano_PbPb2015.C on 22082018) - redundant for nano AODs but necessary for MC

  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");    
//  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany); // I think we loose the possibility to use prefilter?

  AliDielectronVarCuts *varCutsFilter   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCutsFilter = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  // specific cuts
  trkCutsFilter->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkCutsFilter->SetRequireITSRefit(kTRUE);
  trkCutsFilter->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts - nanoAOD filtering
  varCutsFilter->AddCut(AliDielectronVarManager::kPt,           0.2, 8.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster    // did not work on ESD when filtering nanoAODs
  varCutsFilter->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

  //PID Cuts as used in nanoAOD Filtering
//  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
//  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);
//  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
//  pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-4.,4.);

  //Carsten's additional PID cuts: (Physics Forum 12.04.18)
  AliDielectronPID *pidCut_3        = new AliDielectronPID("PIDCuts","PIDCuts");
  pidCut_3->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2, 3.0 , 0. ,100., kFALSE);
  pidCut_3->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
  pidCut_3->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  pidCut_3->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  AliDielectronCutGroup* trackCuts=0x0;

  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");     
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
            
  AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
  double delta = 0.00001;
  AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
  trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
  AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
  trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
  AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
  trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
  AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
  trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
  AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
  trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
  AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
  trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
  SharedClusterCut->AddCut(trackCutsSharedCluster0);
  SharedClusterCut->AddCut(trackCutsSharedCluster2);
  SharedClusterCut->AddCut(trackCutsSharedCluster4);
  SharedClusterCut->AddCut(trackCutsSharedCluster8);
  SharedClusterCut->AddCut(trackCutsSharedCluster16);
  SharedClusterCut->AddCut(trackCutsSharedCluster32);
  
  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);
  trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
    
  trackCuts = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
  //Add nanoAOD prefiltering cuts
  
  //cuts->AddCut(trkFilter);    //not used on ESDs in ConfigLMEE_nano_PbPb2015.C
  cuts->AddCut(trackCuts);
  cuts->AddCut(pidCut_3);
  cuts->Print();
  
  
  trackCuts->AddCut(varCutsFilter);
  trackCuts->AddCut(trkCutsFilter);
  
  trackCuts->AddCut(trackCutsDiel);
  trackCuts->AddCut(trackCutsAOD);
  trackCuts->AddCut(SharedClusterCut);

  trackCuts->AddCut(pidCut_3);
  
  
  die->GetTrackFilter().AddCuts(trackCuts);
}

void SetTPCCorr(AliDielectron *die){
  ::Info("Config_miweber_LMEE_PbPb_woCutLib","starting LMEECutLib::SetEtaCorrectionTPC()\n");
  TString path="alien:///alice/cern.ch/user/s/selehner/recal/recalib_data_tpc_nsigmaele.root";
  gSystem->Exec(TString::Format("alien_cp %s .",path.Data()));
  ::Info("Config_miweber_LMEE_PbPb_woCutLib","Copy TPC correction from Alien: %s",path.Data());
  TFile *infile = TFile::Open("recalib_data_tpc_nsigmaele.root");
  
  TH3D* mean = dynamic_cast<TH3D*>(infile->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(infile->Get("sum_width_correction"));

  if(mean)   ::Info("Config_miweber_LMEE_PbPb_woCutLib","Mean Correction Histo loaded, entries: %f",mean->GetEntries());
  else {
    ::Info("Config_miweber_LMEE_PbPb_woCutLib","Mean Correction Histo not loaded! entries: %f",mean->GetEntries());
    return 0;
  }

    die->SetCentroidCorrFunction(mean, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    die->SetWidthCorrFunction(width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
       
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

  //add MC signal histograms to pair class
  if(die->GetMCSignals()) {
    for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) {
      histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
      // histos->AddClass(Form("Track_%s",die->GetMCSignals()->At(i)->GetName()));
      // histos->AddClass(Form("Track_Legs_%s",die->GetMCSignals()->At(i)->GetName()));
      // histos->AddClass(Form("Pair_%s_MCtruth",die->GetMCSignals()->At(i)->GetName()));
    }
  }
  
  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0,1,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","ZVertex","ZVertex;ZVertex/cm",120,-12.,12.,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","Centrality","Centrality;Centrality/%",202,-1.,100.,AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","CentralityNew","Centrality;Centrality/%",202,-1.,100.,AliDielectronVarManager::kCentralityNew);
  histos->UserHistogram("Event","nEvTPC_eventplaneents",";;ev plane;",AliDielectronHelper::MakeLinBinning(180,  TMath::Pi()/-2.,TMath::Pi()/2.),AliDielectronVarManager::kQnTPCrpH2);
  histos->UserHistogram("Event","nEvZDCA_eventplaneents",";;ev plane;",AliDielectronHelper::MakeLinBinning(180,  TMath::Pi()/-2.,TMath::Pi()/2.),AliDielectronVarManager::kQnZDCArpH1);
  histos->UserHistogram("Event","nEvZDCC_eventplaneents",";;ev plane;",AliDielectronHelper::MakeLinBinning(180,  TMath::Pi()/-2.,TMath::Pi()/2.),AliDielectronVarManager::kQnZDCCrpH1);
  
  histos->UserHistogram("Event","nEvZDCAX_Centrality",";;Centrality;ZDCAX;",202,-1.,100., 200, -1.,1.,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kQnZDCAX);
  histos->UserHistogram("Event","nEvZDCAY_Centrality",";;Centrality;ZDCAY;",202,-1.,100., 200, -1.,1.,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kQnZDCAY);
   histos->UserHistogram("Event","nEvZDCCX_Centrality",";;Centrality;ZDCCX;",202,-1.,100., 200, -1.,1.,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kQnZDCCX);
  histos->UserHistogram("Event","nEvZDCCY_Centrality",";;Centrality;ZDCCY;",202,-1.,100., 200, -1.,1.,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kQnZDCCY);

  histos->UserHistogram("Event","nEvZDCAX_Vx",";;Vx;ZDCAX;",200,0.06,0.09,200, -1.,1.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kQnZDCAX);
  histos->UserHistogram("Event","nEvZDCAY_Vx",";;Vx;ZDCAY;",200,0.06,0.09,200, -1.,1.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kQnZDCAY);
   histos->UserHistogram("Event","nEvZDCCX_Vx",";;Vx;ZDCCX;",200,0.06,0.09,200, -1.,1.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kQnZDCCX);
  histos->UserHistogram("Event","nEvZDCCY_Vx",";;Vx;ZDCCY;",200,0.06,0.09,200, -1.,1.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kQnZDCCY);

  histos->UserHistogram("Event","nEvZDCAX_Vy",";;Vy;ZDCAX;",200,0.32,0.35,200, -1.,1.,AliDielectronVarManager::kYvPrim,AliDielectronVarManager::kQnZDCAX);
  histos->UserHistogram("Event","nEvZDCAY_Vy",";;Vy;ZDCAY;",200,0.32,0.35,200, -1.,1.,AliDielectronVarManager::kYvPrim,AliDielectronVarManager::kQnZDCAY);
   histos->UserHistogram("Event","nEvZDCCX_Vy",";;Vy;ZDCCX;",200,0.32,0.35,200, -1.,1.,AliDielectronVarManager::kYvPrim,AliDielectronVarManager::kQnZDCCX);
  histos->UserHistogram("Event","nEvZDCCY_Vy",";;Vy;ZDCCY;",200,0.32,0.35,200, -1.,1.,AliDielectronVarManager::kYvPrim,AliDielectronVarManager::kQnZDCCY);

  histos->UserHistogram("Event","nEvZDCAX_Vz",";;Vz;ZDCAX;",200,-10.,10.,200, -1.,1.,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kQnZDCAX);
  histos->UserHistogram("Event","nEvZDCAY_Vz",";;Vz;ZDCAY;",200,-10.,10.,200, -1.,1.,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kQnZDCAY);
   histos->UserHistogram("Event","nEvZDCCX_Vz",";;Vz;ZDCCX;",200,-10.,10.,200, -1.,1.,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kQnZDCCX);
  histos->UserHistogram("Event","nEvZDCCY_Vz",";;Vz;ZDCCY;",200,-10.,10.,200, -1.,1.,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kQnZDCCY);
  

  
  histos->UserHistogram("Event","NTPCclsEvent","kNTPCclsEvent;Centrality/%;kNTPCclsEvent",202,-1.,100.,500,0,5000000,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kNTPCclsEvent);
  histos->UserHistogram("Event","NTPCtrkswITSout","kNTPCtrkswITSout;Centrality/%;kNTPCtrkswITSout",202,-1.,100.,500,0,15000,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kNTPCtrkswITSout);


  //add histograms to track class
  histos->UserHistogram("Track","Centrality","Centrality;Centrality (%);#tracks",100,0.,100.,AliDielectronVarManager::kCentralityNew);
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
			"InvMass_QnTPCrpH2","Inv.Mass:QnTPCrpH2;Inv. Mass (GeV/c^{2}); QnTPCrpH2",
			500,0.,5.,200,-TMath::Pi(),TMath::Pi(),
			AliDielectronVarManager::kM, AliDielectronVarManager::kQnTPCrpH2);
  
  histos->UserHistogram("Pair", 
			"InvMass_PairPlaneMagInProZDC","Inv.Mass:PairPlaneMagInProZDC;Inv. Mass (GeV/c^{2}); kPairPlaneMagInProZDC",
			500,0.,5.,200,-1,1,
			AliDielectronVarManager::kM, AliDielectronVarManager::kPairPlaneMagInProZDC);
  
  histos->UserHistogram("Pair", 
			"InvMass_QnTPCDelta","Inv.Mass:QnTPC;Inv. Mass (GeV/c^{2}); kQnDeltaPhiTPCrpH2",
			500,0.,5.,200,-TMath::Pi(),TMath::Pi(),
			AliDielectronVarManager::kM, AliDielectronVarManager::kQnDeltaPhiTPCrpH2);

  histos->UserHistogram("Pair", 
			"InvMass_QnTPCSPV2","Inv.Mass:QnTPCSPV2;Inv. Mass (GeV/c^{2}); kQnTPCrpH2FlowSPV2",
			500,0.,5.,200,-1.,1.,
			AliDielectronVarManager::kM, AliDielectronVarManager::kQnTPCrpH2FlowSPV2);

  histos->UserHistogram("Pair", 
			"InvMass_CosThetaCS","Inv.Mass:CosThetaCS;Inv. Mass (GeV/c^{2}); cos(#theta_{CS})",
			500,0.,5.,200,-1.,1.,
			AliDielectronVarManager::kM, AliDielectronVarManager::kThetaCS);

  histos->UserHistogram("Pair", 
			"InvMass_PhiCS","Inv.Mass:PhiCS;Inv. Mass (GeV/c^{2}); #varphi_{CS}",
			500,0.,5.,200,-TMath::Pi(),TMath::Pi(),
			AliDielectronVarManager::kM, AliDielectronVarManager::kPhiCS);
  
   histos->UserHistogram("Pair", 
                         "InvMass_Centrality","Inv.Mass:Centrality;Inv. Mass (GeV/c^{2}); Centrality (%)",
                         500,0.,5.,10,0.,100.,
                         AliDielectronVarManager::kM, AliDielectronVarManager::kCentralityNew);
   
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

AliDielectronEventCuts *GetEventCuts(Double_t centMin, Double_t centMax, Bool_t reqAliEventCuts, Bool_t reqAliEventCutsCorrelated){

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex SPD && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny); // AOD
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1);

  //centrality cuts if required
  if(centMin > -1 && centMax > -1 && centMin < centMax){
    Bool_t isRun2 = kTRUE;
    eventCuts->SetCentralityRange(centMin,centMax,isRun2);
  }

  // AliEventCuts if required
  if(reqAliEventCuts){
    eventCuts->SetRequireAliEventCuts(reqAliEventCuts,reqAliEventCutsCorrelated);
  }
  
  return eventCuts;
}


void SetupMCsignals(AliDielectron* die){

  Printf("Setting up MC signals...");

  // ##################### "real" pairs from signals (pi0,eta,eta',rho, omega, phi) ##############################
   AliDielectronSignalMC* pi0Sig = new AliDielectronSignalMC("pi0", "pi0Signal"); ///pi0 dalitz pairs 
  pi0Sig->SetLegPDGs(11,-11);
  pi0Sig->SetMotherPDGs(111,111);
  pi0Sig->SetMothersRelation(AliDielectronSignalMC::kSame);
  pi0Sig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  pi0Sig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  pi0Sig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pi0Sig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  pi0Sig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(pi0Sig);

  AliDielectronSignalMC* pi0All = new AliDielectronSignalMC("pi0All", "pi0All"); ///pi0 dalitz pairs (also from secondary)
  pi0All->SetLegPDGs(11,-11);
  pi0All->SetMotherPDGs(111,111);
  pi0All->SetMothersRelation(AliDielectronSignalMC::kSame);
  pi0All->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  pi0All->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pi0All->SetCheckBothChargesMothers(kTRUE,kTRUE);
  pi0All->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(pi0All);


  AliDielectronSignalMC* etaSig = new AliDielectronSignalMC("Eta", "etaSignal"); ///eta dalitz pairs 
  etaSig->SetLegPDGs(11,-11);
  etaSig->SetMotherPDGs(221,221);
  etaSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  etaSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  etaSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  etaSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  etaSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  etaSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(etaSig);


  AliDielectronSignalMC* etaprimeSig = new AliDielectronSignalMC("Etaprime", "etaprimeSignal"); ///etaprime pairs 
  etaprimeSig->SetLegPDGs(11,-11);
  etaprimeSig->SetMotherPDGs(331,331);
  etaprimeSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  etaprimeSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  etaprimeSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  etaprimeSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  etaprimeSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  etaprimeSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(etaprimeSig);


  AliDielectronSignalMC* rhoSig = new AliDielectronSignalMC("Rho", "rhoSignal"); ///rho pairs 
  rhoSig->SetLegPDGs(11,-11);
  rhoSig->SetMotherPDGs(113,113);
  rhoSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  rhoSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  rhoSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  rhoSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  rhoSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  rhoSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(rhoSig);

  AliDielectronSignalMC* omegaSig = new AliDielectronSignalMC("Omega", "omegaSignal"); ///omega pairs 
  omegaSig->SetLegPDGs(11,-11);
  omegaSig->SetMotherPDGs(223,223);
  omegaSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  omegaSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  omegaSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  omegaSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  omegaSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  omegaSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(omegaSig);
  
  AliDielectronSignalMC* phiSig = new AliDielectronSignalMC("Phi", "phiSignal"); ///phi pairs 
  phiSig->SetLegPDGs(11,-11);
  phiSig->SetMotherPDGs(333,333);
  phiSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  phiSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  phiSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  phiSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  phiSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  phiSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(phiSig);
  
  // ##################### "real" pairs from photon conversions in the detector material ##############################
  AliDielectronSignalMC* signalFromResonance_ULS_gammaConv = new AliDielectronSignalMC("signalFromResonance_ULS_gammaConv", "signalFromResonance_ULS_gammaConv");
  signalFromResonance_ULS_gammaConv->SetLegPDGs(11,-11);
  signalFromResonance_ULS_gammaConv->SetMotherPDGs(22,22);
  signalFromResonance_ULS_gammaConv->SetMothersRelation(AliDielectronSignalMC::kSame);
  signalFromResonance_ULS_gammaConv->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // kSecondary means decays in the detector
  signalFromResonance_ULS_gammaConv->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(signalFromResonance_ULS_gammaConv);

  // ##################### combinatorial pairs ##############################
  AliDielectronSignalMC* diEleCombinatiorial = new AliDielectronSignalMC("diEleCombinatiorial", "diEleCombinatiorial");
  diEleCombinatiorial->SetLegPDGs(11,-11);
  diEleCombinatiorial->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleCombinatiorial->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(diEleCombinatiorial);

  AliDielectronSignalMC* diEleCombinatiorialConversion = new AliDielectronSignalMC("diEleCombinatiorialConversion", "diEleCombinatiorialConversion");
  diEleCombinatiorialConversion->SetLegPDGs(11,-11);
  diEleCombinatiorialConversion->SetMotherPDGs(22,0);// 1 leg from photons + 1 leg from everything
  diEleCombinatiorialConversion->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleCombinatiorialConversion->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(diEleCombinatiorialConversion);

    // ##################### HF pairs ##############################
  AliDielectronSignalMC* diEleHF = new AliDielectronSignalMC("diEleHF", "diEleHF");
  diEleHF->SetLegPDGs(11,-11);
  diEleHF->SetMotherPDGs(401,401);
  diEleHF->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleHF->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(diEleHF);

 
  // #################### D-Mesons
  AliDielectronSignalMC* diEleOpenCharmCharged = new AliDielectronSignalMC("DmesonsCharged","di-electrons from open charm D+- mesons no B grandmother");  // dielectrons originating from open charm hadrons
  diEleOpenCharmCharged->SetLegPDGs(11,-11);
  diEleOpenCharmCharged->SetMotherPDGs(401,401);
  diEleOpenCharmCharged->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenCharmCharged->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  // diEleOpenCharmCharged->SetFillPureMCStep(kTRUE);
  diEleOpenCharmCharged->SetCheckStackForPDG(kTRUE);
  diEleOpenCharmCharged->SetPDGforStack(503);
  diEleOpenCharmCharged->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharmCharged->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenCharmCharged);

  AliDielectronSignalMC* diEleOpenCharmNeutral = new AliDielectronSignalMC("DmesonsNeutral","di-electrons from open charm D0 mesons no B grandmother");  // dielectrons originating from open charm hadrons
  diEleOpenCharmNeutral->SetLegPDGs(11,-11);
  diEleOpenCharmNeutral->SetMotherPDGs(405,405);
  diEleOpenCharmNeutral->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenCharmNeutral->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenCharmNeutral->SetCheckStackForPDG(kTRUE);
  diEleOpenCharmNeutral->SetPDGforStack(503);
  diEleOpenCharmNeutral->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharmNeutral->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenCharmNeutral);

  //B meson (3)
  AliDielectronSignalMC* diEleOneOpenB = new AliDielectronSignalMC("B2ee","di-electrons from one B meson");  // dielectrons originating from open charm hadrons
  diEleOneOpenB->SetLegPDGs(11,-11);
  diEleOneOpenB->SetMotherPDGs(401,501);
  diEleOneOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOneOpenB->SetGrandMotherPDGs(501,0);
  diEleOneOpenB->SetCheckMotherGrandmotherRelation(kTRUE,kTRUE);
  diEleOneOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOneOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOneOpenB->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOneOpenB);

  //B meson (1)(1)
  AliDielectronSignalMC* diEleOpenB = new AliDielectronSignalMC("BMesons","di-electrons from B mesons");  // dielectrons originating from open charm hadrons
  diEleOpenB->SetLegPDGs(11,-11);
  diEleOpenB->SetMotherPDGs(501,501);
  diEleOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenB->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenB);

  //B meson (2)(2)
  AliDielectronSignalMC* diEleOpenBtoD = new AliDielectronSignalMC("B2D2ee","di-electrons from B->D-> e");  // dielectrons originating from open charm hadrons
  diEleOpenBtoD->SetLegPDGs(11,-11);
  diEleOpenBtoD->SetMotherPDGs(401,401);
  diEleOpenBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenBtoD->SetGrandMotherPDGs(501,501);
  diEleOpenBtoD->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenBtoD);

  //B meson (1)(2)
  AliDielectronSignalMC* diEleOpenBandBtoD = new AliDielectronSignalMC("B2eAndB2D2e","di-electrons from B->e and B->D->e");  // dielectrons originating from open charm hadrons
  diEleOpenBandBtoD->SetLegPDGs        (11,11);
  diEleOpenBandBtoD->SetMotherPDGs     (401,501);
  diEleOpenBandBtoD->SetGrandMotherPDGs(501,0);
  diEleOpenBandBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenBandBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBandBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBandBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  //do i need this?
  diEleOpenBandBtoD->SetCheckMotherGrandmotherRelation(kTRUE,kFALSE);
  die->AddSignalMC(diEleOpenBandBtoD);
}

