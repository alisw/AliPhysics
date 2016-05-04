void InitHistograms(AliDielectron *die, ULong_t cutDefinition);
void InitCF(AliDielectron* die, ULong_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, ULong_t cutDefinition);
void SetupPairCuts(AliDielectron *die, ULong_t cutDefinition);
void SetupV0Cuts(AliDielectron *die, ULong_t cutDefinition, ULong_t particle = 0, Bool_t exclude = kTRUE);
void SetEtaCorrection();
void SetEtaCorrectionGSI();




int i=-1;
ULong_t kSpdFirst = 1 << ++i;
ULong_t kSpdAny = 1 << ++i;
ULong_t kItsRefit = 1 << ++i;
ULong_t kTrackPt = 1 << ++i;
ULong_t kTpcInclEl = 1 << ++i;
ULong_t kTpcStrictInclEl = 1 << ++i;
ULong_t kTpcExclPion = 1 << ++i;
ULong_t kTpcStrictExclPion = 1 << ++i;
ULong_t kTpcExclProton = 1 << ++i;
ULong_t kTpcStrictExclProton = 1 << ++i;
ULong_t kCombinedTpcTof = 1 << ++i;
ULong_t kCombinedTpcTof2 = 1 << ++i;
ULong_t kCombinedTpcTof3 = 1 << ++i;
ULong_t kCombinedTpcTof4 = 1 << ++i;
ULong_t kTpcExclKaon = 1 << ++i;
ULong_t kTofInclEl = 1 << ++i;
ULong_t kTofExclProton = 1 << ++i;
ULong_t kYcut = 1 << ++i;
ULong_t kYcutStrict = 1 << ++i;
ULong_t kV0electrons = 1 << ++i;
ULong_t kV0pions = 1 << ++i;
ULong_t kTpcInclPion = 1 << ++i;
ULong_t kTofInclPion = 1 << ++i;
ULong_t kTpcExclEl = 1 << ++i;
ULong_t kMCelectrons = 1<< ++i;
ULong_t kMCpions = 1<< ++i;
ULong_t kActiveLengthCut = 1 << ++i;
ULong_t kStrictActiveLengthCut = 1 << ++i;
ULong_t kStrictDcaCut = 1 << ++i;

ULong_t kNoPairing = 1 << ++i;



ULong_t kTpc = kTpcInclEl | kTpcExclPion | kTpcExclProton;
ULong_t kPid = kTpc | kTofInclEl;
ULong_t kStandard = kTpc | kItsRefit | kSpdAny | kTrackPt;
ULong_t kStrict = kTpcStrictInclEl | kTpcStrictExclPion | kTpcStrictExclProton | kItsRefit | kSpdAny | kTrackPt;
ULong_t kSmart = kTpcInclEl | kTpcExclPion | kCombinedTpcTof | kItsRefit | kSpdAny | kTrackPt;
ULong_t kSmart2 = kTpcInclEl | kTpcExclPion | kCombinedTpcTof2 | kItsRefit | kSpdAny | kTrackPt;
ULong_t kSmart3 = kTpcInclEl | kTpcExclPion | kCombinedTpcTof3 | kItsRefit | kSpdAny | kTrackPt;
ULong_t kSmart4 = kTpcInclEl | kTpcExclPion | kCombinedTpcTof4 | kItsRefit | kSpdAny | kTrackPt;



ULong_t kPions = kTpcInclPion |kTofInclPion|kTpcExclEl;


Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
enum { k10b=0, k10c, k10d, k10e, k10f, k10h, k10pp, k11a, k11d, k11h, k12h, k13b, k13c, k13d, k13e, k13f, k15f, k15h };


void ConfigJpsi( AliAnalysisTaskMultiDielectron* task )
{
  
  
  Bool_t usePhysicsSelection = kTRUE;
  Bool_t rejectPileup = kTRUE;
  Bool_t triggerOnV0AND = kFALSE;
  Bool_t setTriggerMaskFromPeriod = kFALSE;
  AliVEvent::EOfflineTriggerTypes triggerMask = AliVEvent::kMB;
  TString triggerClass = "";
  
  
  
  vector< unsigned long > cutsToUse;
  if(!hasMC){

     cutsToUse.push_back( kStandard  );
     cutsToUse.push_back( kStandard | kStrictActiveLengthCut );
   
     cutsToUse.push_back( kStrict  );
     cutsToUse.push_back( kStrict | kStrictActiveLengthCut );
  }

  TString names;
  
  if(hasMC){
    names = "posEta;negEta";
  }
  else{
    names = "StandardCuts;"
            "StandardCuts+StrictActiveLength;"
            "StrictCuts;"
            "StrictCuts+StrictActiveLength;";
  }
  
  TObjArray *arrNames=names.Tokenize(";");

  
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetRequireCompleteDAQ(kTRUE);
  task->SetEventFilter(eventCuts);
  
   
  if( setTriggerMaskFromPeriod ){

    // Get the current train configuration
    TString trainConfig = gSystem->Getenv("CONFIG_FILE");
    TString list = gSystem->Getenv("LIST");
    if( prod.IsNull()) prod=list;

    // selected period
    if(      prod.Contains("LHC10b") ) iPeriod = k10b;
    else if( prod.Contains("LHC10c") ) iPeriod = k10c;
    else if( prod.Contains("LHC10d") ) iPeriod = k10d;
    else if( prod.Contains("LHC10e") ) iPeriod = k10e;
    else if( prod.Contains("LHC10f") ) iPeriod = k10f;
    else if( prod.Contains("LHC10h") ) iPeriod = k10h;
    else if( prod.Contains("LHC11a") ) iPeriod = k11a;
    else if( prod.Contains("LHC11d") ) iPeriod = k11d;
    else if( prod.Contains("LHC11h") ) iPeriod = k11h;
    else if( prod.Contains("LHC12h") ) iPeriod = k12h;
    else if( prod.Contains("LHC13b") ) iPeriod = k13b;
    else if( prod.Contains("LHC13c") ) iPeriod = k13c;
    else if( prod.Contains("LHC13d") ) iPeriod = k13d;
    else if( prod.Contains("LHC13e") ) iPeriod = k13e;
    else if( prod.Contains("LHC13f") ) iPeriod = k13f;
    else if( prod.Contains("LHC15f") ) iPeriod = k15f;
    else if( prod.Contains("LHC15h") ) iPeriod = k15h;



    
    // add special triggers
    switch(iPeriod) {
      case k11d: task->SetTriggerMask(AliVEvent::kEMCEJE+AliVEvent::kEMC7+AliVEvent::kEMCEGA);     break;
      case k11h: task->SetTriggerMask(AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral); break;
      case k12h: task->SetTriggerMask(AliVEvent::kAnyINT); break;                                      
      case k13b: task->SetTriggerMask(AliVEvent::kINT7); break;
      case k13c: task->SetTriggerMask(AliVEvent::kINT7); break;
      case k13d: task->SetTriggerMask(AliVEvent::kAnyINT); break;
      case k13e: task->SetTriggerMask(AliVEvent::kAnyINT); break;
      case k13f: task->SetTriggerMask(AliVEvent::kAnyINT); break;
      case k15f: task->SetTriggerMask(AliVEvent::kINT7); break;
      case k15h: task->SetTriggerMask(AliVEvent::kINT7); break;
    }
  }
  else{
    task->SetTriggerMask( triggerMask );
  }
  if(! triggerClass.IsNull() ) task->SetFiredTriggerName(triggerClass.Data() );
  if(triggerOnV0AND) task->SetTriggerOnV0AND();
  if(usePhysicsSelection)task->UsePhysicsSelection();
  if(rejectPileup) task->SetRejectPileup();

  
  
  for(ULong_t cutNumber=0; cutNumber < cutsToUse.size(); cutNumber++){
    ULong_t cutDefinition = cutsToUse[cutNumber];
    TString name=Form("%02d",cutDefinition);
    if (cutNumber < arrNames->GetEntriesFast()){
      name=arrNames->At(cutNumber)->GetName();
    }
    AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("Track cuts: %s",name.Data()));

    if (hasMC) SetupMCsignals(die);
    // cut setup
    SetupTrackCuts(die, cutDefinition);  
    
    if( !(cutDefinition & kNoPairing)   ) {
      SetupPairCuts(die,cutDefinition);
    }
    
    if( !( cutDefinition & ( kV0electrons | kV0pions) ) ) {
      SetupV0Cuts(die, cutDefinition);
    }   
    else if(  cutDefinition &  kV0electrons  ) {
      SetupV0Cuts(die, cutDefinition, 0, kFALSE);
    }
    else if(  cutDefinition &  kV0pions  ) {
      SetupV0Cuts(die, cutDefinition, 1, kFALSE);
    }
    
    InitHistograms(die,cutDefinition);
    if (hasMC) InitCF(die,cutDefinition);
    if ( cutDefinition & kNoPairing ) die->SetNoPairing();
    
    if(!hasMC &&  !( cutDefinition & kNoPairing  )  ){  
      AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
      mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-8,-5,0,5,8,10");
     // mix->AddVariable(AliDielectronVarManager::kPhi ,"0,3.14,6.3");
      mix->SetDepth(10);
      die->SetMixingHandler(mix);
      AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
      rot->SetConeAnglePhi( TMath::Pi() );
      rot->SetIterations(40);
      die->SetTrackRotator(rot);
    }
  //SetEtaCorrection();
    task->AddDielectron(die);
  }
}
  
/**
*  
*  Setup the track cuts according to given cut definition
*  
**/
void SetupTrackCuts(AliDielectron *die, ULong_t cutDefinition)
{
  
  
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);

  //default quality cuts
  AliDielectronTrackCuts *cut1=new AliDielectronTrackCuts("cut1","cut1");
  if(cutDefinition & kItsRefit){
    cut1->SetRequireITSRefit(kTRUE);
  }
  cut1->SetRequireTPCRefit(kTRUE);
  if( cutDefinition & kSpdFirst ){
    cut1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  }
  else if(cutDefinition & kSpdAny){
    cut1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }
  cuts->AddCut(cut1);
  
  
// track quality cuts
  
  AliDielectronVarCuts *trackQuality = new AliDielectronVarCuts("trackQuality1","trackQuality1");
  
  if(cutDefinition & kStrictDcaCut){
    trackQuality->AddCut(AliDielectronVarManager::kImpactParZ, -1. ,  1. );
    trackQuality->AddCut(AliDielectronVarManager::kImpactParXY,-0.3 ,   0.3 );
    
  }
  else{
    trackQuality->AddCut(AliDielectronVarManager::kImpactParZ, -3. ,  3. );
    trackQuality->AddCut(AliDielectronVarManager::kImpactParXY,-1. ,   1. );
  }
  trackQuality->AddCut(AliDielectronVarManager::kEta,        - .9,   .9);
  trackQuality->AddCut(AliDielectronVarManager::kKinkIndex0,   0.       );
  
  if(cutDefinition & kActiveLengthCut){
    trackQuality->AddCut(AliDielectronVarManager::kTPCGeomLength,  1.02, 2. );
  }
  else if(cutDefinition & kStrictActiveLengthCut){
    trackQuality->AddCut(AliDielectronVarManager::kTPCGeomLength,     1.23, 2. );
  }
  else{
    trackQuality->AddCut(AliDielectronVarManager::kNclsTPC,     70., 160. );
  }
  trackQuality->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.,   4. );
  
  
  if(cutDefinition & kTrackPt ){
    trackQuality->AddCut(AliDielectronVarManager::kPt,         1. ,  1e30);
  }
  if(cutDefinition & kYcut){
    trackQuality->AddCut(AliDielectronVarManager::kYsignedIn, -13., 14.);
  }
  else if(cutDefinition & kYcutStrict){
    trackQuality->AddCut(AliDielectronVarManager::kYsignedIn, -13., 5.);
  }
  
  /*
  if(cutDefinition & kPosEta ){
    trackQuality->AddCut(AliDielectronVarManager::kEta,        0. ,  1.);
  }
  
  
  if(cutDefinition & kNegEta ){
    trackQuality->AddCut(AliDielectronVarManager::kEta,        -1. ,  0.);
  }
  */
  
  
  
  if(cutDefinition & kMCelectrons || cutDefinition & kMCpions){
    AliDielectronCutGroup * mcCut = new AliDielectronCutGroup("mcTruth", "mcTruth");
    if(cutDefinition & kMCelectrons){
      AliDielectronVarCuts *mcEln = new AliDielectronVarCuts("MCelectrons-","MCelectrons-");
      mcEln->AddCut( AliDielectronVarManager::kPdgCode,  11., 11.);
      AliDielectronVarCuts *mcElp = new AliDielectronVarCuts("MCelectrons+","MCelectrons+");
      mcElp->AddCut( AliDielectronVarManager::kPdgCode,  -11., -11.);
      mcCut->AddCut(mcEln);
      mcCut->AddCut(mcElp);
    }  
    if(cutDefinition & kMCpions){
      AliDielectronVarCuts *mcPin = new AliDielectronVarCuts("MCpions-","MCpions-");
      mcPin->AddCut( AliDielectronVarManager::kPdgCode,  -211., -211.);
      AliDielectronVarCuts *mcPip = new AliDielectronVarCuts("MCpions+","MCpions+");
      mcPip->AddCut( AliDielectronVarManager::kPdgCode,  211., 211.);
      mcCut->AddCut(mcPin);
      mcCut->AddCut(mcPip);
    }
  
    cuts->AddCut(mcCut);
  }
  
  
  
  cuts->AddCut(trackQuality);
  
  
// PID cuts 
  
  AliDielectronPID *pid=new AliDielectronPID("PID","PID");
  
 
  
  if(cutDefinition & kTpcInclEl ){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3., 3.);
  }
  else if(cutDefinition & kTpcStrictInclEl ){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.);
  }
  else{
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -999., 999.);
  }
  if(cutDefinition & kTpcExclPion){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion, 3.5, 999.);
  }
  else if(cutDefinition & kTpcStrictExclPion){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion, 4., 999.);
  }
  if(cutDefinition & kTpcExclProton){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton, 3., 999.);
  }  
  else if(cutDefinition & kTpcStrictExclProton){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton, 3.5, 999.);
  }  
  if(cutDefinition & kTofInclEl){
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3., 3.);
  }
  if(cutDefinition & kTpcExclKaon){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kKaon, 3., 999.);
  }
  if(cutDefinition & kTpcInclPion){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -3., 3.);
  }
  if(cutDefinition & kTofInclPion){
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kPion, -3., 3.);
  }
  if(cutDefinition & kTofExclProton){
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kProton, -999., -3.);
  }
  if(cutDefinition & kTpcExclEl){
    pid->AddCut(AliDielectronPID::kTpc,AliPID::kElectron, -999., -3.);
  }
  if(cutDefinition & kCombinedTpcTof){
    AliDielectronCutGroup *pbandwithTOF = new AliDielectronCutGroup("pbandwithTOF", "pbandwithTOF", AliDielectronCutGroup::kCompOR);//NOTE taken out for testing

        // EITHER 1. electron candidates inside proton band and with electron TOFpid, rejected if outside 3sigma TOF
    AliDielectronPID *pidTOF = new AliDielectronPID("PIDTOFthere ","TOF nSigma |e|<3.0 + TPC nSigma |p|<3.");
        //TOF bit already internally required
    pidTOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.);
    pidTOF->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.5,3.5);
    pbandwithTOF->AddCut(pidTOF);

        // OR 2. electron candidates outside proton band accepted if no matching inside with TOF
    AliDielectronPID *pidproton = new AliDielectronPID("PIDTPCproton ","TPC nSigma |P|>3");
    pidproton->AddCut(AliDielectronPID::kTPC,AliPID::kProton,3.5, 999.);
    pbandwithTOF->AddCut(pidproton);
    cuts->AddCut(pbandwithTOF);
  }
  else if(cutDefinition & kCombinedTpcTof2){
    AliDielectronCutGroup *pbandwithTOF = new AliDielectronCutGroup("pbandwithTOF", "pbandwithTOF", AliDielectronCutGroup::kCompOR);//NOTE taken out for testing

        // EITHER 1. electron candidates inside proton band and with electron TOFpid, rejected if outside 3sigma TOF
    AliDielectronPID *pidTOF = new AliDielectronPID("PIDTOFthere ","TOF nSigma |e|<3.0 + TPC nSigma |p|<3.");
        //TOF bit already internally required
    pidTOF->AddCut(AliDielectronPID::kTOF,AliPID::kProton,-999.,-3.);
    pidTOF->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.5,3.5);
    pbandwithTOF->AddCut(pidTOF);

        // OR 2. electron candidates outside proton band accepted if no matching inside with TOF
    AliDielectronPID *pidproton = new AliDielectronPID("PIDTPCproton ","TPC nSigma |P|>3");
    pidproton->AddCut(AliDielectronPID::kTPC,AliPID::kProton,3.5, 999.);
    pbandwithTOF->AddCut(pidproton);
    cuts->AddCut(pbandwithTOF);
  }
  else if(cutDefinition & kCombinedTpcTof3){
    AliDielectronCutGroup *pbandwithTOF = new AliDielectronCutGroup("pbandwithTOF", "pbandwithTOF", AliDielectronCutGroup::kCompOR);//NOTE taken out for testing

        // EITHER 1. electron candidates inside proton band and with electron TOFpid, rejected if outside 3sigma TOF
    AliDielectronPID *pidTOF = new AliDielectronPID("PIDTOFthere ","TOF nSigma |e|<3.0 + TPC nSigma |p|<3.");
        //TOF bit already internally required
    pidTOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.,0.,100.,kFALSE,AliDielectronPID::kIfAvailable);
  //  pidTOF->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.5,3.5);
    pbandwithTOF->AddCut(pidTOF);

        // OR 2. electron candidates outside proton band accepted if no matching inside with TOF
        // OR 2. electron candidates outside proton band accepted if no matching inside with TOF
    AliDielectronVarCuts *p = new AliDielectronVarCuts("PHigh","P>1.8");
    p->AddCut(AliDielectronVarManager::kPt,1.8, 1e30);
    pbandwithTOF->AddCut(p);
    cuts->AddCut(pbandwithTOF);
  }
  else if(cutDefinition & kCombinedTpcTof4){
    AliDielectronCutGroup *pbandwithTOF = new AliDielectronCutGroup("pbandwithTOF", "pbandwithTOF", AliDielectronCutGroup::kCompOR);//NOTE taken out for testing

        // EITHER 1. electron candidates inside proton band and with electron TOFpid, rejected if outside 3sigma TOF
    AliDielectronPID *pidTOF = new AliDielectronPID("PIDTOFthere ","TOF nSigma |e|<3.0 + TPC nSigma |p|<3.");
        //TOF bit already internally required
    pidTOF->AddCut(AliDielectronPID::kTOF,AliPID::kProton,-999.,-3.0.,100.,kFALSE,AliDielectronPID::kIfAvailable);
   // pidTOF->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.5,3.5);
    pbandwithTOF->AddCut(pidTOF);

        // OR 2. electron candidates outside proton band accepted if no matching inside with TOF
    AliDielectronVarCuts *p = new AliDielectronVarCuts("PHigh","P>1.8");
    p->AddCut(AliDielectronVarManager::kPt,1.8,1e30);
    pbandwithTOF->AddCut(p);
    cuts->AddCut(pbandwithTOF);
  }
  
  
  cuts->AddCut(pid);

  //select conversion electrons selected by the tender
  
//   if(cutDefinition & kV0electrons){
//     AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
//     noconv->SetV0DaughterCut(AliPID::kElectron, kFALSE);
//     cuts->AddCut(noconv);
//   } 
//   else if(cutDefinition & kV0pions){
//     AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
//     noconv->SetV0DaughterCut(AliPID::kPion, kFALSE);
//     cuts->AddCut(noconv);
//   }
  
  
}


/**
*
* Setup the pair cuts according to the given cut definition
*
**/


void SetupPairCuts(AliDielectron *die, ULong_t cutDefinition)
{ 
  
  
  
  // add conversion rejection
  if ( !( cutDefinition & ( kV0electrons | kV0pions) ) ){
    AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
    gammaCut->AddCut(AliDielectronVarManager::kM,0.,.1);
    die->GetPairPreFilter().AddCuts(gammaCut);
    die->SetPreFilterUnlikeOnly(0);
    die->SetPreFilterAllSigns(1);
  }
    
  AliDielectronVarCuts *invMassCut=new AliDielectronVarCuts("InvMass","1.6<M");
  invMassCut->AddCut(AliDielectronVarManager::kM,1.6,5.0);
  die->GetPairFilter().AddCuts(invMassCut);

  AliDielectronVarCuts *y09Cut=new AliDielectronVarCuts("y0.9Cut","y0.9Cut");
  y09Cut->AddCut(AliDielectronVarManager::kY,-.9, .9  );
  die->GetPairFilter().AddCuts(y09Cut); 


  
  
}

//______________________________________________________________________________________
void SetupV0Cuts(AliDielectron *die, ULong_t cutDefinition, ULong_t particle, Bool_t exclude )
{
  //
  // Setup the V0 cuts
  //
  // parameter particle: 0 = electrons, 1 = pions

  
  AliDielectronV0Cuts *v0Cuts = new AliDielectronV0Cuts( Form("IsV0daughter%d", particle), Form("IsV0daughter%d", particle));
 // v0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
//  v0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);//to be checked, if properly filled
//  v0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
//  v0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
  //  v0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle,              0.0,   0.1, kFALSE);
  
  
  

  if(particle == 0){
    v0Cuts->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
    v0Cuts->SetPdgCodes(22,11,11);
    v0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
    v0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
    v0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  }
  else if(particle == 1){
//     v0Cuts->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
     v0Cuts->SetPdgCodes(310,211,211);
//     v0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
//     v0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
//     v0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
     v0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.18,   0.22, kFALSE);
     v0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.5,  0.5, kFALSE); // not sure if it works as expected
     v0Cuts->AddCut(AliDielectronVarManager::kM,                             0.486,   0.508, kFALSE);
  }
    
    
  v0Cuts->SetExcludeTracks( exclude );
  
  //  v0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // not sure if it works as expected
  //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
 //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; &&  const Double_t cutQTG2 < 0.04;

  die->GetTrackFilter().AddCuts(v0Cuts);
}


//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, ULong_t cutDefinition)
{
  AliDielectronHistos *histos = new AliDielectronHistos(die->GetName(), die->GetTitle());

  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");

  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (ULong_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
    if(die->GetMCSignals()) {
      for(ULong_t isig=0; isig<die->GetMCSignals()->GetEntriesFast(); isig++) {
        TString sigMCname = die->GetMCSignals()->At(isig)->GetName();

        // mc truth
        if( !(cutDefinition & kNoPairing) ){
          histos->AddClass(Form("Pair_%s_MCtruth",       sigMCname.Data()));
        }
        histos->AddClass(Form("Track_Legs_%s_MCtruth", sigMCname.Data()));
        // mc reconstructed
        if( !(cutDefinition & kNoPairing) ){
          histos->AddClass(Form("Pair_%s",               sigMCname.Data()));
        }
        histos->AddClass(Form("Track_Legs_%s",         sigMCname.Data()));
      }
    }
  




  //add histograms to event class
  histos->AddClass("Event");
  histos->UserHistogram("Event","","", 4,0.,4.,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","","",300,-15.,15.,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","","", 100,-2.,2.,100,-2.,2.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);
  histos->UserHistogram("Event","","",700,0.0,700.,AliDielectronVarManager::kNTrk);
  
  
  
  histos->UserHistogram("Track","","",200,0.2,20.,AliDielectronVarManager::kPt,kTRUE);
  histos->UserHistogram("Track","","",314,0.0,6.28,AliDielectronVarManager::kPhi); 
  histos->UserHistogram("Track","","",314,0.0,6.28,100,-1.,1.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kEta);  
  
  
  histos->UserHistogram("Track","","",100,0.0,10.,AliDielectronVarManager::kTPCchi2Cl); 
  histos->UserHistogram("Track","","",11,-1.0,10.,AliDielectronVarManager::kITSLayerFirstCls);
  /*
  histos->UserHistogram("Track","","",100,0.0,6.28,50,-1.,1.,100,60.,160.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCActiveLength);
  histos->UserHistogram("Track","","",100,0.0,6.28,50,-1.,1.,100,0.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCchi2Cl); 
  histos->UserHistogram("Track","","",100,0.0,6.28,50,-1.,1.,11,-1.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSLayerFirstCls);  
  histos->UserHistogram("Track","","",100,0.0,6.28,50,-1.,1.,60,60.,120.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);  
  histos->UserHistogram("Track","","",100,0.0,6.28,50,-1.,1.,100,0.,100.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kEta,AliDielectronVarManager::kPt);   
  histos->UserHistogram("Track","","",100,0.0,6.28,50,-1.,1.,100,-4.,4.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);  
  histos->UserHistogram("Track","","", 200,0.1,10.,200,20.,120.,AliDielectronVarManager::kP,AliDielectronVarManager::kTPCsignal,kTRUE);
  
  */
  histos->UserHistogram("Track","","",200,0.1,10.,200,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",200,0.1,10.,200,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Track","","",200,0.1,10.,200,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
  histos->UserHistogram("Track","","",200,0.1,10.,200,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
  histos->UserHistogram("Track","","",200,0.1,10.,200,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
  
  histos->UserHistogram("Track","","",200,0.1,10.,200,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaPio,kTRUE);
  histos->UserHistogram("Track","","",200,0.1,10.,200,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaPro,kTRUE);
  histos->UserHistogram("Track","","",200,0.1,10.,200,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaKao,kTRUE);
  histos->UserHistogram("Track","","",200,-2,2.,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","","",162,-1,161,AliDielectronVarManager::kNclsTPC);          
  histos->UserHistogram("Track","","",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","","",200,-10.,10.,AliDielectronVarManager::kImpactParZ);

  
  histos->UserHistogram("Track","","",314,0.0,6.28,200,0.1,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kP,kFALSE,kTRUE);
  histos->UserHistogram("Track","","",40,-1.0,1.0,200,0.1,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kP,kFALSE,kTRUE);
  histos->UserHistogram("Track","","",100,-20.,20.,200,0.1,10.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kP,kFALSE,kTRUE);
  histos->UserHistogram("Track","","",160,0.,160.,200,0.1,10.,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kP,kFALSE,kTRUE);
  histos->UserHistogram("Track","","",160,0.,160.,200,0.1,10.,AliDielectronVarManager::kTPCsignalN,AliDielectronVarManager::kP,kFALSE,kTRUE);
  histos->UserHistogram("Track","","",160,0.,160.,200,0.1,10.,AliDielectronVarManager::kNFclsTPCr,AliDielectronVarManager::kP,kFALSE,kTRUE);
  histos->UserHistogram("Track","","",101,0.,1.01,200,0.1,10.,AliDielectronVarManager::kNFclsTPCrFrac,AliDielectronVarManager::kP,kFALSE,kTRUE);
  histos->UserHistogram("Track","","",101,0.,1.01,200,0.1,10.,AliDielectronVarManager::kNFclsTPCfCross,AliDielectronVarManager::kP,kFALSE,kTRUE);
  histos->UserHistogram("Track","","",101,0.0,1.01,200,0.1,10.,AliDielectronVarManager::kTPCsignalNfrac,AliDielectronVarManager::kP,kFALSE,kTRUE);
  
  
  
  
  
  histos->UserHistogram("Track","","",314,0.0,6.28,100,20.,120.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",314,0.0,6.28,200,-10.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",314,0.0,6.28,200,-15.,15.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","","",314,0.0,6.28,200,-10.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaPro);
  histos->UserHistogram("Track","","",314,0.0,6.28,200,-10.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTOFnSigmaPro);
  
  histos->UserHistogram("Track","","",314,0.0,6.28,100,-20.,20.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kYsignedIn);
  histos->UserHistogram("Track","","",314,0.0,6.28,160,0.,160.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","","",314,0.0,6.28,160,0.,160.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignalN);
  histos->UserHistogram("Track","","",314,0.0,6.28,160,0.,160.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCr);
  
  histos->UserHistogram("Track","","",314,0.0,6.28,101,0.,1.01,AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCrFrac);
  histos->UserHistogram("Track","","",314,0.0,6.28,101,0.,1.01,AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","","",314,0.0,6.28,101,0.,1.01,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignalNfrac);
  
  
  histos->UserHistogram("Track","","",40,-1.0,1.0,100,20.,120.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",40,-1.0,1.0,100,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",40,-1.0,1.0,100,-15.,15., AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","","",40,-1.0,1.0,100,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPro);
  histos->UserHistogram("Track","","",40,-1.0,1.0,100,-15.,15.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaPro);
  histos->UserHistogram("Track","","",40,-1.0,1.0,100,-20.,20.,AliDielectronVarManager::kEta,AliDielectronVarManager::kYsignedIn);
  histos->UserHistogram("Track","","",40,-1.0,1.0,160,0.,160., AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","","",40,-1.0,1.0,160,0.,160., AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignalN);
  histos->UserHistogram("Track","","",40,-1.0,1.0,160,0.,160., AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCr);
  
  histos->UserHistogram("Track","","",40,-1.0,1.0,101,0.,1.01, AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCrFrac);
  histos->UserHistogram("Track","","",40,-1.0,1.0,101,0.,1.01, AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","","",40,-1.0,1.0,101,0.,1.01, AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignalNfrac);
  
  
  
  
  
  
  
  histos->UserHistogram("Track","","",100,-20.,20.,160,0.,160.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","","",100,-20.,20.,160,0.,160.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCsignalN);
  histos->UserHistogram("Track","","",100,-20.,20.,160,0.,160.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kNFclsTPCr);
  
  histos->UserHistogram("Track","","",100,-20.,20.,101,0.,1.01,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kNFclsTPCrFrac);
  histos->UserHistogram("Track","","",100,-20.,20.,101,0.,1.01,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","","",100,-20.,20.,101,0.,1.01,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCsignalNfrac);
  
  
  
  
  
  //---------------------------------------------------------------
  
  histos->UserHistogram("Track","","",100,-20.,20.,100,20.,120.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",160,0.,160.,100,20.,120.,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",160,0.,160.,100,20.,120.,AliDielectronVarManager::kTPCsignalN,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",160,0.,160.,100,20.,120.,AliDielectronVarManager::kNFclsTPCr,AliDielectronVarManager::kTPCsignal);
  
  histos->UserHistogram("Track","","",101,0.,1.01,100,20.,120.,AliDielectronVarManager::kNFclsTPCrFrac,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",101,0.,1.01,100,20.,120.,AliDielectronVarManager::kNFclsTPCfCross,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",101,0.0,1.01,100,20.,120.,AliDielectronVarManager::kTPCsignalNfrac,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",165,0.,165.,100,20.,120.,AliDielectronVarManager::kTPCActiveLength,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",200,0.3,1.3,100,20.,120.,AliDielectronVarManager::kTPCGeomLength,AliDielectronVarManager::kTPCsignal);
  
  
  //---------------------------------------------------------------
  
  histos->UserHistogram("Track","","",100,-20.,20.,100,-10.,10.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",160,0.,160.,100,-10.,10.,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",160,0.,160.,100,-10.,10.,AliDielectronVarManager::kTPCsignalN,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",160,0.,160.,100,-10.,10.,AliDielectronVarManager::kNFclsTPCr,AliDielectronVarManager::kTPCnSigmaEle);
  
  histos->UserHistogram("Track","","",101,0.0,1.01,100,-10.,10.,AliDielectronVarManager::kNFclsTPCrFrac,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",101,0.0,1.01,100,-10.,10.,AliDielectronVarManager::kNFclsTPCfCross,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",101,0.0,1.01,100,-10.,10.,AliDielectronVarManager::kTPCsignalNfrac,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",165,0.,165.,100,-10.,10.,AliDielectronVarManager::kTPCActiveLength,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",200,0.3,1.3,100,-10.,10.,AliDielectronVarManager::kTPCGeomLength,AliDielectronVarManager::kTPCnSigmaEle);
  
  
  //---------------------------------------------------------------
  
  histos->UserHistogram("Track","","",100,-20.,20.,100,-15.,15.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","","",160,0.,160.,100,-15.,15.,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","","",160,0.,160.,100,-15.,15.,AliDielectronVarManager::kTPCsignalN,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","","",160,0.,160.,100,-15.,15.,AliDielectronVarManager::kNFclsTPCr,AliDielectronVarManager::kTPCnSigmaPio);
  
  histos->UserHistogram("Track","","",101,0.0,1.01,100,-15.,15.,AliDielectronVarManager::kNFclsTPCrFrac,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","","",101,0.0,1.01,100,-15.,15.,AliDielectronVarManager::kNFclsTPCfCross,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","","",101,0.0,1.01,100,-15.,15.,AliDielectronVarManager::kTPCsignalNfrac,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","","",165,0.,165.,100,-15.,15.,AliDielectronVarManager::kTPCActiveLength,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","","",200,0.3,1.3,100,-15.,15.,AliDielectronVarManager::kTPCGeomLength,AliDielectronVarManager::kTPCnSigmaPio);
  //---------------------------------------------------------------
  
  histos->UserHistogram("Track","","",100,-20.,20.,100,-15.,15.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTOFnSigmaPro);
  histos->UserHistogram("Track","","",160,0.,160.,100,-15.,15.,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kTOFnSigmaPro);
  histos->UserHistogram("Track","","",160,0.,160.,100,-15.,15.,AliDielectronVarManager::kTPCsignalN,AliDielectronVarManager::kTOFnSigmaPro);
  histos->UserHistogram("Track","","",160,0.,160.,100,-15.,15.,AliDielectronVarManager::kNFclsTPCr,AliDielectronVarManager::kTOFnSigmaPro);
  
  histos->UserHistogram("Track","","",101,0.0,1.01,100,-15.,15.,AliDielectronVarManager::kNFclsTPCrFrac,AliDielectronVarManager::kTOFnSigmaPro);
  histos->UserHistogram("Track","","",101,0.0,1.01,100,-15.,15.,AliDielectronVarManager::kNFclsTPCfCross,AliDielectronVarManager::kTOFnSigmaPro);
  histos->UserHistogram("Track","","",101,0.0,1.01,100,-15.,15.,AliDielectronVarManager::kTPCsignalNfrac,AliDielectronVarManager::kTOFnSigmaPro);
  histos->UserHistogram("Track","","",165,0.,165.,100,-15.,15.,AliDielectronVarManager::kTPCActiveLength,AliDielectronVarManager::kTOFnSigmaPro);
  histos->UserHistogram("Track","","",200,0.3,1.3,100,-15.,15.,AliDielectronVarManager::kTPCGeomLength,AliDielectronVarManager::kTOFnSigmaPro);
  
  
  
  
  
  
  histos->UserHistogram("Track","","", 165,0.,165.,AliDielectronVarManager::kTPCActiveLength);
  histos->UserHistogram("Track","","", 300,-1.,6.5,165,0.,165.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCActiveLength);
  histos->UserHistogram("Track","","",200,-20.,20.,165,0.,165.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCActiveLength);
  
  
  histos->UserHistogram("Track","","",200,0.3,1.3,AliDielectronVarManager::kTPCGeomLength);
  histos->UserHistogram("Track","","",100,-20.,20.,200,0.3,1.3,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCGeomLength);
  
  
  histos->UserHistogram("Track","","",165, 0.,165.,200, 0.3,1.3,AliDielectronVarManager::kTPCActiveLength,AliDielectronVarManager::kTPCGeomLength);
  
  
  if( !(cutDefinition & kNoPairing) ){
  
        //Pair classes
      // to fill also mixed event histograms loop until 10
    ULong_t pairClasses [5] = { AliDielectron::kEv1PP,  AliDielectron::kEv1PM, AliDielectron::kEv1MM, AliDielectron::kEv1MEv2P, AliDielectron::kEv1PMRot};
      
      
    for (ULong_t i=0; i<5; ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(   pairClasses[i]   ) ) );
    }

      //add MC signal histograms to track and pair class
  
    
    //add histograms to Pair classes
    histos->UserHistogram("Pair","","",300,0.,6.,AliDielectronVarManager::kM);
    histos->UserHistogram("Pair","","",300,0., 6.,100,0.,10., AliDielectronVarManager::kM,AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","","",100,0., 6.,200,-2.,2., AliDielectronVarManager::kM,AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","","",200,0,20.,AliDielectronVarManager::kPt,kTRUE);
    histos->UserHistogram("Pair","","",200,-2.,2.,AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","","",315,0.,3.15,AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","","", 200,-2.,2.,AliDielectronVarManager::kEta);
    histos->UserHistogram("Pair","","",100,0.,6.,315,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);

  }
  

  
  
  // add histograms for track PID performance study  
/*  

  const ULong_t dim = 4;
  ULong_t bins1[dim] = {100, 100, 90, 100};
  Double_t mins1[dim] = {0., -1., 70., 20.};
  Double_t maxs1[dim] = {10., 1., 160., 120.};
  ULong_t vars1[dim] = {AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kNclsTPC, AliDielectronVarManager::kTPCsignal  };
  histos->UserSparse("Track", dim, bins1, mins1, maxs1, vars1);
  
  
  ULong_t bins2[dim] = {100, 100, 90, 100};
  Double_t mins2[dim] = {0, -1., 70., -4.};
  Double_t maxs2[dim] = {10, 1., 160., 10.};
  ULong_t vars2[dim] = {AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kNclsTPC, AliDielectronVarManager::kTPCnSigmaPio};
  histos->UserSparse("Track", dim, bins2, mins2, maxs2, vars2);
  

  ULong_t bins3[dim] = {100, 100, 90, 100};
  Double_t mins3[dim] = {0, -1., 70., -10.};
  Double_t maxs3[dim] = {10, 1., 160., 4.};
  ULong_t vars3[dim] = {AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kNclsTPC, AliDielectronVarManager::kTPCnSigmaEle};
  histos->UserSparse("Track", dim, bins3, mins3, maxs3, vars3);*/


  die->SetHistogramManager(histos);
  
}

void InitCF(AliDielectron* die, ULong_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 1.0, 1.3, 2.0, 3.0, 5., 7.0, 10.0, 100.0");
  
  cf->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
  cf->AddVariable(AliDielectronVarManager::kY,"-1,-0.9,-0.8,-0.3,0,0.3,0.9,1.0");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);

  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 1.0, 1.1, 1.2, 1.3, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-1.,-0.9,-0.8,0,0.8,0.9,1.0",kTRUE);
  
  cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,0,6,kTRUE);

  
  if (hasMC){
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,10000,-5000.5,4999.5,kTRUE);

    //only steps for efficiencies
    cf->SetStepsForMCtruthOnly();
  }
  
  //only in this case write MC truth info
  if (cutDefinition==0){
    cf->SetStepForMCtruth();
  }
  
  cf->SetStepsForSignal();
  
  die->SetCFManagerPair(cf);
}

void SetupMCsignals(AliDielectron *die){

  
  AliDielectronSignalMC* inclusiveJpsi = new AliDielectronSignalMC("inclusiveJpsi","Inclusive J/psi");
  inclusiveJpsi->SetLegPDGs(11,-11);
  inclusiveJpsi->SetMotherPDGs(443,443);
  inclusiveJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  inclusiveJpsi->SetFillPureMCStep(kTRUE);
  inclusiveJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  inclusiveJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(inclusiveJpsi);
  
  
  
  AliDielectronSignalMC* promptJpsi = new AliDielectronSignalMC("promptJpsi","Prompt J/psi");   // prompt J/psi (not from beauty decays)
  promptJpsi->SetLegPDGs(11,-11);
  promptJpsi->SetMotherPDGs(443,443);
  promptJpsi->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsi->SetFillPureMCStep(kTRUE);
  promptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(promptJpsi);
  
  
  AliDielectronSignalMC* nonpromptJpsi = new AliDielectronSignalMC("nonpromptJpsi","Nonprompt J/psi");   // nonprompt J/psi (from beauty decays)
  nonpromptJpsi->SetLegPDGs(11,-11);
  nonpromptJpsi->SetMotherPDGs(443,443);
  nonpromptJpsi->SetGrandMotherPDGs(503,503);   // from beauty hadrons
  nonpromptJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  nonpromptJpsi->SetFillPureMCStep(kTRUE);
  nonpromptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  nonpromptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  nonpromptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(nonpromptJpsi);
  
}


void SetEtaCorrectionGSI()
{

}

void SetEtaCorrection()
{
  
  if (AliDielectronPID::GetEtaCorrFunction()) return;
  
    
  TString list=gSystem->Getenv("LIST");
  if( list.IsNull() ) list = gSystem->Getenv("PERIOD_NAME");
  if( list.IsNull() ){
      printf("Cannot get run list for eta correction!");
      return;
  }
  Bool_t correctionFound = kFALSE;
  TF1*func = new TF1("etaCorrection","pol4",0,1);
  if( list.Contains("LHC10c") ){
    func->SetParameters(0.980352-0.00906396,0.0596654,0.0111099,0.0286887);
    correctionFound = kTRUE;
  }
  if( list.Contains("LHC10d") ){
    func->SetParameters(0.979027,0.00365582,0.0590777,-0.00462389,0.0195373);
    correctionFound = kTRUE;
  }
  if( list.Contains("LHC10e") ){
    func->SetParameters(0.976058,0.00317814,0.0686704,-0.00372511,0.0154494);
    correctionFound = kTRUE;
  }
  if( correctionFound ){
    printf("Using local Eta Correction: ");
    func->Print();
    AliDielectronPID::SetEtaCorrFunction(func);
  }
  
  if(!correctionFound){
      printf("No Eta correction function found."):
  }
  
}