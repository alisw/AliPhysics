namespace ConfDef {

void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

void SetEtaCorrection(AliDielectron *die);

  enum {kTPCDefault,kTPCStrict,kTPCTOFlow,kTPCTOF,kQA};

  // enum {kTPCDefault,kQA}; 
}

TString names=("TPCDefault;TPCStrict;kTPCTOFlow;kTPCTOF;QA");
//TString names=("TPCDefault;QA");

TObjArray *arrNames=names.Tokenize(";");

const Int_t nDie=arrNames->GetEntries(); // loop cutDefinition < nDie in Addtask

void ConfigDefault(AliAnalysisTaskMultiDielectron *task)
{
  //
  // Setup the instance of AliDielectron
  //
  
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //Do we have an AOD handler?
  Bool_t isAOD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliAODInputHandler::Class() ? kTRUE : kFALSE);

  //===================================================================
  // EVENT CUTS
  //==================================================================
  //set trigger
  task->SetTriggerMask(AliVEvent::kINT7);

  task->UsePhysicsSelection();

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  //  eventCuts->SetRequireCompleteDAQ(kTRUE);
  task->SetEventFilter(eventCuts);

  //======================================================================


  //Loop over all cutDefinitions
  for (Int_t cutDefinition=0; cutDefinition<nDie; ++cutDefinition){
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()),Form("Track cuts: %s",name.Data()));
  //  ConfDef::SetEtaCorrection(die);

     
  die->SetUseKF(kFALSE);
  
  if (hasMC) ConfDef::SetupMCsignals(die);

  //PID QA
  //  ConfDef::SetupV0Cuts(die,cutDefinition);

  // cut setup
  ConfDef::SetupTrackCuts(die,cutDefinition);  
  //
  // ConfDef::SetupPairCuts(die,cutDefinition);
  //
  //PID QA
  ConfDef::SetupV0Cuts(die,cutDefinition);



  ConfDef::InitHistograms(die,cutDefinition);
  //
  if (hasMC) ConfDef::InitCF(die,cutDefinition);
  //
  //  if (cutDefinition==ConfDef::kQA) die->SetNoPairing();
  //
  die->SetNoPairing();

  /*
AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
  mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-8,-5,0,5,8,10");
//   mix->AddVariable(AliDielectronVarManager::kPhi ,"0,3.14,6.3");
  mix->SetDepth(10);
 die->SetMixingHandler(mix);
  */
  //
  // ConfDef::SetEtaCorrection();
  task->AddDielectron(die);


  }//end loop

}

namespace ConfDef{
//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);

  
  //default quality cuts
    AliDielectronTrackCuts *cut1=new AliDielectronTrackCuts("cut1","cut1");
    //cut1->SetRequireITSRefit(kTRUE); //ITS refit for conversion elecs is probably a too hard cut, danger: pile up, also: TPC-TRD matching eff for tpc only tracks maybe worse??
  cut1->SetRequireTPCRefit(kTRUE);
  //cut1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); // no ITS requirement, danger: pile up
  cuts->AddCut(cut1);
  
   //pt and kink mother
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");
  pt->AddCut(AliDielectronVarManager::kP,1.,1.e30); //conversion legs reconstructed above 50MeV
  //inefficiency induced by pT cut: now we only accept conversions in which both legs are above the pT threshold,
  //but conversions are likely to be kinematically asymmetric. better: use no threshold here and use the cfcontainer for studies afterwards
  //  pt->AddCut(AliDielectronVarManager::kP,1.5,1.e30); //for conversion elecs for TPC PiD checks: p (total momentum)> 1.5GeV (avoid proton contamination)
  pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.); //no daughter cut?

  
  //DCA cuts

  //IMHO(DW): makes sense for e+-??, only efficiency loss, no purity gain, what else?
  // new IMHO: may make sense to make conversions more primary-like
  
  pt->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  pt->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
  
  if (cutDefinition==kTPCTOFlow){
  AliDielectronPID *pid = new AliDielectronPID("PID","PID");
  pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.,0.,1.5);
  cuts->AddCut(pid);
  }
  if (cutDefinition==kTPCTOF){

    AliDielectronPID *pid = new AliDielectronPID("PID","PID");
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.,0.,1.e10);
    cuts->AddCut(pid);
  }
  
  if (cutDefinition==kTPCDefault){

     pt->AddCut(AliDielectronVarManager::kP,1.5,1.e30);
  }
  
  
  if (cutDefinition==kTPCStrict){

    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-2.,3.);

    
  }
  
  
  cuts->AddCut(pt);

  AliDielectronVarCuts *an = new AliDielectronVarCuts("ancuts","ancuts");

  an->AddCut(AliDielectronVarManager::kEta,-0.9,0.9); //p-Pb cuts
  an->AddCut(AliDielectronVarManager::kNclsTPC,0,70,kTRUE); // exclude this range
  //TPC XI square cut???

  
  cuts->AddCut(an);
  
  /*
  //exclude conversion electrons selected by the tender
  AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  //noconv->SetV0DaughterCut(AliPID::kElectron,kFALSE); //include conversion electrons?
  cuts->AddCut(noconv);
  */
}

//______________________________________________________________________________________
  void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition)
  {
    //
    // Setup the V0 cuts
    //

    // Quality cuts (add the gamma filter to the cut group)
    TIter next(die->GetTrackFilter().GetCuts());
    AliAnalysisCuts *cuts;
    while((cuts = (AliAnalysisCuts*)next())) {
      if(cuts->IsA() == AliDielectronCutGroup::Class())  break;
    }

    AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");
    gammaV0Cuts->SetPdgCodes(22,11,11);


    // add default PID cuts (defined in AliDielectronPID)
    // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)

    gammaV0Cuts->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
    //gammaV0Cuts->SetDefaultPID(16, AliDielectronV0Cuts::kAny); // TPC+- 3.5, TOF +-3 if available 
    
    //gammaV0Cuts->SetDefaultPID(16); // TPC+- 3.5, TOF +-3 if available


    //   gammaV0Cuts->SetDefaultPID(7); // TPC+-10, TOF+-3 required in 0<p<1.5
    //gammaV0Cuts->SetDefaultPID(14); // TPC+-10, TOF+-3 if available, TRD 90% eleEff if available
    //gammaV0Cuts->SetDefaultPID(15); // TPC+-10, TOF+-3 if availbale, TRD 90% eleEff if available, TRD Chi2<2 
    //gammaV0Cuts->SetDefaultPID(8); // TOF+-5 if available
    
    gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
    //  gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle,              0.0,   0.1, kFALSE);
    // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.03, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kPhivPair,                         2.8,   3.2, kFALSE);
    //  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                      -0.35,  0.35, kFALSE); //does weird stuff
    //could also have a look at ArmAlpha, a cut on this should in principle improve the purity, however it is not clear if the cut works (-> see JB)
    //INCLUDE CONVERSIONS!!!!!!
    gammaV0Cuts->SetExcludeTracks(kFALSE); //exclude oder include conversion e+-
    // gammaV0Cuts->Print();

    //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
    //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; &&  const Double_t cutQTG2 < 0.04;

    if(cuts)
      ((AliDielectronCutGroup*)cuts)->AddCut(gammaV0Cuts);
    else
      die->GetTrackFilter().AddCuts(gammaV0Cuts);
  }
  
  
//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //

  /*
  // add conversion rejection
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,.05);
  die->GetPairPreFilter().AddCuts(gammaCut);
  die->SetPreFilterUnlikeOnly();
  */
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
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
  }

  //add histograms to event class
  if (cutDefinition==0) {
    histos->AddClass("Event");
    histos->UserHistogram("Event","","",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","","",
	100,-2.,2.,100,-2.,2.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);
    histos->UserHistogram("Event","","",100,0.,100,AliDielectronVarManager::kCentralityNew);

    histos->UserHistogram("Event","","",600,0.,6000.,AliDielectronVarManager::kRefMultTPConly);
    
    histos->UserHistogram("Event","","", 500,0.,10000.,AliDielectronVarManager::kNacc);

    histos->UserHistogram("Event","","", 500,0.,40000.,
			  AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","","",
			  10,0.,100., 500,0.,40000.,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kMultV0);    

    histos->UserHistogram("Event","","",
			  200,0.,10000., 200,0.,40000.,AliDielectronVarManager::kNacc,AliDielectronVarManager::kMultV0);
    
  }

  //add histograms to Track classes
  histos->UserHistogram("Track","","",200,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","","",
      400,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","","",
      100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",
      100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Track","","",
      100,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);

  histos->UserHistogram("Track","","",
      100,-2,2,144,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","","",
      160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","","",
      100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","","",
      150,-15,15,160,-0.5,159.5,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCsignalN);

  histos->UserHistogram("Track","","",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","","",600,-3.,3.,AliDielectronVarManager::kImpactParZ);

  // new histos 
  //histos->UserHistogram("Track","nameInRootF","hist title (on canvas); X-Axis-Label; Y-Axis-Label",
  histos->UserHistogram("Track","","",60,-4.5,1.5,AliDielectronVarManager::kKinkIndex0);
  //
  // 1dim histos:
  //

  histos->UserHistogram("Track","","",7,-0.5,6.5.,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","","",50,0.,25.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","","",7,-0.5,6.5.,AliDielectronVarManager::kITSLayerFirstCls);
  
  histos->UserHistogram("Track","","",161,-0.5,160.5.,AliDielectronVarManager::kNclsTRD); //how many clusters are possible?
  histos->UserHistogram("Track","","",7,-0.5,6.5,AliDielectronVarManager::kTRDntracklets);
  histos->UserHistogram("Track","","",7,-0.5,6.5,AliDielectronVarManager::kTRDpidQuality);
  histos->UserHistogram("Track","","",61,-0.5,60.5,AliDielectronVarManager::kTRDchi2);
  histos->UserHistogram("Track","","",50,0.,1.,AliDielectronVarManager::kTRDprobEle);
  histos->UserHistogram("Track","","",50,0.,1.,AliDielectronVarManager::kTRDprobPio);
  histos->UserHistogram("Track","","",50,0.,1.,AliDielectronVarManager::kTRDprob2DEle);
  histos->UserHistogram("Track","","",50,0.,1.,AliDielectronVarManager::kTRDprob2DPio);
  histos->UserHistogram("Track","","",100,0.,100.,AliDielectronPID::kTRDeleEff);
  histos->UserHistogram("Track","","",100,0.,100.,AliDielectronPID::kTRDeleEff2D);
  histos->UserHistogram("Track","","",300,-7.,7.,AliDielectronVarManager::kTRDphi);
  histos->UserHistogram("Track","","",50,0.,1.,AliDielectronVarManager::kTRDpidEffLeg);
  histos->UserHistogram("Track","","",201,-0.5,200.5,AliDielectronVarManager::kTRDsignal);
  //TOF histos

  histos->UserHistogram("Track","","",
			250,0.0,5.0,300,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","","",
			100,0.2,20.,50,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);


  //
  // 2dim histos:
  //
  histos->UserHistogram("Track","","",400,0.2,20.,6,-4.5,1.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kKinkIndex0);
  //chi2
  histos->UserHistogram("Track","","",400,0.2,20.,61,0.,60.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTRDchi2);
  histos->UserHistogram("Track","","",500,-2.,2.,61,-0.5,60.5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTRDchi2);
  histos->UserHistogram("Track","","",350,0.,7.,61,-0.5,60.5,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTRDchi2);
  //Nclusters
  histos->UserHistogram("Track","","",350,0.,7.,101,-0.5,100.5,AliDielectronVarManager::kPhi,AliDielectronVarManager::kNclsTRD);
  histos->UserHistogram("Track","","",400,0.2,20.,101,-0.5,100.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kNclsTRD);
  histos->UserHistogram("Track","","",300,-2.,2.,101,-0.5,100.5,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsTRD);
  //various other 2dhis
  histos->UserHistogram("Track","","",400,0.2,20.,400,0.2,20.,AliDielectronVarManager::kPt,AliDielectronVarManager::kPOut);

  histos->UserHistogram("Track","","",200,-1,1.,200,-10,10,AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",100,0,100.,200,-10,10,AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",200,0,6.3.,200,-10,10,AliDielectronVarManager::kPhi, AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",200,0,200.,200,-10,10,AliDielectronVarManager::kNclsTPC, AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",300,-15.,15.,200,-10,10,AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",200,0.,10000.,200,-10,10,AliDielectronVarManager::kNacc, AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",200,0.,40000.,200,-10,10,AliDielectronVarManager::kMultV0, AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",10,0,100., 18,-0.9,0.9 ,100,-5,5,AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  

  

  
  
  //add histograms to Pair classes
  histos->UserHistogram("Pair","","",
      301,-.01,6.01,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","","",
      100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","","",
      100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","","",100,-1,1.,100,0.,1.,AliDielectronVarManager::kArmAlpha,AliDielectronVarManager::kArmPt);
  histos->UserHistogram("Pair","","",
			100,0.,3.15,AliDielectronVarManager::kPhivPair);

  histos->UserHistogram("Pair","","",100,0.,3.15,100, 0, 6,AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kM);


  
  die->SetHistogramManager(histos);

}

void InitCF(AliDielectron* die, Int_t cutDefinition)
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

  cf->AddVariable(AliDielectronVarManager::kP,"0.0, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 2.0, 2.1, 2.2, 2.5, 2.6, 2.7, 3.0, 3.1, 3.2, 3.5, 3.6, 3.7, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-1.,-0.9,-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kCentralityNew,20,0,100);
  //cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,40,-20*0.2.,20*.2, kTRUE); //0.2Sigma Steps
  
  //cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,0,6,kTRUE);


  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
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

}


void SetEtaCorrection(AliDielectron *die)
{
  //Apply post calibration here
  /*
  if (AliDielectronPID::GetCentroidCorrFunction()) return;
  //if (AliDielectronPID::GetEtaCorrFunction()) return;
  //TString list=gSystem->Getenv("LIST");


  // TH2F *CentroidCorrMap =0x0;
  //TH2F *WidthCorrMap = 0x0;
  
  //try how good unsmoothed/unfitted correction maps already are by using the histos
  TFile f("$TRAIN_ROOT/dweiser_ff/jpsiPbPb/PostCalibMaps.root");
  if (!f.IsOpen()){
    printf("WARNING!!! Could not open TPC post calibration file!!! No postcalibration available! WARNING!!!");
    return;
  }

  //get pointers to correction histos
  TH2F *CentroidCorrMap = f.Get("CorrMapCentroid");
  TH2F *WidthCorrMap = f.Get("CorrMapWidth");

  if (!CentroidCorrMap || !WidthCorrMap) printf("WARNING!!! At least one correction map could not be opened!!!!");
  
  //set histos as corr histo
  if (CentroidCorrMap!=NULL && WidthCorrMap!=NULL){
  die->SetCentroidCorrFunction(CentroidCorrMap,AliDielectronVarManager::kCentrality, AliDielectronVarManager::kEta );
  die->SetWidthCorrFunction(WidthCorrMap,AliDielectronVarManager::kCentrality, AliDielectronVarManager::kEta );
  }

  
  */
}
}
