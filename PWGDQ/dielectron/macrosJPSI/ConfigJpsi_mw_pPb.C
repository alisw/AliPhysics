void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);


TString namesDieleData=("basicQ+SPDfirst+pt>1+PID; basicQ+SPDany+pt>1+PID; basicQ+ITS012+pt>1+PID; basicQ+SPDany+pt>1+noexclPIDforcontrolpurpose");//; basicQ+SPDany+pt>1+p>1.2+PIDrequirementsHFEnoexclPIDforcontrolpurpose; basicQ+SPDany+pt>1+p>1.2+PIDrequirementsHFE+PID; basicQ+SPDany+p>1.+PID");
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
  if(cutDefinition == 2 || cutDefinition == 3 || cutDefinition == 4 || cutDefinition == 5|| cutDefinition == 6) InitCFDieleData(diele, cutDefinition, isAOD);//first 2 cut sets in 3rd included

  AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
  if(cutDefinition == 0 || cutDefinition == 1 || cutDefinition == 2 || cutDefinition == 3 || cutDefinition == 4 || cutDefinition == 5|| cutDefinition == 6) rot->SetConeAnglePhi(TMath::Pi());
  //else if(cutDefinition == 2) rot->SetConeAnglePhi(TMath::Pi()/180.*135.);
  //rot->SetIterations(10);
  rot->SetIterations(20);
  diele->SetTrackRotator(rot);

  /*  AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
  mix->SetMixType(AliDielectronMixingHandler::kAll);
  mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
  //NOTE: other mixing classes probably needed!!!
  mix->SetDepth(10);
  diele->SetMixingHandler(mix);*/
 
  
  return diele;
}

//______________________________________________________________________________________
void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setup the track cuts
  //
  
  //if (!isAOD) {//NOTE: seems to work, see AliDielectronTrackCuts method IsSelected at the beginning, to be checked with AODs
    //exclude conversion electrons selected by the tender
    AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
    noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
    diele->GetTrackFilter().AddCuts(noconv);
    // }
    AliDielectronTrackCuts *trackCuts=new AliDielectronTrackCuts("ITSandgeneral_trackCuts","ITSandgeneral_trackCuts");
    //ITS related cuts
    if (cutDefinition==0)
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    else if ((cutDefinition==1)||(cutDefinition ==3)||(cutDefinition ==4)||(cutDefinition ==5)||(cutDefinition ==6))
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
      trackCuts->SetRequireTPCRefit(kTRUE);
      trackCuts->SetRequireITSRefit(kTRUE);
      diele->GetTrackFilter().AddCuts(trackCuts);
  
  //Pt cut ----------------------------------------------------------
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("trackkineandTPCQ","trackkine_and_TPCQ");
  if ((cutDefinition==0)||(cutDefinition==1)||(cutDefinition ==3)||(cutDefinition ==4)||(cutDefinition ==5))
    pt->AddCut(AliDielectronVarManager::kPt,1.,1e30);
  else if (cutDefinition==6) pt->AddCut(AliDielectronVarManager::kP,1.,1e30);
  pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
  // ITS 0 1 2 : only for third variant
  if (cutDefinition==2) pt->AddCut(AliDielectronVarManager::kITSLayerFirstCls,0.,2.5);
  //AOD additions since there are no AliESDtrackCuts -----------------
  //
  // TPC #clusteres cut
  pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
  if((cutDefinition==4)||(cutDefinition==5))  pt->AddCut(AliDielectronVarManager::kNclsTPC,120.,160.);
  pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
  pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
  //TODO: DCA cuts to be investigated!!! NOTE: why?? (mwinn, 15.01.2013)
  pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
  pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  
  diele->GetTrackFilter().AddCuts(pt);
    
  // PID cuts --------------------------------------------------------
  if(cutDefinition ==0 || cutDefinition ==1 ||cutDefinition ==2||cutDefinition ==6){
    AliDielectronPID *pid = new AliDielectronPID("PID10","TPC nSigma |e|<3 + |Pi|>3.5 + P>3");
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.5,3.5,0.,0.,kTRUE);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-20.,3.,0.,0.,kTRUE);
    diele->GetTrackFilter().AddCuts(pid);
  }
  if(cutDefinition ==3){
    AliDielectronPID *pid = new AliDielectronPID("PID10","TPC nSigma |e|<3");
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
    diele->GetTrackFilter().AddCuts(pid);
    /* //taken out at the 04.02.2013 in order to have some feeling for PID cuts in QAplots...
    AliDielectronVarCuts *pidsubs = new AliDielectronVarCuts("pidSubs","pidsubs cut");
     pidsubs->AddCut(AliDielectronVarManager::kP,1.2,1e30);
     pidsubs->AddCut(AliDielectronVarManager::kTPCsignal,70.,110.);
     diele->GetTrackFilter().AddCuts(pidsubs);*/
  }
  if(cutDefinition ==4||cutDefinition ==5){
    
    AliDielectronPID *pid = new AliDielectronPID("PID10","TPC nSigma |e|<3");
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
    diele->GetTrackFilter().AddCuts(pid);
    AliDielectronVarCuts *hfe = new AliDielectronVarCuts("hfeTPCQ","hfeTPCQ");
    hfe->AddCut(AliDielectronVarManager::kTPCsignalN,80.,160.);//TPC PID Cluster
    diele->GetTrackFilter().AddCuts(hfe);
    if(cutDefinition ==5){
      AliDielectronPID *pidexclusion = new AliDielectronPID("PIDhfe","TPC nSigma |Pi|>3.5 + P>3");
      pidexclusion->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.5,3.5,0.,0.,kTRUE);
      pidexclusion->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-20.,3.,0.,0.,kTRUE);
      diele->GetTrackFilter().AddCuts(pidexclusion);
    }
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
  for (Int_t i=0; i<3/*for mixing until 10*/; ++i){
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
  histos->UserHistogram("Track","TOFPIDBit","TOFPIDBit;TOFPIDBit;#tracks",2,-0.5,1.5,AliDielectronVarManager::kTOFPIDBit);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","TPCsignalN","Number of PID Clusters TPC;TPC PID number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  histos->UserHistogram("Track","nClsoverfindablecluster","Number of found Clusters TPC over findably ;TPC number cluster over findable;#tracks",160,0.0,1.1,AliDielectronVarManager::kNFclsTPCrFrac);
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                        100,-1,1,144,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        200,0.2,20.,100,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaPi_P","TPC number of sigmas Kaons;PIN [GeV];TPC number of sigmas Pions;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaPro_P","TPC number of sigmas Protons;PIN [GeV];TPC number of sigmas Protons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;PIN [GeV];TPC number of sigmas Kaons;#tracks",
                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
  histos->UserHistogram("Track","TOFbeta_P","TOF beta;P [GeV];TOF beta;#tracks",
                        200,0.2,20.,100,0.,1.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaEle_P","TOF number of sigmas Electrons;P [GeV];TOF number of sigmas Electrons;#tracks",
                        200,0.2,20.,200,-20.,20.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TRDnCls","Number of Clusters TRD;TRD number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTRD);
  histos->UserHistogram("Track","TRDntracklets","Number of tracklets TRD;TRD number tracklets;#tracks",7,-0.5,6.5,AliDielectronVarManager::kTRDntracklets);
  histos->UserHistogram("Track","TRDprobEle_P","TRD electron prob.;P [GeV];TRD electron prob.;#tracks",
                        200,0.2,20.,100,.0,1.,AliDielectronVarManager::kP,AliDielectronVarManager::kTRDprobEle,kTRUE);
  histos->UserHistogram("Track","TRDprobEle2D_P","TRD electron prob. 2D;P [GeV];TRD electron prob. 2D;#tracks",
                        200,0.2,20.,100,.0,1.,AliDielectronVarManager::kP,AliDielectronVarManager::kTRDprob2DEle,kTRUE);

  


      
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        125,0.,125*.04,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","Pt","Pt;Pt;#pairs",
                        200,0.,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","OpeningAngletransverse","Opening angle transverse;angle",
                        100,0.,3.15,AliDielectronVarManager::kDeltaPhi);
  histos->UserHistogram("Pair","Chi2NDF","chisquareNDF;chisquare/ndf;#pairs",
                        100,0.,30.,AliDielectronVarManager::kChi2NDF);
  histos->UserHistogram("Pair","distanceXY","distancelegsXY;distanceXY[cm];#pairs",
			100,0.,.0001,AliDielectronVarManager::kLegDistXY);
  histos->UserHistogram("Pair","distance","distancelegs;distance[cm];#pairs",
                        100,0.,.0001,AliDielectronVarManager::kLegDist);
  histos->UserHistogram("Pair","pseudoproperdecaylength","pseudoproperdecaylength;pseudoproperdecaylength[cm];#pairs",
			100,0.,.5,AliDielectronVarManager::kPseudoProperTime);
  histos->UserHistogram("Pair","Armenteros-Podolanski","Armenteros-Podolanski;ArmAlpha;ArmPt[GeV];#tracks",
			100,-10.0,10.,100,0.,2.,AliDielectronVarManager::kArmAlpha,AliDielectronVarManager::kArmPt,kTRUE);
 
 

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
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.5, 3.8, 4.2, 4.6, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0");
  cf->AddVariable(AliDielectronVarManager::kY,"-5,-1,-0.9,-0.8,-0.7,-0.5,-0.3,0.3,0.5,0.7,0.8,0.9,1.0,5");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPseudoProperTime,150,-0.3,0.3);
  cf->AddVariable(AliDielectronVarManager::kPseudoProperTimeErr,600,0.,0.3);
  cf->AddVariable(AliDielectronVarManager::kPairType,"-0.5,0.5,1.5,2.5,9.5,10.5");
  cf->AddVariable(AliDielectronVarManager::kArmAlpha,"-10.,-5.,-2.,-.1.,0.,1.,2.,5.,10.");
  cf->AddVariable(AliDielectronVarManager::kArmPt,"0.0,0.4,0.6,0.8,1.0,1.5,2.0,3.0");
  

  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 70, 75, 80, 85, 90, 100, 120, 160",kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kTPCsignalN,"0, 50, 60, 70, 80, 90, 100, 160",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCchi2Cl,"0, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.",kTRUE);  
  cf->AddVariable(AliDielectronVarManager::kEta, "-1.0, -0.88, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.88, 1.0", kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-2.5,-2,-1.5,-1,-0.5,4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.,3.5,4.,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.,3.5,4.,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,-0.5,5.5,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kP,"0.0, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2,1.5,2.0,3.0,4.0,5.0, 100.0",kTRUE);
  //event variables
  //cf->AddVariable(AliDielectronVarManager::kNaccTrcklts,"0.0, 9.0, 17.0, 25.0, 36.0, 55.0, 500.0");
  cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd10,101,-0.5,100.5);
  //cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd10Corr,101,-0.5,100.5);
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

