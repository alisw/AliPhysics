void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition);
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition);
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition);
void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition);
void SetupV0Cuts(AliDielectron *diele, Int_t cutDefinition);
void SetSignals(AliDielectron *diele);


TString namesDieleData=("basicQ+SPDfirst+pt>1+PID;basicQ+SPDany+pt>1+PID;basicQ+SPDany+pt>1+noexclPIDforcontrolpurpose;nocuts");
//has to introduce PID cuts....
//TString namesDieleData=("basicQ+SPDfirst+pt>1+PID");

TObjArray *arrNamesDieleData=namesDieleData.Tokenize(";");

const Int_t nDie=arrNamesDieleData->GetEntries();

AliDielectron* ConfigJpsi_mw_pPb_MC(Int_t cutDefinition) 
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
  diele->SetEstimatorFilename("$ALICE_PHYSICS/PWGDQ/dielectron/files/estimators.root");
  //diele->SetEstimatorFilename("estimators.root");
  // cut setup
  SetupTrackCutsDieleData(diele, cutDefinition);
  SetupPairCutsDieleData(diele, cutDefinition);
  SetupV0Cuts(diele, cutDefinition);	

   // Set MC signals
  SetSignals(diele);
  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielelectron framework histograms will be filled
  //
  InitHistogramsDieleData(diele, cutDefinition);
  
  // the last definition uses no cuts and only the QA histograms should be filled!, now for all cuts
  if(cutDefinition == 1 ||cutDefinition == 2 || cutDefinition == 3 || cutDefinition == 4) InitCFDieleData(diele, cutDefinition);//first 2 cut sets in 3rd included
 
  
  return diele;
}

//______________________________________________________________________________________
void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //
  
  //NOTE: seems to work, see AliDielectronTrackCuts method IsSelected at the beginning, to be checked with AODs
    //exclude conversion electrons selected by the tender
 /* if(!cutDefinition==4){   
    AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
    noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
    diele->GetTrackFilter().AddCuts(noconv);
  }*/
    // }
    AliDielectronTrackCuts *trackCuts=new AliDielectronTrackCuts("ITSandgeneral_trackCuts","ITSandgeneral_trackCuts");
    //ITS related cuts
    if (cutDefinition==0)
      trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst); //does not work in MC???
      trackCuts->SetRequireTPCRefit(kTRUE);
      trackCuts->SetRequireITSRefit(kTRUE);
      diele->GetTrackFilter().AddCuts(trackCuts);
  
  //Pt cut ----------------------------------------------------------
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("trackkineandTPCQ","trackkine_and_TPCQ");
  if ((cutDefinition==0)){
    pt->AddCut(AliDielectronVarManager::kPt,1.,1e30);
  }else{ pt->AddCut(AliDielectronVarManager::kP,.8,1e30);}
  pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
  // SPDany
  if (cutDefinition==2) pt->AddCut(AliDielectronVarManager::kITSLayerFirstCls,-0.5,1.5);
  //AOD additions since there are no AliESDtrackCuts -----------------
  //
  // TPC #clusteres cut
  if (!(cutDefinition==1) ) pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);//does not work in MC???
   if (!(cutDefinition==4) )pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);//to be checked
  if (!(cutDefinition==1) ) pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);//does not work in MC???
  //TODO: DCA cuts to be investigated!!! NOTE: why?? (mwinn, 15.01.2013)
  if (!(cutDefinition==4) ) pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
  if (!(cutDefinition==4) ) pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  
  diele->GetTrackFilter().AddCuts(pt);
    
  // PID cuts --------------------------------------------------------
  if(cutDefinition ==0 || cutDefinition ==1 ||cutDefinition ==2){
    AliDielectronPID *pid = new AliDielectronPID("PID10","TPC nSigma |e|<3 + |Pi|>3.5 + P>3");
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.5,3.5,0.,0.,kTRUE);//does not work in MC???
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-20.,3.,0.,0.,kTRUE);//does not work in MC???
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
}

//______________________________________________________________________________________
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //
  // conversion rejection
  //Double_t gCut = 0.05;             // default
   Double_t gCut = 0.100;             // default
   if(!cutDefinition==4){
     AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
     gammaCut->AddCut(AliDielectronVarManager::kM,0.,gCut);
     diele->GetPairPreFilter().AddCuts(gammaCut);
  diele->SetPreFilterUnlikeOnly();
   }
 
   //Invariant mass and rapidity selection
   if(!cutDefinition==4){
     AliDielectronVarCuts *pairCut=new AliDielectronVarCuts("|Y|<.9","|Y|<.9");
     // pairCut->AddCut(AliDielectronVarManager::kM,2.,4.);
     pairCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
     diele->GetPairFilter().AddCuts(pairCut);
   }
}
//______________________________________________________________________________________

void SetupV0Cuts(AliDielectron *diele, Int_t cutDefinition)
{
  //
  // Setup the V0 cuts
  //
 if(!cutDefinition==4){  

  AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");
  gammaV0Cuts->SetPdgCodes(22,11,11);
  gammaV0Cuts->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
  gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);//to be checked, if properly filled
  gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
  //  gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle,              0.0,   0.1, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  //  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // not sure if it works as expected
  gammaV0Cuts->SetExcludeTracks(kTRUE);
  gammaV0Cuts->Print();
  
  //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
 //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; &&  const Double_t cutQTG2 < 0.04;

  diele->GetTrackFilter().AddCuts(gammaV0Cuts);
	}

}

//______________________________________________________________________________________

void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition)
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
 
  
  
  //add histograms to event class
  if (cutDefinition==0) {
    histos->AddClass("Event");
    histos->AddClass("Event_noCuts"); 
   histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    // nAcc
   histos->UserHistogram("Event","NAccRaw","Accepted raw SPD tracklets, |y|<1; nTrackl; #Entries",301,-0.5,300.5, AliDielectronVarManager::kNaccTrckltsEsd10); 
   // histos->UserHistogram("Event","NAccCorr","Accepted corr SPD tracklets, |y|<1; nTrackl; #Entries",101,-0.5,100.5, AliDielectronVarManager::kNaccTrckltsEsd10Corr); 
   // nAcc vs Zvtx
   // histos->UserHistogram("Event","NAccRaw_vs_Zvtx","Accepted raw SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,101,-0.5,100.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10);
   // histos->UserHistogram("Event","NAccCorr_vs_Zvtx","Accepted corr SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,101,-0.5,100.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10Corr);

 
   // no event cuts 
    histos->UserHistogram("Event_noCuts","VtxZ","Vertex Z;Z[cm]",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    // nAcc
    histos->UserHistogram("Event_noCuts","NAccRaw","Accepted raw SPD tracklets, |y|<1; nTrackl; #Entries",301,-0.5,300.5, AliDielectronVarManager::kNaccTrckltsEsd10);
    histos->UserHistogram("Event_noCuts","NAccRaw |y|<1 eta phi","Eta; Phi; nTrackl; #Entries",100,-1,1,144,0,6.285,301,-0.5,300.5,AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi, AliDielectronVarManager::kNaccTrckltsEsd10);
   //  histos->UserHistogram("Event_noCuts","NAccCorr","Accepted corr SPD tracklets, |y|<1; nTrackl; #Entries",101,-0.5,100.5, AliDielectronVarManager::kNaccTrckltsEsd10Corr);
   // nAcc vs Zvtx
   //  histos->UserHistogram("Event_noCuts","NAccRaw_vs_Zvtx","Accepted raw SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,101,-0.5,100.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10);
   // histos->UserHistogram("Event_noCuts","NAccCorr_vs_Zvtx","Accepted corr SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,101,-0.5,100.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10Corr);

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
  // histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;PIN [GeV];TPC number of sigmas Kaons;#tracks",
  //                        200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
   histos->UserHistogram("Track","TOFbeta_P","TOF beta;P [GeV];TOF beta;#tracks",
                      200,0.2,20.,100,0.,1.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
  // histos->UserHistogram("Track","TOFnSigmaEle_P","TOF number of sigmas Electrons;P [GeV];TOF number of sigmas Electrons;#tracks",
  //                      200,0.2,20.,200,-20.,20.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
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
  //  histos->UserHistogram("Pair","InvMass_NaccCorr_PtJpsi","Inv.Mass - NaccCorr - PtJpsi;Inv. Mass [GeV];NaccCor; pTJpsi[GeV/c]", 125,0.,125*.04,101,-0.5,100.5,100,0.,10., AliDielectronVarManager::kM,AliDielectronVarManager::kNaccTrckltsEsd10Corr, AliDielectronVarManager::kPt);
 
  diele->SetHistogramManager(histos);
}


void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition)
{//number of dimensions also not the problem
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(diele->GetName(),diele->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0,1.0,1.3,2.0, 3.0,5.0, 7.0,10.0,100.0");
  cf->AddVariable(AliDielectronVarManager::kY,"-1,-0.9,-0.8,-0.5,-0.3,0,0.3,0.5,0.8,0.9,1.0");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPseudoProperTime,150,-0.3,0.3);
  cf->AddVariable(AliDielectronVarManager::kPseudoProperTimeErr,600,0.,0.3);
  cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);
  // cf->AddVariable(AliDielectronVarManager::kArmAlpha,"-10.,-5.,-2.,-.1.,0.,1.,2.,5.,10.");
  //cf->AddVariable(AliDielectronVarManager::kArmPt,"0.0,0.4,0.6,0.8,1.0,1.5,2.0,3.0");
  

  //leg variables
  //TODO: add variables used for cuts, debug-tree??
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 1.0, 1.1, 1.2, 1.3, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 60, 65, 70, 75, 80, 85, 90, 100, 120, 160",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignalN,"0, 50, 60, 70, 80, 90, 100, 160",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCchi2Cl,"0, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.",kTRUE);  
  cf->AddVariable(AliDielectronVarManager::kEta, "-1.0, -0.9, -0.8,0.0,0.8, 0.9, 1.0", kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-2.5,-2,-1.5,-1,-0.5,4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.,3.5,4.,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.,3.5,4.,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,-0.5,5.5,kTRUE);  
  cf->AddVariable(AliDielectronVarManager::kP,"0.0, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2,1.5,2.0,3.0,4.0,5.0, 100.0",kTRUE);
  //event variables
  //cf->AddVariable(AliDielectronVarManager::kNaccTrcklts,"0.0, 9.0, 17.0, 25.0, 36.0, 55.0, 500.0");
  // cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd10,101,-0.5,100.5);
  //  cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd10Corr,101,-0.5,100.5);
  //cf->AddVariable(AliDielectronVarManager::kZvPrim,"-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.");
    Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
    if (hasMC){
      cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
      cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
      cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,10000,-5000.5,4999.5,kTRUE);

      //only steps for efficiencies//seems not to help with problem!!!!, not only MC truth step, but also MCreconstructed with Signals, no MC reconstructed in total
      //cf->SetStepsForMCtruthOnly();
    }
  
  //only in this case write MC truth info//NOTE perhaps only wirking for first one????
    // if (cutDefinition==2){
  if (hasMC) cf->SetStepForMCtruth();
  // }
  cf->SetStepsForSignal();

  diele->SetCFManagerPair(cf);
  
}

void SetSignals(AliDielectron *diele){
  AliDielectronSignalMC* promptJpsi = new AliDielectronSignalMC("promptJpsi","Prompt J/psi");   // prompt J/psi (not from beauty decays)
  promptJpsi->SetLegPDGs(11,-11);
  promptJpsi->SetMotherPDGs(443,443);//???
  promptJpsi->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);//options kUndefined, kDifferent
  promptJpsi->SetFillPureMCStep(kTRUE);
  promptJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);//kDirect produces empty Histograms...//
  promptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  diele->AddSignalMC(promptJpsi);


  AliDielectronSignalMC* beautyJpsi = new AliDielectronSignalMC("producedbeautyJpsi","produced beauty hadron -> J/psi");  // J/psi->e+e- from beauty hadron decays
  beautyJpsi->SetLegPDGs(11,-11);
  beautyJpsi->SetMotherPDGs(443,443);
  beautyJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);//test
  beautyJpsi->SetMotherSources(AliDielectronSignalMC::kPrimary,AliDielectronSignalMC::kPrimary);//added lines!!!
  beautyJpsi->SetGrandMotherSources(AliDielectronSignalMC::kPrimary,AliDielectronSignalMC::kPrimary);//added lines!!!!
  beautyJpsi->SetGrandMotherPDGs(503,503);
  beautyJpsi->SetFillPureMCStep(kTRUE); //does not help to take it out!!!
  beautyJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);//possibilities: kFinalState, kDirect, kPrimary, kDontCare, kSecondary, kNoCocktail
  // beautyJpsi->SetJpsiRadiative(AliDielectronSignalMC::kIsNotRadiative);//tested line, works, but is not solve problem!!!
  beautyJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  diele->AddSignalMC(beautyJpsi);


 AliDielectronSignalMC* background1 = new AliDielectronSignalMC("backgroundnolegjpsi","backgroundnolegjpsi");  //not J/psi->e+e-
 //background1->SetLegPDGs(11,-11);
  background1->SetMotherPDGs(443,443, kTRUE,kTRUE );
  background1->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  background1->SetGrandMotherPDGs(503,503);
  background1->SetFillPureMCStep(kFALSE);
  background1->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  background1->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(background1);

AliDielectronSignalMC* background2 = new AliDielectronSignalMC("backgroundonelegjpsi","backgroundonelegjpsi");  //not J/psi->e+e-
// background2->SetLegPDGs(11,-11);
  background2->SetMotherPDGs(443,443, kTRUE,kFALSE );
  background2->SetMothersRelation(AliDielectronSignalMC::kUndefined);
  background2->SetFillPureMCStep(kFALSE);
  background2->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  background2->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(background2);
  /*
  AliDielectronSignalMC* beautyMesonJpsi = new AliDielectronSignalMC("beautyMesonJpsi","beauty meson -> J/psi");  // J/psi->e+e- from beauty hadron decays
  beautyMesonJpsi->SetLegPDGs(11,-11);
  beautyMesonJpsi->SetMotherPDGs(443,443);
  beautyMesonJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  beautyMesonJpsi->SetGrandMotherPDGs(500,500);
  beautyMesonJpsi->SetFillPureMCStep(kTRUE);
  beautyMesonJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  beautyMesonJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  beautyMesonJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  beautyMesonJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  diele->AddSignalMC(beautyMesonJpsi);


  // physical backgrounds (electrons from other sources)
  AliDielectronSignalMC* diEleContinuum = new AliDielectronSignalMC("diEleContinuum","di-electron continuum");     // all di-electrons originating in the collision
  diEleContinuum->SetLegPDGs(11,-11);
  diEleContinuum->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleContinuum->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(diEleContinuum);

  AliDielectronSignalMC* diEleCharm = new AliDielectronSignalMC("diEleCharm","di-electrons from charm");  // dielectrons originating from charm hadrons (not neccessary from same mother)
  diEleCharm->SetLegPDGs(11,-11);
  diEleCharm->SetMotherPDGs(403,403);
  diEleCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diele->AddSignalMC(diEleCharm);

  AliDielectronSignalMC* diEleOpenCharm = new AliDielectronSignalMC("diEleOpenCharm","di-electrons from open charm");  // dielectrons originating from open charm hadrons
  diEleOpenCharm->SetLegPDGs(11,-11);
  diEleOpenCharm->SetMotherPDGs(402,402);
  diEleOpenCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diele->AddSignalMC(diEleOpenCharm);

  AliDielectronSignalMC* diEleOpenCharmJpsi = new AliDielectronSignalMC("diEleOpenCharmJpsi","1 leg from open charm + 1 jpsi leg");  // 1 leg from open charm + 1 leg from jpsi
  diEleOpenCharmJpsi->SetLegPDGs(11,-11);
  diEleOpenCharmJpsi->SetMotherPDGs(402,443);
  diEleOpenCharmJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenCharmJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharmJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diele->AddSignalMC(diEleOpenCharmJpsi);

  AliDielectronSignalMC* muonMuonPairs = new AliDielectronSignalMC("muonMuonPairs","muon+muon pairs");   // dimuon pairs, all sources (both from collision and from material interaction)
  muonMuonPairs->SetLegPDGs(13,13);
  muonMuonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(muonMuonPairs);

  // background from secondary electrons
  AliDielectronSignalMC* secondaryElectrons = new AliDielectronSignalMC("secondaryElectrons","Secondary electrons");   // all di-electrons from secondary electrons (interaction with detector)
  secondaryElectrons->SetLegPDGs(11,-11);
  secondaryElectrons->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  secondaryElectrons->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(secondaryElectrons);

  AliDielectronSignalMC* primarySecElePairs = new AliDielectronSignalMC("primarySecElePairs","Primary+Secondary electron pairs");  // primary-secondary pairs
  primarySecElePairs->SetLegPDGs(11,-11);
  primarySecElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  primarySecElePairs->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kSecondary);
  diele->AddSignalMC(primarySecElePairs);

  AliDielectronSignalMC* conversionElePairs = new AliDielectronSignalMC("conversionElePairs","conversion electron pairs");      // pairs made from conversion (may be also from 2 different conversions)
  conversionElePairs->SetLegPDGs(11,-11);
  conversionElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  conversionElePairs->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  conversionElePairs->SetMotherPDGs(22,22);
  diele->AddSignalMC(conversionElePairs);

  // misidentification
  AliDielectronSignalMC* allEleMisIdPairs = new AliDielectronSignalMC("allEleMisIdPairs","all electron+misid. pairs");  // one true electron + a mis-id electron (all sources included)
  allEleMisIdPairs->SetLegPDGs(11,11,kFALSE,kTRUE);
  allEleMisIdPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(allEleMisIdPairs);

  AliDielectronSignalMC* allMisIdMisIdPairs = new AliDielectronSignalMC("allMisIdMisIdPairs","all misid.+misid. pairs");  // mis-id + mis-id
  allMisIdMisIdPairs->SetLegPDGs(11,11,kTRUE,kTRUE);
  allMisIdMisIdPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(allMisIdMisIdPairs);

  AliDielectronSignalMC* elePionPairs = new AliDielectronSignalMC("elePionPairs","electron+pion pairs");    // true electron + mis-id pion
  elePionPairs->SetLegPDGs(11,211);
  elePionPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(elePionPairs);

  AliDielectronSignalMC* eleKaonPairs = new AliDielectronSignalMC("eleKaonPairs","electron+kaon pairs");   // true electron + mis-id kaon
  eleKaonPairs->SetLegPDGs(11,321);
  eleKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(eleKaonPairs);

  AliDielectronSignalMC* eleProtonPairs = new AliDielectronSignalMC("eleProtonPairs","Electron+proton pairs");  // true electron + mis-id proton
  eleProtonPairs->SetLegPDGs(11,2212);
  eleProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(eleProtonPairs);

  AliDielectronSignalMC* piPiPairs = new AliDielectronSignalMC("piPiPairs","pion+pion pairs");    // mis-id pion + mis-id pion
  piPiPairs->SetLegPDGs(211,211);
  piPiPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(piPiPairs);

  AliDielectronSignalMC* piKaonPairs = new AliDielectronSignalMC("piKaonPairs","pion+kaon pairs");  // mis-id pion + mis-id kaon
  piKaonPairs->SetLegPDGs(211,321);
  piKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(piKaonPairs);

  AliDielectronSignalMC* piProtonPairs = new AliDielectronSignalMC("piProtonPairs","pion+proton pairs");  // mis-id pion + mis-id proton
  piProtonPairs->SetLegPDGs(211,2212);
  piProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(piProtonPairs);

  AliDielectronSignalMC* kaonKaonPairs = new AliDielectronSignalMC("kaonKaonPairs","kaon+kaon pairs");  // mis-id kaon + mis-id kaon
  kaonKaonPairs->SetLegPDGs(321,321);
  kaonKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(kaonKaonPairs);

  AliDielectronSignalMC* kaonProtonPairs = new AliDielectronSignalMC("kaonProtonPairs","kaon+proton pairs");   // mis-id kaon + mis-id proton
  kaonProtonPairs->SetLegPDGs(321,2212);
  kaonProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(kaonProtonPairs);

  AliDielectronSignalMC* protonProtonPairs = new AliDielectronSignalMC("protonProtonPairs","proton+proton pairs");  // mis-id proton + mis-id proton
  protonProtonPairs->SetLegPDGs(2212,2212);
  protonProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(protonProtonPairs);

  AliDielectronSignalMC* muonAllPairs = new AliDielectronSignalMC("muonAllPairs","muon+everything pairs");        // mis-id muon + something else (electron, pion, kaon, proton)
  muonAllPairs->SetLegPDGs(13,13,kFALSE,kTRUE);
  muonAllPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diele->AddSignalMC(muonAllPairs);
  */
  
}
