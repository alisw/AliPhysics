 void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition);
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition);
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition);
void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition);
void SetupV0Cuts(AliDielectron *diele, Int_t cutDefinition);

//blabla
//blabla
TString namesDieleData=("basicQ+ITS012+p>.8+nostrongexclusionPID;basicQ+ITS012+p>.8+noexclusionPID+nomassprefilter;basicQ+ITS012+p>.8+noexclusionPID+nomassprefilter_differentmixing");//;basicQ+ITS012+p>.8+noexclusionPID+nomassprefilter+activevolumecut");//;noITS");//;basicQ+ITS012+p>.8+noexclusionPID+nomassprefilterTOFifavail");//;basicQ+SPDany+pt>1+PID;basicQ+ITS012+pt>1+PID;basicQ+SPDany+pt>1+noexclPIDforcontrolpurpose;basicQ+SPDfirst+p>1.+PID;basicQ+SPDfirst+noPIDcuts+nopairing");// basicQ+SPDany+pt>1+p>1.2+PIDrequirementsHFEnoexclPIDforcontrolpurpose;basicQ+SPDany+pt>1+p>1.2+PIDrequirementsHFE+PID;basicQ+SPDany+p>1.+PID");
//TString namesDieleData=("basicQ+SPDfirst+pt>1+PID");

TObjArray *arrNamesDieleData=namesDieleData.Tokenize(";");

const Int_t nDie=arrNamesDieleData->GetEntries();

AliDielectron* ConfigJpsi2ee(Int_t cutDefinition) 
{
  //
  // Setup the instance of AliDielectron
  //
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);;
  if (cutDefinition<arrNamesDieleData->GetEntriesFast()){
    name=arrNamesDieleData->At(cutDefinition)->GetName();
  }
  AliDielectron *diele = new AliDielectron(Form("%s",name.Data()),
                             Form("Track cuts: %s",name.Data()));
  //setter to enable/disable KF-package
  diele->SetUseKF(kFALSE);
  // estimators filename
  //NOTE: what does this mean?: estimator for pp multiplicity, not needed for instance for my pA-purpose(mwinn 16.1.2012)..
  //diele->SetEstimatorFilename("$ALICE_PHYSICS/PWGDQ/dielectron/files/estimators.root");
  //diele->SetEstimatorFilename("TRAIN_ROOT/mwinn_jpsiCorr/estimators_pPb.root");
  //  diele->SetEstimatorFilename("$TRAIN_ROOT/mwinn_jpsiCorr/estimators_pPb.root");
  //cout <<"TRAIN_ROOT" << TRAIN_ROOT << endl;
  //diele->SetEstimatorFilename("estimators.root");
  // cut setup
  SetupTrackCutsDieleData(diele, cutDefinition);
  SetupPairCutsDieleData(diele, cutDefinition);
  SetupV0Cuts(diele, cutDefinition);
  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielelectron framework histograms will be filled
  //
    InitHistogramsDieleData(diele, cutDefinition);
  
  // the last definition uses no cuts and only the QA histograms should be filled!, now for all cuts

    InitCFDieleData(diele, cutDefinition);//first 2 cut sets in 3rd included
  //for last CF-Container has to set no pairing condition...


    //outcommented at the 26.11.2013
    //  AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
    // rot->SetConeAnglePhi(TMath::Pi());
    //rot->SetIterations(20);
    // diele->SetTrackRotator(rot);
  
  //outcommented for central train!! 10.09.2013
  
    AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
    mix->SetMixType(AliDielectronMixingHandler::kAll);
    if(cutDefinition==3)  mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 ,0, 2.5, 5., 7.5 , 10.");
    mix->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd10, "0.0,10.0, 20.0, 30.0, 40.0, ,50.0, 60.0, 80.0, 10000.0");//0-10,10-20,20-30,30-40,40-50,50-60, 60-80,80
    //NOTE: other mixing classes probably needed!!!
    mix->SetDepth(20);//need more jobs in one run in order to satisfy Pool-depth requirements, default settings not enough...
    diele->SetMixingHandler(mix);
    
   //no QA-class by Julian (attention introduced on he 15.06. 2013, will produce crashes with earlier Aliroot versions)
   diele->SetCutQA();  
   //   if(cutDefinition==3)diele->SetNoPairing();
  return diele;
}

//______________________________________________________________________________________
void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition /*, Bool_t isAOD*/)
{
  //
  // Setup the track cuts
  //
  
  
  
  AliDielectronTrackCuts *ITSrefit=new AliDielectronTrackCuts("ITSrefit","ITSrefit");
  ITSrefit->SetRequireITSRefit(kTRUE);
  diele->GetTrackFilter().AddCuts(ITSrefit);
  
  AliDielectronTrackCuts *TPCrefit=new AliDielectronTrackCuts("TPCrefit","TPCrefit");
  TPCrefit->SetRequireTPCRefit(kTRUE);
  diele->GetTrackFilter().AddCuts(TPCrefit);  

  //Pt cut ----------------------------------------------------------
  AliDielectronVarCuts *p = new AliDielectronVarCuts("mom","mom");
  p->AddCut(AliDielectronVarManager::kPt,0.8,1e30);
  
  diele->GetTrackFilter().AddCuts(p);
  
  
  AliDielectronVarCuts *kink = new AliDielectronVarCuts("kink","kink");
  kink->AddCut(AliDielectronVarManager::kKinkIndex0,.000001,1e30,kTRUE);
  diele->GetTrackFilter().AddCuts(kink);
  

  AliDielectronVarCuts *sPD = new AliDielectronVarCuts("SPDfirstlayercut","SPDfirstlayercut");
  sPD->AddCut(AliDielectronVarManager::kITSLayerFirstCls,-0.5,1.5);//change 30.5.2013
  if(cutDefinition==2){
     sPD->AddCut(AliDielectronVarManager::kITSLayerFirstCls,-0.5,0.5);//change 30.5.2013
  }
  diele->GetTrackFilter().AddCuts(sPD);


 // TPC #clusteres cut
 AliDielectronVarCuts *clustercut = new AliDielectronVarCuts("60cluster","60cluster");
 
 clustercut->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);//ATTENTION: wider in order to vary on CFContainer level

  diele->GetTrackFilter().AddCuts(clustercut);
  
  AliDielectronVarCuts *chisquarecut = new AliDielectronVarCuts("chisquarecut4","chisquarecut4");
  chisquarecut->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
  //cut to be added for comparison with AODs: kITSchi2Cl
  //NOTE next cut added in order to be compatible with AODs!!!!!
  // chisquarecut->AddCut(AliDielectronVarManager::kITSchi2Cl,0.,36.);
  diele->GetTrackFilter().AddCuts(chisquarecut);

  AliDielectronVarCuts *etacut = new AliDielectronVarCuts("etacut","etacut");
  
  etacut ->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);//ATTENTION: wider in order to vary on CFContainer level
  diele->GetTrackFilter().AddCuts(etacut);
  
  AliDielectronVarCuts *dcacut = new AliDielectronVarCuts("dcacut","dcacut");
  //TODO: DCA cuts to be investigated!!! NOTE: why?? (mwinn, 15.01.2013)
  dcacut->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
  dcacut->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  //for comparison with AODs cut at -2,2 for z-vertex coordinate
  //NOTE next cut added in order to be compatible with AODs!!!!!
  //  dcacut->AddCut(AliDielectronVarManager::kImpactParZ,-2.,2.);
  diele->GetTrackFilter().AddCuts(dcacut);
  
    
  // PID cuts --------------------------------------------------------
    AliDielectronPID *pid = new AliDielectronPID("PID10","TPC nSigma e>-3, e<4, pion>2.5,proton>2.5");
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);//ATTENTION looser, postprocessing by CFContainer machinerie
    //problem with variable Cuts for nsigma electron in order to express proton and pion exclusion as electron cuts, and how to integrate them into CFContainer....., if not feasible different dielectron objects...
    if(cutDefinition!=3) pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.5,3.0,0.,0.,kTRUE);//attention changed at the 19th of may to 2.0 instead of 1.5
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-20.,3.0,0.,0.,kTRUE);//attention changed at the 19th of may, instead of 1.5, also important for trending
   

  if(cutDefinition==3){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-5.0,3.0,0.,0.,kTRUE);
    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
    esdTrackCuts->SetMinLengthActiveVolumeTPC();
    diele->GetTrackFilter().AddCuts(esdTrackCuts);
  }
  
  diele->GetTrackFilter().AddCuts(pid);

}

//______________________________________________________________________________________
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition /*, Bool_t isAOD*/)
{
  //
  // Setup the pair cuts
  //
  // conversion rejection
  //Double_t gCut = 0.05;             // default
  Double_t gCut = 0.2;             // default for pPb was 0.100, reduced here for studies, in CFContainer as variable
  if(cutDefinition==0){//changed at the 11.12.2013!!!!!!!!, default without prefilter cut
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("prefiltercut 100 MeV","prefiltercut 100 MeV");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,gCut);
  //add to this cut further restrictions via 
  /*  AliDielectronPairLegCuts* daughtercut = new AliDielectronPairLegCuts("cuts on daughter particles for prefilter","cuts on daughter particles for prefilter");
  AliDielectronPID *pID_pre = new AliDielectron("PIDcut prefilter", "PIDcut prefilter");
  //has to invert logic: exclude range
  //second: has to bet connecction between the daughtercut and the invariant mass cut...
  pID_pre->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.5,3.5,0.,0.,kTRUE);//attention changed at the 19th of may to 2.0 instead of 1.5
  pID_pre->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-20.,3.0,0.,0.,kTRUE);
  pID_pre->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.0,3.0);
  AliDielectronVarCuts *otherthanPID = new AliDielectronVarCuts("everything except of PID","everything except of PID");
  otherthanPID->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
  otherthanPID->AddCut(AliDielectronVarManager::kITSLayerFirstCls,-0.5,0.5);
  otherthanPID->AddCut(AliDielectronVarManager::kPt,1.0,1e30);
  daughtercut->GetLeg1Filter().AddCuts(pID_pre);
  daughtercut->GetLeg2Filter().AddCuts(pID_pre);
  daughtercut->GetLeg1Filter().AddCuts( otherthanPID);
  daughtercut->GetLeg2Filter().AddCuts( otherthanPID);*/
  //....some lines missing...
  //  AliDielectronVarCut
  diele->GetPairPreFilter().AddCuts(gammaCut);
  //  diele->GetPairPrefilter().AddCuts(daughtercut);
  //  diele->SetPreFilterUnlikeOnly();//outcommented 11.06.2013
  }
  
  //Invariant mass and rapidity selection
  AliDielectronVarCuts *pairCut=new AliDielectronVarCuts("|Y|-cut","|Y|-cut");
  
  pairCut->AddCut(AliDielectronVarManager::kY,-.9,.9);
  
  diele->GetPairFilter().AddCuts(pairCut);
}
//______________________________________________________________________________________
void SetupV0Cuts(AliDielectron *diele, Int_t cutDefinition)
{
  //
  // Setup the V0 cuts
  //

  
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
  // gammaV0Cuts->Print();
  
  //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
 //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; &&  const Double_t cutQTG2 < 0.04;

  diele->GetTrackFilter().AddCuts(gammaV0Cuts);
}


//______________________________________________________________________________________


void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition /*, Bool_t isAOD*/)
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
  for (Int_t i=0; i<10 /*for mixing until 10*/; ++i){ //
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
   histos->UserHistogram("Event","CentralityV0M","CentralityV0M;centrality_{V0M}(%)",100,0.0,100.0,AliDielectronVarManager::kCentrality);
   histos->UserHistogram("Event","CentralityCL1","CentralityCL1;centrality_{CL1}(%)",100,0.0,100.0,AliDielectronVarManager::kCentralitySPD);
   histos->UserHistogram("Event","CentralityV0A","CentralityV0A;centrality_{V0A}(%)",100,0.0,100.0,AliDielectronVarManager::kCentralityV0A);
   histos->UserHistogram("Event","CentralityV0C","CentralityV0C;centrality_{V0C}(%)",100,0.0,100.0,AliDielectronVarManager::kCentralityV0C);
   histos->UserHistogram("Event","CentralityZNA","CentralityZNA;centrality_{ZNA}(%)",100,0.0,100.0,AliDielectronVarManager::kCentralityZNA);
   
   histos->UserHistogram("Event","CentralityCL1_vs_CentralityZNA","CentralityCL1_vs_CentralityZNA;centrality_{CL1}(%);centrality_{ZNA}(%)",100,0.0,100.0,100.0,0.0,100.0,AliDielectronVarManager::kCentralitySPD,AliDielectronVarManager::kCentralityZNA);
   histos->UserHistogram("Event","CentralityCL1_vs_CentralityV0A","CentralityCL1_vs_CentralityV0A;centrality_{CL1}(%);centrality_{V0A}(%)",100,0.0,100.0,100.0,0.0,100.0,AliDielectronVarManager::kCentralitySPD,AliDielectronVarManager::kCentralityV0A);
   histos->UserHistogram("Event","CentralityCL1_vs_CentralityV0A","CentralityCL1_vs_CentralityV0M;centrality_{CL1}(%);centrality_{V0M}(%)",100,0.0,100.0,100.0,0.0,100.0,AliDielectronVarManager::kCentralitySPD,AliDielectronVarManager::kCentrality);
   histos->UserHistogram("Event","CentralityV0A_vs_AmplitudeV0A","CentralityV0A_vs_AmplitudeV0A",100,0.0,100.0,401,-0.5,400.5,AliDielectronVarManager::kCentralityV0A,AliDielectronVarManager::kMultV0A);
   histos->UserHistogram("Event","CentralityV0A_vs_AmplitudeV0Aeq","CentralityV0A_vs_AmplitudeV0Aeq",100,0.0,100.0,401,-0.5,400.5,AliDielectronVarManager::kCentralityV0A,AliDielectronVarManager::kEqMultV0A);
  histos->UserHistogram("Event","AmplitudeV0A_vs_AmplitudeV0Aeq","AmplitudeV0A_vs_AmplitudeV0Aeq",401,-0.5,400.5,401,-0.5,400.5,AliDielectronVarManager::kMultV0A,AliDielectronVarManager::kEqMultV0A);
   histos->UserHistogram("Event","CentralityCL1_vs_SPDTrckltsCorr10","CentralityCL1_vs_SPDTrckltsCorr10",100,0.0,100.0,200.0,-0.5,200.5,AliDielectronVarManager::kCentralityV0A,AliDielectronVarManager::kNaccTrckltsEsd10Corr);


    // nAcc
   //    histos->UserHistogram("Event","NAccRaw","Accepted raw SPD tracklets, |y|<1; nTrackl; #Entries",101,-0.5,100.5, AliDielectronVarManager::kNaccTrckltsEsd10); 
    //   histos->UserHistogram("Event","NAccCorr","Accepted corr SPD tracklets, |y|<1; nTrackl; #Entries",101,-0.5,100.5, AliDielectronVarManager::kNaccTrckltsEsd10Corr); 
   // nAcc vs Zvtx
     histos->UserHistogram("Event","NAccRaw_vs_Zvtx","Accepted raw SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,301,-0.5,300.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrcklts10);
 histos->UserHistogram("Event","NAccRaw_vs_Zvtx_16","Accepted raw SPD tracklets vs Z vtx, |y|<1.6; Zvtx[cm]; nTrackl ",300,-15.,15.,301,-0.5,300.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrcklts);
 histos->UserHistogram("Event","NAccRaw_vs_Zvtx_05","Accepted raw SPD tracklets vs Z vtx, |y|<0.5; Zvtx[cm]; nTrackl ",300,-15.,15.,301,-0.5,300.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd05);

  histos->UserHistogram("Event","NAccCorr_vs_Zvtx","Accepted corr. SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10Corr);
 histos->UserHistogram("Event","NAccCorr_vs_Zvtx_16","Accepted corr raw SPD tracklets vs Z vtx, |y|<1.6; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd16Corr);
 histos->UserHistogram("Event","NAccCorr_vs_Zvtx_05","Accepted corr raw SPD tracklets vs Z vtx, |y|<0.5; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd05Corr);
 histos->UserHistogram("Event","NAccV0A_vs_Zvtx","V0A amplitude vs. Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,401,-0.5,400.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kMultV0A); 

 
 histos->UserHistogram("Event","NAccTPCITS_vs_Zvtx_05","Accepted raw combined tracks vs Z vtx, |y|<0.5; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsTpcEsd05);
 histos->UserHistogram("Event","NAccTPCITS_vs_Zvtx","Accepted raw combined tracks vs Z vtx, |y|<1.0; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsTpcEsd10);
 histos->UserHistogram("Event","NAccTPCITS_vs_Zvtx_16","Accepted raw combined tracks vs Z vtx, |y|<1.6; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsTpcEsd16);

histos->UserHistogram("Event","NAccITSsa_vs_Zvtx_05","Accepted raw ITS tracks +tracklets vs Z vtx, |y|<0.5; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsPureEsd05);
 histos->UserHistogram("Event","NAccITSsa_vs_Zvtx","Accepted raw ITS tracks +tracklets vs Z vtx, |y|<1.0; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsPureEsd10);
 histos->UserHistogram("Event","NAccTPCITSsa_vs_Zvtx_16","Accepted raw ITS tracks +tracklets vs Z vtx, |y|<1.6; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsPureEsd16);




 histos->UserHistogram("Event","NAccTPCITSCorr_vs_Zvtx_05","Accepted corr combined tracks vs Z vtx, |y|<0.5; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsTpcEsd05Corr);
 histos->UserHistogram("Event","NAccTPCITSCorr_vs_Zvtx","Accepted corr combined tracks vs Z vtx, |y|<1.0; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsTpcEsd10Corr);
 histos->UserHistogram("Event","NAccTPCITSCorr_vs_Zvtx_16","Accepted corr combined tracks vs Z vtx, |y|<1.6; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsTpcEsd16Corr);

 histos->UserHistogram("Event","NAccCorr_vs_Zvtx","Accepted corr SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,101,-0.5,100.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10Corr);

 histos->UserHistogram("Event","V0ACentrality_vs_NaccTrckltsEsd10Corr","Accepted corr SPD tracklets vs V0A centrality, |y|<1; Zvtx[cm]; nTrackl ",100,0.0,100.,201,-0.5,200.5, AliDielectronVarManager::kCentralityV0A,AliDielectronVarManager::kNaccTrckltsEsd10Corr);

   
 histos->UserHistogram("Event","NaccRaw_vs_NaccCorr","Accepted raw SPD tracklets vs acc. corrected SPD tracklets, |y|<1; nTrackl(raw); nTrackl(corr)",201,-0.5,200.5,201,-0.5,200.5, AliDielectronVarManager::kNaccTrcklts10,AliDielectronVarManager::kNaccTrckltsEsd10Corr);
   histos->UserHistogram("Event","NaccRaw_vs_NaccCorr_05","Accepted raw SPD tracklets vs acc. corrected SPD tracklets, |y|<0.5;  nTrackl(raw); nTrackl(corr)",201,-0.5,200.5,201,-0.5,200.5, AliDielectronVarManager::kNaccTrckltsEsd05,AliDielectronVarManager::kNaccTrckltsEsd05Corr);
   histos->UserHistogram("Event","NaccRaw_vs_NaccCorr_16","Accepted raw SPD tracklets vs acc. corrected SPD tracklets, |y|<1.6;  nTrackl(raw); nTrackl(corr)",201,-0.5,200.5,201,-0.5,200.5, AliDielectronVarManager::kNaccTrcklts,AliDielectronVarManager::kNaccTrckltsEsd16Corr);
    

 
   // no event cuts 
   histos->UserHistogram("Event_noCuts","VtxZ","Vertex Z;Z[cm]",300,-15.,15.,AliDielectronVarManager::kZvPrim);
   histos->UserHistogram("Event_noCuts","CentralityV0M","CentralityV0M;centrality_{V0M}(%)",100,0.0,100.0,AliDielectronVarManager::kCentrality);
   histos->UserHistogram("Event_noCuts","CentralityCL1","CentralityCL1;centrality_{CL1}(%)",100,0.0,100.0,AliDielectronVarManager::kCentralitySPD);
   histos->UserHistogram("Event_noCuts","CentralityV0A","CentralityV0A;centrality_{V0A}(%)",100,0.0,100.0,AliDielectronVarManager::kCentralityV0A);
   histos->UserHistogram("Event_noCuts","CentralityV0C","CentralityV0C;centrality_{V0C}(%)",100,0.0,100.0,AliDielectronVarManager::kCentralityV0C);
   histos->UserHistogram("Event_noCuts","CentralityZNA","CentralityZNA;centrality_{ZNA}(%)",100,0.0,100.0,AliDielectronVarManager::kCentralityZNA);
    // nAcc
   //  histos->UserHistogram("Event_noCuts","NAccRaw","Accepted raw SPD tracklets, |y|<1; nTrackl; #Entries",101,-0.5,100.5, AliDielectronVarManager::kNaccTrckltsEsd10);
    histos->UserHistogram("Event_noCuts","NAccCorr","Accepted corr SPD tracklets, |y|<1; nTrackl; #Entries",201,-0.5,200.5, AliDielectronVarManager::kNaccTrckltsEsd10Corr);
   // nAcc vs Zvtx
   histos->UserHistogram("Event_noCuts","NAccRaw_vs_Zvtx","Accepted raw SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,201,-0.5,200.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10);
   histos->UserHistogram("Event_noCuts","CentralityCL1_vs_CentralityZNA","CentralityCL1_vs_CentralityZNA;centrality_{CL1}(%);centrality_{ZNA}(%)",200,0.0,200.0,100.0,0.0,100.0,AliDielectronVarManager::kCentralitySPD,AliDielectronVarManager::kCentralityZNA);
   histos->UserHistogram("Event_noCuts","CentralityCL1_vs_CentralityV0A","CentralityCL1_vs_CentralityV0A;centrality_{CL1}(%);centrality_{V0A}(%)",100,0.0,100.0,100.0,0.0,100.0,AliDielectronVarManager::kCentralitySPD,AliDielectronVarManager::kCentralityV0A);
  histos->UserHistogram("Event_noCuts","CentralityCL1_vs_CentralityV0A","CentralityCL1_vs_CentralityV0M;centrality_{CL1}(%);centrality_{V0M}(%)",100,0.0,100.0,100.0,0.0,100.0,AliDielectronVarManager::kCentralitySPD,AliDielectronVarManager::kCentrality);

   // }

   //    histos->UserHistogram("Event_noCuts","NAccCorr_vs_Zvtx","Accepted corr SPD tracklets vs Z vtx, |y|<1; Zvtx[cm]; nTrackl ",300,-15.,15.,101,-0.5,100.5, AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccTrckltsEsd10Corr);
   
  }
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",400,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","TOFPIDBit","TOFPIDBit;TOFPIDBit;#tracks",2,-0.5,1.5,AliDielectronVarManager::kTOFPIDBit);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","TPCnClsFracS","fraction of shared Clusters TPC; fraction of shared TPC clusters;#tracks",110,-0.1,1.1,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCsignalN","Number of PID Clusters TPC;TPC PID number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  histos->UserHistogram("Track","TPCchi2Cl","TPCchi2Cl",100,0,8,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","nClsoverfindablecluster","Number of found Clusters TPC over findably ;TPC number cluster over findable;#tracks",160,0.0,1.1,AliDielectronVarManager::kNFclsTPCrFrac);
  //  histos->UserHistogram("Track","TPCnCls vs PT","Number of Clusters TPC vs PT;TPC number clusteres;pT [GeV/c];#tracks",159,0.,159.,200,0.2,20,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kPt,kTRUE);
  //need more histograms: TPC-cluster vs. mult., TPC-cluster vs. z-Vertex-position, TPC nclus found/crossed, average nsigma electron vs mult. , vs. z-vertex-position, number of selected tracks per event pos, neg. number of selected tracks per event (distribution and vs. run number -> could be already done ), number of dielectron pairs vs. phi and vs rapidity, pT vs mass, average pT of pairs vs run, average invariant mass vs run (the later two can be easily done via report)
  //  histos->UserHistogram("Track","TPCnCls vs PT","Number of Clusters TPC vs PT;TPC number clusteres;pT [GeV/c];#tracks",159,0.,159.,200,0.2,20,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kPt,kTRUE);
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
  //  histos->UserHistogram("Track","TOFbeta_P","TOF beta;P [GeV];TOF beta;#tracks",
  //                    200,0.2,20.,100,0.,1.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaEle_P","TOF number of sigmas Electrons;P [GeV];TOF number of sigmas Electrons;#tracks",
                        200,0.2,20.,200,-20.,20.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);
  // histos->UserHistogram("Track","TRDnCls","Number of Clusters TRD;TRD number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTRD);
  //  histos->UserHistogram("Track","TRDntracklets","Number of tracklets TRD;TRD number tracklets;#tracks",7,-0.5,6.5,AliDielectronVarManager::kTRDntracklets);
  // histos->UserHistogram("Track","TPCnSigmaEle_P_eta","total Momentum[GeV/c] - eta - TPC number of sigmas Electrons;PIN [GeV/c];Eta; TPC number of sigmas Electrons;#tracks",
  //			200,0.2,20.,20,-1.0,1,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  //  histos->UserHistogram("Track","TRDprobEle_P","TRD electron prob.;P [GeV];TRD electron prob.;#tracks",
  //                    200,0.2,20.,100,.0,1.,AliDielectronVarManager::kP,AliDielectronVarManager::kTRDprobEle,kTRUE);
  // histos->UserHistogram("Track","TRDprobEle2D_P","TRD electron prob. 2D;P [GeV];TRD electron prob. 2D;#tracks",
  //                    200,0.2,20.,100,.0,1.,AliDielectronVarManager::kP,AliDielectronVarManager::kTRDprob2DEle,kTRUE);

  


      
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        125,0.,125*.04,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        110,-1.1,1.1,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","Pt","Pt;Pt;#pairs",
                        200,0.,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","OpeningAngletransverse","Opening angle transverse;angle",
                        100,0.,3.15,AliDielectronVarManager::kDeltaPhi);
  histos->UserHistogram("Pair","Chi2NDF","chisquareNDF;chisquare/ndf;#pairs",
                        100,0.,30.,AliDielectronVarManager::kChi2NDF);
  // histos->UserHistogram("Pair","distanceXY","distancelegsXY;distanceXY[cm];#pairs",
  //			100,0.,.0001,AliDielectronVarManager::kLegDistXY);
  // histos->UserHistogram("Pair","distance","distancelegs;distance[cm];#pairs",
  //                    100,0.,.0001,AliDielectronVarManager::kLegDist);
  //  histos->UserHistogram("Pair","pseudoproperdecaylength","pseudoproperdecaylength;pseudoproperdecaylength[cm];#pairs",
  //			100,0.,.5,AliDielectronVarManager::kPseudoProperTime);
  //  histos->UserHistogram("Pair","Armenteros-Podolanski","Armenteros-Podolanski;ArmAlpha;ArmPt[GeV];#tracks",
  //			100,-10.0,10.,100,0.,2.,AliDielectronVarManager::kArmAlpha,AliDielectronVarManager::kArmPt,kTRUE);
  //   histos->UserHistogram("Pair","InvMass_OpeningAngle_Pt","InvMass - OpeningAngle - Pt;Inv. Mass [GeV/c^2]; Opening angle;pTdielectron[GeV/c]",
  //			125,0.,125*.04, 100,0.,3.15, 100, 0.,10.,AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle,AliDielectronVarManager::kPt);
   //  histos->UserHistogram("Pair","InvMass_OpeningAngletransverse_Pt","InvMass - OpeningAngletransverse - Pt;Inv. Mass [GeV/c^2]; Opening angle transverse;pTdielectron[GeV/c]",
   //			125,0.,125*.04, 100,0.,3.15, 100, 0.,10.,AliDielectronVarManager::kM, AliDielectronVarManager::kDeltaPhi,AliDielectronVarManager::kPt,kTRUE );

  diele->SetHistogramManager(histos);
}


void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition /*, Bool_t isAOD*/)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(diele->GetName(),diele->GetTitle());
  
  //pair variables
  //  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 1.0, 1.3, 1.5, 2.0, 3.0, 5.0, 6.0, 7.0, 8.0, 10.0");
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 1.3, 1.5, 4.5, 10.0");
  //  cf->AddVariable(AliDielectronVarManager::kY,"-0.9,-0.8,-0.5,0.5,0.8,0.9");
  // cf->AddVariable(AliDielectronVarManager::kDeltaEta,20,0,TMath::Pi());//added for private train to be checked what about memory consumption
  cf->AddVariable(AliDielectronVarManager::kM,500,0.,500*.01); //20Mev Steps
  //  cf->AddVariable(AliDielectronVarManager::kPseudoProperTime,150,-0.3,0.3);
  // cf->AddVariable(AliDielectronVarManager::kPseudoProperTimeErr,600,0.,0.3);
  cf->AddVariable(AliDielectronVarManager::kPairType,"-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5");//ATTENTION: to be changed for ME
  // cf->AddVariable(AliDielectronVarManager::kArmAlpha,"-10.,-5.,-2.,-.1.,0.,1.,2.,5.,10.");
  //cf->AddVariable(AliDielectronVarManager::kArmPt,"0.0,0.4,0.6,0.8,1.0,1.5,2.0,3.0");
 

  //leg variables
  //first three variables new in order to emulate AODs
  // cf->AddVariable(AliDielectronVarManager::kImpactParZ,"-2.0,2.0", kTRUE);
  //  cf->AddVariable(AliDielectronVarManager::kITSchi2Cl,"0.,36.,100.0", kTRUE);
  //  cf->AddVariable(AliDielectronVarManager::kNFclsTPCfCross,"0.0, 0.8, 1.01", kTRUE);  
  //  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-5.0,-3.0,3.0",kTRUE);
  //
  cf->AddVariable(AliDielectronVarManager::kPt,"0.8, 1.0, 1.1",kTRUE);//ATTENTION OMITTED
    //cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 60, 65, 70, 75, 80, 90, 100, 120, 160",kTRUE);
  // cf->AddVariable(AliDielectronVarManager::kTPCsignalN,"0, 50, 60, 70, 80, 90, 100,120, 160",kTRUE);
  //  cf->AddVariable(AliDielectronVarManager::kTPCchi2Cl,"0, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.",kTRUE);  
  // cf->AddVariable(AliDielectronVarManager::kEta, "-0.9, -0.8, -0.5,  0.0, 0.5, 0.8, 0.9", kTRUE);
  //  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-3.0, -2,  3.",kTRUE);
  //  cf->AddVariable(AliDielectronVarManager::kTOFnSigmaEle,"-3.0, -2.5,-2,-1.5,-1, -0.5, 3., 4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.,3.5, 4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.,3.5, 4.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,2,-0.5,1.5,kTRUE);
  //  cf->AddVariable(AliDielectronVarManager::kPIn,"0.0,0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2,1.5,2.0,3.0,4.0,5.0, 100.0",kTRUE);
  //event variables
  //cf->AddVariable(AliDielectronVarManager::kNaccTrcklts,"0.0, 9.0, 17.0, 25.0, 36.0, 55.0, 500.0");
  cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd10Corr, "0.0,0.5,13.5, 18.5, 24.5, 29.5,35.5, 44.5, 59.5,74.5, 99.5, 129.5, 199.50");
  //  cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd16Corr, "0.0,20.8, 28.8, 38.4, 48.0,57.6, 72.0, 112.0, 160.0, 208.0, 320.00"); //to be done in a separate train run...
  //  cf->AddVariable(AliDielectronVarManager::kNaccTrckltsEsd05Corr, "0.0,6.5, 9.0, 12.0, 15.0,18.0, 22.5, 35.0, 50.0, 65.0, 100.00");//to be done in a separate train run
  cf->AddVariable(AliDielectronVarManager::kCentrality, "0.0, 20.0, 40.0, 60.0,  100.0");
  cf->AddVariable(AliDielectronVarManager::kCentralitySPD, "0.0, 20.0, 40.0, 60.0, 80.0, 100.0");
  cf->AddVariable(AliDielectronVarManager::kCentralityV0A, "0.0, 20.0, 40.0, 60.0, 80.0, 100.0");
  //cf->AddVariable(AliDielectronVarManager::kMultV0A, "0.0, 80.0, 160.0, 240.0, 320.0, 400.0");//not needed: take correspondence from centrality...
  cf->AddVariable(AliDielectronVarManager::kCentralityZNA, "0.0, 20.0, 40.0, 60.0, 80.0, 100.0");
  cf->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5.,2.5,  0., 2.5, 5., 7.5, 10.");
  //  cf->AddVariable(AliDielectronVarManager::kImpactParXY,"-1.0,-0.6,-0.4,-0.2,0.2,0.4,0.6,1.0");


  /*if (!isAOD){
    Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
    if (hasMC){
      cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
      cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
      cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,10000,-5000.5,4999.5,kTRUE);
      }
  }*/

  
  //only in this case write MC truth info
  //if (cutDefinition==0){
  //  cf->SetStepForMCtruth();
  //}

  diele->SetCFManagerPair(cf);
  
}

