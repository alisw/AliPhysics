
void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

AliESDtrackCuts *SetupESDtrackCuts();

TString names=("no_cuts;track_quality;tq+mcPIDele;tq+4#sigma_dEdx_E");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();

AliDielectron* ConfigJpsi2ee(Int_t cutDefinition=1)
{
  //
  // Setup the instance of AliDielectron
  //
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die =
    new AliDielectron(Form("%s",name.Data()),
                      Form("Track cuts: %s",name.Data()));
  
  // cut setup
  SetupTrackCuts(die,cutDefinition);
  SetupPairCuts(die,cutDefinition);
  
  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  //  dielectron framework histograms will be filled
  //
  InitHistograms(die,cutDefinition);
  InitCF(die,cutDefinition);
  
  return die;
}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //
  
  //ESD quality cuts
  if (cutDefinition>0){
    die->GetTrackFilter().AddCuts(SetupESDtrackCuts());
  }
  
  //MC pid
  if (cutDefinition==2){
    //MC pid
    AliDielectronVarCuts *mcpid=new AliDielectronVarCuts("legMCpid","Leg MC pid");
    mcpid->SetCutOnMCtruth();
    mcpid->SetCutType(AliDielectronVarCuts::kAny);
    mcpid->AddCut(AliDielectronVarManager::kPdgCode, 11);
    mcpid->AddCut(AliDielectronVarManager::kPdgCode, -11);
    die->GetTrackFilter().AddCuts(mcpid);
  }
  
  if (cutDefinition>=3){
    //ESD pid cuts (TPC nSigma)
    AliDielectronVarCuts *pid = new AliDielectronVarCuts("TPCnSigma","TPC nSigma cut");
    pid->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -4., 4.);
    die->GetTrackFilter().AddCuts(pid);
  }
  
}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //
  
  
  if (cutDefinition>0){
    //MC mother in y +-.9
    AliDielectronVarCuts *etaMC=new AliDielectronVarCuts("|MC Y|<.9","|MC Y|<.9");
    etaMC->AddCut(AliDielectronVarManager::kY,-.8,.8);
    etaMC->SetCutOnMCtruth();
    die->GetPairFilter().AddCuts(etaMC);
    
    //reject conversions
    AliDielectronVarCuts *openingAngleCut=new AliDielectronVarCuts("OpeningAngle","Opening angle > .35rad");
    openingAngleCut->AddCut(AliDielectronVarManager::kOpeningAngle,.35,4.);
    die->GetPairFilter().AddCuts(openingAngleCut);
    
  }

  if (cutDefinition==3){
    //
    //ESD pid cuts (TPC nSigma Pions)
    //
    AliDielectronPairLegCuts *tpcNSigmaPi = new AliDielectronPairLegCuts("|n#sigma#pi|>2","|n#sigma#pi|>2");
    
    AliDielectronVarCuts *pidPiM = new AliDielectronVarCuts("TPCnSigmaPi-","TPC nSigma- pi cut");
    pidPiM->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -10., -2.);
    
    AliDielectronVarCuts *pidPiP = new AliDielectronVarCuts("TPCnSigmaPi+","TPC nSigma+ pi cut");
    pidPiP->AddCut(AliDielectronVarManager::kTPCnSigmaPio, 2., 10.);
    
    tpcNSigmaPi->GetLeg1Filter().AddCuts(pidPiM);
    tpcNSigmaPi->GetLeg1Filter().AddCuts(pidPiP);
    
    tpcNSigmaPi->GetLeg2Filter().AddCuts(pidPiM);
    tpcNSigmaPi->GetLeg2Filter().AddCuts(pidPiP);
    
    die->GetPairFilter().AddCuts(tpcNSigmaPi);
    
    //
    //TRD pid quality cuts
    //
    AliDielectronVarCuts *trdNtrklet = new AliDielectronVarCuts("TRDpidQuality","TRD pid quality");
    trdNtrklet->AddCut(AliDielectronVarManager::kTRDpidQuality, 3.5, 8.);
    
    AliDielectronPairLegCuts *trdNtrkletAny = new AliDielectronPairLegCuts("TRDntrlt>=4 any","TRDntrlt>=4 any");
    trdNtrkletAny->GetLeg1Filter().AddCuts(trdNtrklet);
    trdNtrkletAny->GetLeg2Filter().AddCuts(trdNtrklet);
    trdNtrkletAny->SetCutType(AliDielectronPairLegCuts::kAnyLeg);
    die->GetPairFilter().AddCuts(trdNtrkletAny);
    
    AliDielectronPairLegCuts *trdNtrkletBoth = new AliDielectronPairLegCuts("TRDntrlt>=4 both","TRDntrlt>=4 both");
    trdNtrkletBoth->GetLeg1Filter().AddCuts(trdNtrklet);
    trdNtrkletBoth->GetLeg2Filter().AddCuts(trdNtrklet);
    die->GetPairFilter().AddCuts(trdNtrkletBoth);
  }
}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts()
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
  
  esdTrackCuts->SetMaxDCAToVertexZ(3.0);
  esdTrackCuts->SetMaxDCAToVertexXY(1); 
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  
  esdTrackCuts->SetMinNClustersTPC(75);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  
  return esdTrackCuts;
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram classes
  AliDielectronHistos *histos=
    new AliDielectronHistos(die->GetName(),
                            die->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  
  //Event class
  histos->AddClass("Event");
  
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
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",
                        1,0.,1.,AliDielectronVarManager::kNevents);
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        400,1e-2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  histos->UserHistogram("Track","TPCnSigma_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        400,1e-2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        500,0.,4.,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-2.,2.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMass_OpeningAngle","Opening angle:Inv.Mass;Inv. Mass [GeV];angle",
                        100,0.,4.,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  
  die->SetHistogramManager(histos);
}

void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,50,0,10);
  cf->AddVariable(AliDielectronVarManager::kEta,40,-2,2);
  cf->AddVariable(AliDielectronVarManager::kY,40,-2,2);
  cf->AddVariable(AliDielectronVarManager::kM,100,0,4);
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  cf->AddVariable(AliDielectronVarManager::kOpeningAngle,10,0,3.15);
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kEta,40,-2,2,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kImpactParXY,50,-.1,.1,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,12,-3,3,kTRUE);
  
  //no cuts
  //only in this case write MC truth info
  if (cutDefinition==0){
    cf->SetStepForMCtruth();
    cf->SetStepForNoCutsMCmotherPid();
    cf->SetStepForAfterAllCuts(kFALSE);
  }
  
  if (cutDefinition==1){
    cf->SetStepForNoCutsMCmotherPid();
  }
  
  if (cutDefinition>1){
    cf->SetStepsForEachCut();
  }
  
  if (cutDefinition>2){
    cf->SetStepsForCutsIncreasing();
  }
  
  die->SetCFManagerPair(cf);
  
}
