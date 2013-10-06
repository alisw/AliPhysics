namespace ConfDef {

void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

enum {kDefault,kQA};
}
//
TString names=("default");
TObjArray *arrNames=names.Tokenize(";");

const Int_t nDie=arrNames->GetEntries();

AliDielectron* ConfigDefault(Int_t cutDefinition)
{
  //
  // Setup the instance of AliDielectron
  //
  
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die =
    new AliDielectron(Form("%s",name.Data()),
                      Form("Track cuts: %s",name.Data()));

  if (hasMC) ConfDef::SetupMCsignals(die);
  // cut setup
  ConfDef::SetupTrackCuts(die,cutDefinition);  
  //
  ConfDef::SetupPairCuts(die,cutDefinition);
  //
  ConfDef::InitHistograms(die,cutDefinition);
  //
  ConfDef::InitCF(die,cutDefinition);
  //
//   AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
//   mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-8,-5,0,5,8,10");
//   mix->AddVariable(AliDielectronVarManager::kPhi ,"0,3.14,6.3");
//   mix->SetDepth(10);
//  die->SetMixingHandler(mix);
  //
  return die;
}

namespace ConfDef{
//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);

  //default quality cuts
  AliDielectronTrackCuts *cut1=new AliDielectronTrackCuts("cut1","cut1");
  cut1->SetRequireITSRefit(kTRUE);
  cut1->SetRequireTPCRefit(kTRUE);
  cut1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  cuts->AddCut(cut1);

  //pt and kink mother
  AliDielectronVarCuts *cut2 = new AliDielectronVarCuts("cut2Cut","cut2 cut");
  cut2->AddCut(AliDielectronVarManager::kPt,            0.8 ,     1.e30 );
  cut2->AddCut(AliDielectronVarManager::kImpactParZ,   -3.  ,     3.    );
  cut2->AddCut(AliDielectronVarManager::kImpactParXY,  -1.  ,     1.    );
  cut2->AddCut(AliDielectronVarManager::kEta,          -0.9 ,     0.9   );
  cut2->AddCut(AliDielectronVarManager::kNclsTPC,      70.  ,  160.     );
  cut2->AddCut(AliDielectronVarManager::kTPCchi2Cl,     0.  ,    4.     );
  //   cut2->AddCut(AliDielectronVarManager::kITSchi2Cl,     0.  ,   36.     ); // removes all tracks in AOD!! why???
  cut2->AddCut(AliDielectronVarManager::kKinkIndex0,    0.              );
  cut2->AddCut(AliDielectronVarManager::kTPCnSigmaEle,  -3. ,    3.     );
  cut2->AddCut(AliDielectronVarManager::kTPCnSigmaPio,   3.5, 1000.     );
  cut2->AddCut(AliDielectronVarManager::kTPCnSigmaPro,   3. , 1000.     );
  cuts->AddCut(cut2);
  
  
  //exclude conversion electrons selected by the tender
//   AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
//   noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
//   cuts->AddCut(noconv);
}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //
  
  
  // add conversion rejection
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,.05);
  die->GetPairPreFilter().AddCuts(gammaCut);
  die->SetPreFilterUnlikeOnly();
  
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
  histos->UserHistogram("Track","","",
                        1000,0.,0.,AliDielectronVarManager::kKinkIndex0);
  
  histos->UserHistogram("Track","","",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","","",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  
  //add histograms to Pair classes
  histos->UserHistogram("Pair","","",
                        301,-.01,6.01,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","","",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","","",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);

  die->SetHistogramManager(histos);
  
}

//______________________________________________________________________________________
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
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 1.0, 1.1, 1.2, 1.3, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-1.,-0.9,-0.8,0,0.8,0.9,1.0",kTRUE);
  
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

//______________________________________________________________________________________
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

}
