void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void AddMCSignals(AliDielectron* die);

void SetupCuts(AliDielectron *die, Int_t cutDefinition);

AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);
AliDielectronPID *SetPIDcuts(Int_t cutDefinition);
const AliDielectronEventCuts *GetEventCuts();

Bool_t isRandomRejTask=kFALSE;//needed for InitHistograms() //dont change!!!

const TString  RndmPtExpr="exp(-x/3.)";
const Double_t RndmPtMin=0.2;
const Double_t RndmPtMax=10.0;
const Double_t RndmEtaMax=0.8;
const Int_t nTestpartPerEle=500; // number of testparticles used per final analysis electron in an event.


TString names("pt200_cut1;pt200_cut2;pt200_cut3;pt200_cut4;pt200_cut5;pt200_cut6;pt200_cut7;pt200_cut8;pt200_cut9;pt200_cut10;pt200_cut11;pt200_cut12;pt200_cut13;pt200_cut14;pt200_cut15;pt200_cut16;pt200_cut17;pt200_cut18;pt200_cut19;pt200_cut20");

     Bool_t kRot = 0;
     Bool_t kMix = 1;

     Bool_t randomizeDau = kFALSE;
     Bool_t hasMC;
     
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();
Bool_t MCenabled=kFALSE;
const Int_t nPF = 50;

AliDielectron* Config_tbroeker_lowmass(Int_t cutDefinition=1, Bool_t isRandomRej=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  isRandomRejTask=isRandomRej;
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if(cutDefinition<arrNames->GetEntriesFast())
    name=arrNames->At(cutDefinition)->GetName();

  AliDielectron *die = new AliDielectron(Form("%s",name.Data()),Form("Track cuts: %s",name.Data()));
  if(hasMC) AddMCSignals(die);
  
	if(kRot){
	  AliDielectronTrackRotator *rot = new AliDielectronTrackRotator;
	  rot->SetConeAnglePhi(TMath::Pi());
	  rot->SetIterations(10);
	  die->SetTrackRotator(rot);
	}//kRot
 
	if(kMix){
  	AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
	  mix->SetMixType(AliDielectronMixingHandler::kAll);
  	mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
    mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
  	mix->SetDepth(20);
	  die->SetMixingHandler(mix);
	}//kMix


	// set track cuts
  SetupCuts(die,cutDefinition);

  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //

  InitHistograms(die,cutDefinition);
//  InitCF(die,cutDefinition);

  die->SetNoPairing(kFALSE);
  
  return die;
}

//______________________________________________________________________________________
void SetupCuts(AliDielectron *die, Int_t cutDefinition)
{
  // Setup the track cuts

  AliDielectronV0Cuts *noconv = new AliDielectronV0Cuts("IsGamma","IsGamma");
  // which V0 finder you want to use
  noconv->SetV0finder(AliDielectronV0Cuts::kAll);  // kAll(default), kOffline or kOnTheFly
  // add some pdg codes (they are used then by the KF package and important for gamma conversions)
  noconv->SetPdgCodes(22,11,11); // mother, daughter1 and 2
  // add default PID cuts (defined in AliDielectronPID)
  // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
  //noconv->SetDefaultPID(16, AliDielectronV0Cuts::kAny);
  // add the pair cuts for V0 candidates
  noconv->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.00, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.00, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kR,                             3.0,  90.00, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kM,                             0.0,   0.10, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);

  // selection or rejection of V0 tracks
  noconv->SetExcludeTracks(kTRUE);

  //pairing with TLorentzVector
  die->SetUseKF(kFALSE);
           
  AliDielectronVarCuts* pairCutsInvM  = new AliDielectronVarCuts("pairCutsInvM","pairCutsInvM");
  AliDielectronVarCuts* pairCutsOpAng = new AliDielectronVarCuts("pairCutsOpAng","pairCutsOpAng");  
        
  if(cutDefinition < nPF){

    die->GetTrackFilter().AddCuts(SetupPreFilterESDtrackCuts(cutDefinition));
    die->GetPairPreFilterLegs().AddCuts(noconv);
	//pairPrefilter
//  	AliAnalysisCuts* pairPreCuts=0x0;

//    pairCutsInvM ->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
    pairCutsOpAng->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.050); 
    
//    AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
//    pairCutsCG->AddCut(pairCutsInvM);
//    pairCutsCG->AddCut(pairCutsOpAng);
//    pairPreCuts = pairCutsCG;
    
    die->GetPairPreFilter().AddCuts(pairCutsOpAng);
    
	  //FinalTrackCuts after prefiltering
  	die->GetPairPreFilterLegs().AddCuts(SetPIDcuts(cutDefinition));
	  die->GetPairPreFilterLegs().AddCuts(SetupESDtrackCuts(cutDefinition));
    die->GetPairPreFilterLegs().AddCuts(noconv);

    die->SetPreFilterUnlikeOnly(kTRUE);
  }
  else{
    die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));
    die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));
    die->GetTrackFilter().AddCuts(noconv);
  }
//  AliDielectronVarCuts* finalPairCutsOpAng = new AliDielectronVarCuts("finalPairCutsOpAng","finalPairCutsOpAng");
//  finalPairCutsOpAng->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.050,kTRUE); 
//  die->GetPairFilter().AddCuts(finalPairCutsOpAng);        
}
//______________________________________________________________________________________
//-----------------------------------pid------------------------------------------------

AliDielectronPID *SetPIDcuts(Int_t cutDefinition){
  
  AliDielectronPID *pid = new AliDielectronPID();

  if(cutDefinition == 0){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }	
  if(cutDefinition == 1){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 2){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 3){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -4. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3.5,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 4){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3.5,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 5){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }	
  if(cutDefinition == 6){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 7){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 8){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3.5,0. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 9){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 10){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 11){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 12){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 13){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 14){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 15){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 16){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 17){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 18){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 19){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition >= 20){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    if(cutDefinition != 21 && cutDefinition != 23 && cutDefinition != 24 && cutDefinition != 25) 
      pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    if(cutDefinition == 20 || cutDefinition == 21){
      pid->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -1.5,1.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
      pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,  -1.5,1.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    }
    if(cutDefinition == 22 || cutDefinition == 23){
      pid->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -1.8,1.8,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
      pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,  -1.8,1.8,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    }
  }
 return pid;
}

//______________________________________________________________________________________
//-----------------------------------track cuts-----------------------------------------
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition){

  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts();

  //global
  fesdTrackCuts->SetPtRange( 0.2 , 100. );
  fesdTrackCuts->SetEtaRange( -0.8 , 0.8 );
  fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fesdTrackCuts->SetDCAToVertex2D(kFALSE);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.);
  fesdTrackCuts->SetMaxDCAToVertexXY(1.);
 
  fesdTrackCuts->SetRequireTPCRefit(kTRUE);
  fesdTrackCuts->SetRequireITSRefit(kTRUE);

  if(cutDefinition == 0){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 1){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(130);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 2){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(80);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 3){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(120);
    fesdTrackCuts->SetMinNCrossedRowsTPC(130);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 4){
    fesdTrackCuts->SetMinNClustersITS(6);
    fesdTrackCuts->SetMaxChi2PerClusterITS(2.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(80);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 5){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 6){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(80);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 7){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(120);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 8){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(120);
    fesdTrackCuts->SetMinNCrossedRowsTPC(80);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 9){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 10){
    fesdTrackCuts->SetMinNClustersITS(5);
//    fesdTrackCuts->SetMaxChi2PerClusterITS(6);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(120);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 11){
    fesdTrackCuts->SetMinNClustersITS(5);
//    fesdTrackCuts->SetMaxChi2PerClusterITS(100);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(120);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 12){
    fesdTrackCuts->SetMinNClustersITS(6);
//    fesdTrackCuts->SetMaxChi2PerClusterITS(100);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 13){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 14){
    fesdTrackCuts->SetMinNClustersITS(6);
    fesdTrackCuts->SetMaxChi2PerClusterITS(2.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 15){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 16){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
//    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 17){
    fesdTrackCuts->SetMinNClustersITS(3);
//    fesdTrackCuts->SetMaxChi2PerClusterITS(4);
//    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 18){
    fesdTrackCuts->SetMinNClustersITS(4);
//    fesdTrackCuts->SetMaxChi2PerClusterITS(6);    
//    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 19){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(120);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition >= 20){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }

  return fesdTrackCuts;
}
//______________________________________________________________________________________
//-------------------------------prefilter pid------------------------------------------
//AliDielectronPID *SetPreFilterPIDcuts(Int_t cutDefinition){
//  AliDielectronPID *pid = new AliDielectronPID();
//    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
//    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
//    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
//    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);   
// return pid;
//}

//______________________________________________________________________________________
//-------------------------------prefilter track cuts-----------------------------------
AliESDtrackCuts *SetupPreFilterESDtrackCuts(Int_t cutDefinition){

  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts;

  //global
  fesdTrackCuts->SetPtRange( 0.08 , 100. );
  fesdTrackCuts->SetEtaRange( -1.1 , 1.1 );
  fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fesdTrackCuts->SetDCAToVertex2D(kFALSE);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.);
  fesdTrackCuts->SetMaxDCAToVertexXY(1.);
 
  //ITS cuts
  fesdTrackCuts->SetRequireITSRefit(kTRUE);
  fesdTrackCuts->SetMinNClustersITS(3);
  fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);

  return fesdTrackCuts;
  
}
//______________________________________________________________________________________
//-------------------------------adding MC signals-----------------------------------
void AddMCSignals(AliDielectron* die){
  
  AliDielectronSignalMC* lmee_same = new AliDielectronSignalMC("lmee_same","dielectrons same mother");
  lmee_same->SetLegPDGs(11,-11);
  lmee_same->SetCheckBothChargesLegs(kTRUE,kTRUE);
  lmee_same->SetLegSources(AliDielectronSignalMC::kFinalState,AliDielectronSignalMC::kFinalState);
  lmee_same->SetMothersRelation(AliDielectronSignalMC::kSame);
  lmee_same->SetMotherPDGs(22,22,kTRUE,kTRUE);
  
  AliDielectronSignalMC* lmee_diff = new AliDielectronSignalMC("lmee_diff","dielectrons different mother");
  lmee_diff->SetLegPDGs(11,-11);
  lmee_diff->SetCheckBothChargesLegs(kTRUE,kTRUE);
  lmee_diff->SetLegSources(AliDielectronSignalMC::kFinalState,AliDielectronSignalMC::kFinalState);
  lmee_diff->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  lmee_diff->SetMotherPDGs(22,22,kTRUE,kTRUE);
                    
  AliDielectronSignalMC* conv = new AliDielectronSignalMC("conversions","conversions same mother");
  conv->SetLegPDGs(11,-11);
  conv->SetCheckBothChargesLegs(kTRUE,kTRUE);
  conv->SetLegSources(AliDielectronSignalMC::kSecondary,AliDielectronSignalMC::kSecondary);
  conv->SetMothersRelation(AliDielectronSignalMC::kSame);
  conv->SetMotherPDGs(22,22,kFALSE,kFALSE);
    
  die->AddSignalMC(lmee_same);
  die->AddSignalMC(lmee_diff);
  die->AddSignalMC(conv);
}
//______________________________________________________________________________________
//-------------------------------Histogram definition-----------------------------------
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
  if(die->GetMCSignals()) for(Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
      
  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0,1,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","ZVertex","ZVertex;ZVertex/cm",120,-12.,12.,AliDielectronVarManager::kZvPrim);

  //add histograms to track class
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","P","P;P [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kP);
  histos->UserHistogram("Track","PIn","PIn;PIn [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track","Eta_phi","Eta vs Phi;Eta;Phi",90,-0.9,0.9,160,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Pt_phi","Pt vs Phi;Pt;Phi [GeV];#tracks",500,0.,5.,320,0.,6.4,AliDielectronVarManager::kPt,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","ImpParXY","ImpParXY; ImpParXY ;#tracks",100,-5.,5.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","ImpParZ","ImpParZ; ImpParZ ;#tracks",100,-5.,5.,AliDielectronVarManager::kImpactParZ);
 
  histos->UserHistogram("Track","NClusterTPC","NClusterTPC; NClusterTPC ;#tracks",200,-0.5,199.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","CrossedRows","CrossedRows; CrossedRows ;#tracks",200,-0.5,199.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","CrossedRowsOverFindable","CrRowsOverFindable; CrRows/FindableCls ;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCchi2perCls","TPCchi2perCls; TPCchi2perCls ;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  
  histos->UserHistogram("Track","NClusterITS","NClusterITS; NClusterITS ;#tracks",8,-0.5,7.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","ITSchi2perCls","ITSchi2perCls; ITSchi2perCls ;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","NSharedClusterITS","NSharedClusterITS; NSharedClusterITS ;#tracks",8,-0.5,7.5,AliDielectronVarManager::kNclsSITS);
  
//   histos->UserHistogram("Track","ITSdEdx_P","dEdx;P [GeV];ITS signal (arb units) vs Momentum;Mom;ITSsignal",     200,0.,10.,150,  0.,150. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
//   histos->UserHistogram("Track","TPCdEdx_P","dEdx;P [GeV];TPC signal (arb units) vs Momentum;Mom;TPCsignal",     200,0.,10.,150,  0.,150. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
//   histos->UserHistogram("Track","TOFbeta_Mom","kTOFbeta vs Momentum;Mom;TOFbeta"                           ,     200,0.,10.,120,  0.,  1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);

   histos->UserHistogram("Track","TPCnSigma_MomEle","TPC number of sigmas Electrons vs Momentum;p;TPCsigmaEle", 300,0.,6.,200,-10., 10. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
   histos->UserHistogram("Track","ITSnSigma_MomEle","ITS number of sigmas Electrons vs Momentum;p;ITSsigmaEle", 300,0.,6.,200,-10., 10. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
   histos->UserHistogram("Track","TOFnSigma_MomEle","TOF number of sigmas Electrons vs Momentum;p;TOFsigmaEle", 300,0.,6.,200,-10., 10. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
 
//   histos->UserHistogram("Track","TPCnSigma_EtaEle","TPC number of sigmas Electrons vs Eta;Eta;TPCsigmaEle", 200,-2.,2.,400,-20., 20. ,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
//   histos->UserHistogram("Track","ITSnSigma_EtaEle","ITS number of sigmas Electrons vs Eta;Eta;ITSsigmaEle", 200,-2.,2.,400,-20., 20. ,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
//   histos->UserHistogram("Track","TOFnSigma_EtaEle","TOF number of sigmas Electrons vs Eta;Eta;TOFsigmaEle", 200,-2.,2.,400,-20., 20. ,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle);
 
//   histos->UserHistogram("Track","TPCnSigma_PhiEle","TPC number of sigmas Electrons vs Phi;Phi;TPCsigmaEle", 200,0.,7.,400,-20., 20. ,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle);
//   histos->UserHistogram("Track","ITSnSigma_PhiEle","ITS number of sigmas Electrons vs Phi;Phi;ITSsigmaEle", 200,0.,7.,400,-20., 20. ,AliDielectronVarManager::kPhi,AliDielectronVarManager::kITSnSigmaEle);
//   histos->UserHistogram("Track","TOFnSigma_PhiEle","TOF number of sigmas Electrons vs Phi;Phi;TOFsigmaEle", 200,0.,7.,400,-20., 20. ,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTOFnSigmaEle);

  //add histograms to pair classes
    TVectorD *binsMee  = AliDielectronHelper::MakeLinBinning(400, 0., 4.);
  TVectorD *binsPtee = AliDielectronHelper::MakeArbitraryBinning("
  0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
  0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,
  0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,
  0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,
  0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,
  0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,
  1.00,1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,
  1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,
  2.00,2.05,2.10,2.15,2.20,2.25,2.30,2.35,2.40,2.45,
  2.50,2.60,2.70,2.80,2.90,3.00,3.10,3.20,3.30,3.40,
  3.50,3.60,3.70,3.80,3.90,4.00,4.10,4.20,4.30,4.40,
  4.50,5.00,5.50,6.00,6.50,7.00
  ");
  TVectorD *binsPhiV = AliDielectronHelper::MakeLinBinning( 30, 0., TMath::Pi());

  histos->UserHistogram("Pair","InvMass_PairPt_PhivPair","InvMass:PairPt:PhivPair;Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                        binsMee,binsPtee,binsPhiV,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);  
//  histos->UserHistogram("Pair","InvMass_PairPt_PhivPair","InvMass:PairPt:PhivPair;Inv. Mass [GeV];Pair Pt [GeV];PhiV",
//                        600,0.,6., 600,0.,6., 30,0.,TMath::Pi(),
//                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair); 
//  histos->UserHistogram("Pair",
//                        "InvMass_pPt","Inv.Mass_PairPt;m_{ee} (GeV/c^{2});p_{T,pair} (GeV/c)",
//                        500,0.,5.,250,0.,5.,
//                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
//  histos->UserHistogram("Pair",
//                        "Eta_phi_pair","Eta vs Phi (pair);Eta;Phi",
//                        50,-1.,1.,80,0.,6.4,
//                        AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
//  histos->UserHistogram("Pair",
//                        "InvMass_PhivPair","InvMass_PhivPair;InvMass;PhivPair",
//                         50, 0. , 0.5, 160 , 0., 3.2,
//                         AliDielectronVarManager::kM , AliDielectronVarManager::kPhivPair );
//  histos->UserHistogram("Pair",
//		                   	"InvMass_OpAngle","InvMass_OpAngle;Invariant Mass;Opening angle",
//		                   	100, 0., 0.5 ,160, 0., 3.2,
//		                   	AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair",
                        "Y","Y;counts;Y",
                        60, -1.2 , 1.2, 
                        AliDielectronVarManager::kY);
                        
//  histos->UserHistogram("Pair","DCA_lin_sqr","#it{DCA}_{ee,lin} vs. #it{DCA}_{ee,sqr}",100,0.,10.,100,0.,10., AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPairLinDCAsigXY);
//
//  //3d histos for pair dca
//  histos->UserHistogram("Pair","InvMass_DCAsigma_pPt","",
//                        100,0.,4, 100,0.,20., 100,0.,5.,
//                        AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPt);
//
//  histos->UserHistogram("Pair","InvMass_DCAsigmaL_pPt","",
//                        100,0.,4, 100,0.,20., 100,0.,5.,
//                        AliDielectronVarManager::kM, AliDielectronVarManager::kPairLinDCAsigXY, AliDielectronVarManager::kPt);
//
//  histos->UserHistogram("Pair","InvMass_DCAsigma_DCAsigmaL","",
//                        100,0.,4, 100,0.,20., 100,0.,20.,
//                        AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPairLinDCAsigXY);
//                        
  if(!isRandomRejTask && cutDefinition < nPF){ 
    // Legs of rejected Pairs. Both charges together. One track can and will make multiple entries.
    histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(1))); // not TrackClassName, see 'AliDielectron::FillHistogramsPair(...)'
    //Create Classes for Rejected Tracks/Pairs:
//    histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(1)));
    histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(0))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
    histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(1)));
    histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(2)));
    
//    histos->UserHistogram("RejPair",
//                          "InvMass_pPt","Inv.Mass_PairPt;m_{ee} (GeV/c^{2});p_{T,pair} (GeV/c)",
//                          500,0.,5.,250,0.,5.,
//                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
//    histos->UserHistogram("RejPair",
//                          "OpAngle_InvMass","InvMass_openingAngle;Invariant Mass;opening angle",
//                          100, 0., 0.2, 100, 0. ,0.2,
//                          AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
//    histos->UserHistogram("RejPair",
//                          "InvMass_PhivPair","InvMass_PhivPair;InvMass;PhivPair",
//                          50, 0. , 0.5, 160 , 0., 3.2,
//                          AliDielectronVarManager::kM , AliDielectronVarManager::kPhivPair );
//    histos->UserHistogram("RejPair",
//                          "Y","Y;counts;Y",
//                          60, -1.2 , 1.2, 
//                          AliDielectronVarManager::kY);                                
    histos->UserHistogram("RejTrack","Pt",";Pt [GeV];#tracks",500,0,10.,AliDielectronVarManager::kPt); 
//    histos->UserHistogram("RejTrack","Pt_charge",";Pt [GeV];#tracks",3,-1.5,1.5,500,0,10.,AliDielectronVarManager::kCharge,AliDielectronVarManager::kPt);
//    histos->UserHistogram("RejTrack","Eta_phi","Eta vs Phi;Eta;Phi",90,-0.9,0.9,160,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
//    histos->UserHistogram("RejTrack","TPCnSigma_MomEle","TPC number of sigmas Electrons vs Momentum;Mom;TPCsigmaEle", 200,0.,10.,300,-30., 30. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
//    
    histos->UserHistogram("Track_Legs","Pt",";Pt [GeV];#tracks",500,0,10.,AliDielectronVarManager::kPt);
//    histos->UserHistogram("Track_Legs","Pt_charge","",3,-1.5,1.5,500,0,10.,AliDielectronVarManager::kCharge,AliDielectronVarManager::kPt);
//    histos->UserHistogram("Track_Legs","Eta_phi","Eta vs Phi;Eta;Phi",90,-0.9,0.9,160,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
//    histos->UserHistogram("Track_Legs","TPCnSigma_MomEle","TPC number of sigmas Electrons vs Momentum;Mom;TPCsigmaEle", 200,0.,10.,300,-30., 30. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);

  }
  die->SetHistogramManager(histos);

}

const AliDielectronEventCuts *GetEventCuts(){

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1); 
  
  return eventCuts;
}