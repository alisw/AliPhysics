TString names=(
               "PID1_SPDfirst_PrefAllm60t100_pt200;"
               "PID1_SPDorSDD_PrefAllm60t100_pt200;"
               "PID1_SPDfirst_PrefAllm60t100_pt400;"
               "PID1_SPDorSDD_PrefAllm60t100_pt400;"
               );
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();
//
// _____ for Random Rejection task (AddTask_reichelt_RandomRejection.C)
//TF1    PtFunc = TF1("MyPtFunc", "exp(-x/3.)", 0.4, 3.5); // not working
const TString  RndmPtExpr="exp(-x/3.)";
const Double_t RndmPtMin=0.2; // pt and eta ranges need to cover at least the kinematic range of final analysis electrons
const Double_t RndmPtMax=3.5;
const Double_t RndmEtaMax=0.8;
const Int_t nTestpartPerEle=100; // number of testparticles used per final analysis electron in an event.
// _____
const Bool_t randomizeDau=kFALSE;

AliDielectron* Config_reichelt_LMEEPbPb2011(Int_t cutSet, Bool_t hasMC=kFALSE, Bool_t isESD=kFALSE, Bool_t isRandomRej=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  
  LMEECutLib*  LMcutlib = new LMEECutLib();
  LMcutlib->SetIsESDTask(isESD);
  LMcutlib->SetIsRandomRejTask(isRandomRej);
  
  // task name
  TString name=Form("%02d",cutSet);
  if (cutSet<arrNames->GetEntriesFast())  name=arrNames->At(cutSet)->GetName();
  
  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("AliDielectron with cuts: %s",name.Data()));
  die->SetHasMC(hasMC);
  
  //
  // Setup Analysis Selection
  // __________________________________________________
  // SOME POSSIBLE SETTINGS
  // shown together with useful or mandatory cutlib setting
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // deactivate pairing to check track cuts or run with loose pid cuts:
  //   die->SetNoPairing();
  //   LMcutlib->SetIsQATask(kFALSE);
  // apply correct Pre-Filter scheme, if necessary:
  //   die->SetPreFilterAllSigns();   LMcutlib->SetDoRejectionStep(kTRUE);
  //   die->SetPreFilterUnlikeOnly(); LMcutlib->SetDoRejectionStep(kTRUE);
  //   (prefilter will automatically be deactivated in 'AliAnalysisTaskRandomRejection.cxx')
  // set up internal train for efficient systematics study:
  //  in first cutSet do:
  //   die->SetDontClearArrays(kTRUE); // keep pair arrays for next cutSet(s). they will not be modified by those cutSets!
  //  in subsequent cutSets do:
  //   die->SetEventProcess(kFALSE); // use existing pair arrays and study additional cuts. see AliDielectron::FillHistogramsFromPairArray()
  //  HOWEVER: there must be a bug in the Mixed-Events, where some pair legs cannot be found. Those pair type histograms have much less entries!
  // --------------------------------------------------
  cout << endl;
  cout<<"_________________________ "<<"SET UP cutSet NR: "<<cutSet<<" -----> "<<name.Data()<<" _________________________"<<endl;
  
  // --------------------------------------------------
  // common settings:
  // --------------------------------------------------
  LMcutlib->selectedCentrality  = LMEECutLib::kPbPb2011_10to50;
  // prefilter settings:
  die->SetPreFilterAllSigns();  LMcutlib->SetDoRejectionStep(kTRUE);
  LMcutlib->selectedPIDPre      = LMEECutLib::kPbPb2011PID_TPCITSif_2;
  LMcutlib->selectedQualityPre  = LMEECutLib::kPbPb2011TRK_FilterBit0;
  LMcutlib->selectedKineCutsPre = LMEECutLib::kKineCut_pt50_eta090;
  LMcutlib->selectedPairCutsPre = LMEECutLib::kPairCut_mee60_theta100;
  // ana settings:
  LMcutlib->selectedPairCutsAna = LMEECutLib::kPairCut_theta20;
  // post-PID corrections must be done later, because they depend on kine cuts, which may still change.
  //
  // --------------------------------------------------
  // additional settings:
  // --------------------------------------------------
  switch (cutSet) {
    case 0:
      LMcutlib->selectedPIDAna      = LMEECutLib::kPbPb2011PID_ITSTPCTOFif_1;
      LMcutlib->selectedQualityAna  = LMEECutLib::kPbPb2011TRK_SPDfirst_1;
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt200_eta080;
      break;
    case 1:
      LMcutlib->selectedPIDAna      = LMEECutLib::kPbPb2011PID_ITSTPCTOFif_1;
      LMcutlib->selectedQualityAna  = LMEECutLib::kPbPb2011TRK_SPDorSDD_1;
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt200_eta080;
      break;
    case 2:
      LMcutlib->selectedPIDAna      = LMEECutLib::kPbPb2011PID_ITSTPCTOFif_1;
      LMcutlib->selectedQualityAna  = LMEECutLib::kPbPb2011TRK_SPDfirst_1;
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt400_eta080;
      break;
    case 3:
      LMcutlib->selectedPIDAna      = LMEECutLib::kPbPb2011PID_ITSTPCTOFif_1;
      LMcutlib->selectedQualityAna  = LMEECutLib::kPbPb2011TRK_SPDorSDD_1;
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt400_eta080;
      break;
    default:
      cout << " =============================== " << endl;
      cout << " ==== INVALID CONFIGURATION ==== " << endl;
      cout << " cutSet = " << cutSet << endl;
      cout << " =============================== " << endl;
      return 0x0;
  }
  LMcutlib->SetITSSigmaEleCorrection(die, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kEta);
  LMcutlib->SetTPCSigmaEleCorrection(die, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kEta);
  
  //
  // Further configuration of task
  //
  // add centrality selection to event cuts
  die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts() );
  // switch off KF Particle
  die->SetUseKF(kFALSE);

  // --------------------------------------------------
  // with Rejection Step (Prefilter)
  // --------------------------------------------------
  // in case of internal train, use the logic of "without rejection step",
  // so that the analysis track cuts will be given to the correct filter (not to fPairPreFilterLegs).
  // --------------------------------------------------
  if (LMcutlib->GetDoRejectionStep() && die->DoEventProcess()) 
  {
    //if (isESD) { // now done within CutLib
    //  die->GetTrackFilter().AddCuts( LMcutlib->GetESDTrackCutsAna() );
    //  die->GetPairPreFilterLegs().AddCuts( LMcutlib->GetESDTrackCutsAna() ); // this is redundant!?
    //}
    // set initial track filter.
    die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsPre() );
    
    // set Prefilter. "remove all tracks from the Single track arrays that pass the cuts in this filter" (comment in AliDielectron.cxx)
    // cuts = REJECTION!!!
    die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre() ); 
    
    // "apply leg cuts after the pre filter" (comment in AliDielectron.cxx)
    die->GetPairPreFilterLegs().AddCuts( LMcutlib->GetTrackCutsAna() );
    
    // set Pairfilter.
    // cuts = SELECTION!!!
	  die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna() );
    
	}
  // --------------------------------------------------
  // without Rejection Step
  // --------------------------------------------------
	else 
  {
    //if (isESD) { // now done within CutLib
    //  die->GetTrackFilter().AddCuts( LMcutlib->GetESDTrackCutsAna() );
    //}
	  die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna() );
	  die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna() );
	}
  // --------------------------------------------------
  
  
  AliDielectronTrackRotator *rot= 0x0;
  //To save time and as it is not 100% test, rotation switched off
  /*AliDielectronTrackRotator *rot= LMcutlib->GetTrackRotator();
   die->SetTrackRotator(rot);
   */
  AliDielectronMixingHandler *mix=LMcutlib->GetMixingHandler();
  die->SetMixingHandler(mix);
  
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //
  LMcutlib->InitHistograms(die,cutSet);
  
  // the last definition uses no cuts and only the QA histograms should be filled!
  //  LMcutlib->InitCF(die,cutSet);
  
  return die;
}
