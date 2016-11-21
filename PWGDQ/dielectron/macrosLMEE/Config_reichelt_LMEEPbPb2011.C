TString names=(
               "cut16_SPDorSDD14_PID16_pt200;"
               "cut16_ITScuts_pt200;"
               "cut16_ITScuts_Pref_pt200;"
               "cut16_SPDorSDD14_PID16_pt400;"
               "cut16_ITScuts_pt400;"
               "cut16_ITScuts_Pref_pt400;"
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
const Int_t nTestpartPerEle=10; // number of testparticles used per final analysis electron in an event.
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
  LMcutlib->selectedPIDPre      = LMEECutLib::kPbPb2011PID_TPCITSif_2;
  LMcutlib->selectedQualityPre  = LMEECutLib::kPbPb2011TRK_FilterBit0;
  LMcutlib->selectedKineCutsPre = LMEECutLib::kKineCut_pt50_eta090;
  LMcutlib->selectedPairCutsPre = LMEECutLib::kPairCut_mee40_theta80;
  // ana settings:
  //LMcutlib->selectedPairCutsAna = LMEECutLib::kPairCut_theta50; // cut can now be realized via a 4D histogram.
  // post-PID corrections must be done later, because they depend on kine cuts, which may still change.
  // --------------------------------------------------
  // settings for THIS CUT VARIATION:
  // --------------------------------------------------
  LMcutlib->selectedPIDAna      = LMEECutLib::kCut16;
  LMcutlib->selectedQualityAna  = LMEECutLib::kCut16;
  //LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt200_eta080;
  //
  // --------------------------------------------------
  // specific settings for each cutset:
  // --------------------------------------------------
  switch (cutSet) {
    case 0:
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt200_eta080;
      break;
    case 1:
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt200_eta080;
      LMcutlib->SetUseITScutsESD(kTRUE);
      break;
    case 2:
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt200_eta080;
      die->SetPreFilterAllSigns();  LMcutlib->SetDoRejectionStep(kTRUE);
      LMcutlib->SetUseITScutsESD(kTRUE);
      break;
    case 3:
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt400_eta080;
      break;
    case 4:
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt400_eta080;
      LMcutlib->SetUseITScutsESD(kTRUE);
      break;
    case 5:
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt400_eta080;
      die->SetPreFilterAllSigns();  LMcutlib->SetDoRejectionStep(kTRUE);
      LMcutlib->SetUseITScutsESD(kTRUE);
      break;
    default:
      cout << " =============================== " << endl;
      cout << " ==== INVALID CONFIGURATION ==== " << endl;
      cout << " cutSet = " << cutSet << endl;
      cout << " =============================== " << endl;
      return 0x0;
  }
  
  if (hasMC) {
    LMcutlib->SetITSSigmaEleCorrectionMC(die, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
    LMcutlib->SetTPCSigmaEleCorrectionMC(die, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
  } else {
    LMcutlib->SetITSSigmaEleCorrection(die, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kEta);
    LMcutlib->SetTPCSigmaEleCorrection(die, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kEta);
  }
  
  if (hasMC) {
    LMcutlib->AddMCSignals(die, 2);
  }
  
  //
  // Further configuration of task
  //
  // add centrality selection to event cuts
  die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(hasMC) );
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
	  die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna() );
	  die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna() );
	}
  // --------------------------------------------------
  
  if (isRandomRej && !(LMcutlib->GetDoRejectionStep())) return 0x0; // avoid adding unneeded cutsets to random rejection task.
  
  if (!hasMC) { // mixing crashes for MC task
    AliDielectronMixingHandler *mix=LMcutlib->GetMixingHandler();
    die->SetMixingHandler(mix);
  }
  
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //
  LMcutlib->InitHistograms(die,cutSet);
  
  return die;
}
