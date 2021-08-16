AliAnalysisTask *AddTaskJHOCFAMaster10h (TString taskName = "JHOCFAMaster10h", double ptMin = 0.5, Bool_t removebadarea = kFALSE,
  Bool_t applyHMOcut = kTRUE, Bool_t saveCatalystQA = kFALSE, Bool_t saveHMOQA = kFALSE, Int_t Nsets = 6, TString configArray = "0 1 6 7 8 9",
  Int_t Ncombi = 6, TString combiArray = "2 3 4 2 3 5 2 3 6 2 4 5 2 4 6 3 4 5")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

// Configuration of the wagons.
  const int maxNsets = 10;   // Maximum number of configurations (= max number of mapfiles available).
  Int_t iConfig = -1;
  Int_t index = 0;
  TString configNames[maxNsets];
  std::stringstream sConfig;
  sConfig << configArray;

  while (sConfig >> iConfig) {
    if (index >= Nsets) {break;}
    switch(iConfig) { // Hardcoded names to prevent typo in file names.
      case 0 :  // Default: tpconly.
        configNames[index] = "tpconly";
        break;
      case 1 :  // Syst: hybrid.
        configNames[index] = "hybrid";
        break;
      case 2 :  // Syst: nqq. TBI
        configNames[index] = "nqq";
        break;
      case 3 :  // Syst: pileup. TBI
        configNames[index] = "pileup";
        break;
      case 4 :  // Syst: SPD. TBI??
        configNames[index] = "SPD";
        break;
      case 5 :  // Syst: subA. TBI
        configNames[index] = "subA";
        break;
      case 6 :  // Syst: V0M.
        configNames[index] = "V0M";
        break;
      case 7 :  // Syst: ztx < 8.
        configNames[index] = "zvtx";
        break;
      case 8 :  // Syst: ztx < 6.
        configNames[index] = "zvtx6";
        break;
      case 9 :  // Syst: ztx < 12.
        configNames[index] = "zvtx12";
        break;
      default :
        printf("ERROR: Invalid configuration index. Skipping this element.\n");
        Nsets--;  // To ensure the array is still the right length.
        index--;  // To come back to the situation before the mess.
    }
    index++;
  }

// Loading of the correction map.
  TString MAPfilenames[maxNsets];
  TString MAPdirname = "alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
  AliJCorrectionMapTask *cMapTask = new AliJCorrectionMapTask ("JCorrectionMapTask");

  cMapTask->EnableEffCorrection ("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC10h-LHC11a10a_bis-0-Lists.root"); // Efficiency correction.
  for (Int_t i = 0; i < Nsets; i++) {
    MAPfilenames[i] = Form ("%sPhiWeights_LHC10h_Error_pt%02d_s_%s.root", MAPdirname.Data(), Int_t (ptMin * 10), configNames[i].Data());  // Azimuthal correction.
    cMapTask->EnablePhiCorrection (i, MAPfilenames[i]); // i is index for set file correction ->SetPhiCorrectionIndex(i);
  }
  mgr->AddTask((AliAnalysisTask *) cMapTask); // Loaded correction map added to the analysis manager.

// Setting of the general parameters.
  UInt_t configTrigger = AliVEvent::kMB; // Minimum bias trigger for LHC10h.
  Int_t TPConlyFilter = 128;  // Filterbit for TPConly tracks.
  Int_t hybridCut = 768;

// Configuration of the catalyst tasks for the data selection.
  AliJCatalystTask *fJCatalyst[maxNsets];  // One catalyst needed per configuration.
  for (Int_t i = 0; i < Nsets; i++) {
    fJCatalyst[i] = new AliJCatalystTask(Form("JCatalystTask_%s_s_%s", taskName.Data(), configNames[i].Data()));
    cout << "Setting the catalyst: " << fJCatalyst[i]->GetJCatalystTaskName() << endl;

    if (applyHMOcut) {fJCatalyst[i]->AddFlags(AliJCatalystTask::FLUC_CUT_OUTLIERS);}
    fJCatalyst[i]->SelectCollisionCandidates(configTrigger);
    fJCatalyst[i]->SetSaveAllQA(saveCatalystQA);
    fJCatalyst[i]->SetSaveHMOhist(saveHMOQA);
    fJCatalyst[i]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
    fJCatalyst[i]->SetInitializeCentralityArray();

    if (strcmp(configNames[i].Data(), "V0M") == 0) {    
      fJCatalyst[i]->SetCentDetName("V0M"); // V0M for syst
    }
    else {
      fJCatalyst[i]->SetCentDetName("CL1"); // SPD clusters in the default analysis.
    }

    if (strcmp(configNames[i].Data(), "hybrid") == 0) {
      fJCatalyst[i]->SetTestFilterBit(hybridCut);
    }
    else {
      fJCatalyst[i]->SetTestFilterBit(TPConlyFilter); // default
    }

    if (strcmp(configNames[i].Data(), "zvtx6") == 0) {    
      fJCatalyst[i]->SetZVertexCut(6.0);
    }
    else if (strcmp(configNames[i].Data(), "zvtx12") == 0) {    
      fJCatalyst[i]->SetZVertexCut(12.0);
    }
    else if (strcmp(configNames[i].Data(), "zvtx") == 0) {    
      fJCatalyst[i]->SetZVertexCut(8.0);
    } // Else: do nothing, default is 10.

    fJCatalyst[i]->SetPtRange(ptMin, 5.0);
    fJCatalyst[i]->SetEtaRange(-0.8, 0.8);
    fJCatalyst[i]->SetPhiCorrectionIndex (i);
    fJCatalyst[i]->SetRemoveBadArea(removebadarea);

    mgr->AddTask((AliAnalysisTask *)fJCatalyst[i]);
  }

// Configuration of the analysis task itself.
  AliJHOCFATask *myTask[maxNsets];
  for (Int_t i = 0; i < Nsets; i++){
    myTask[i] = new AliJHOCFATask(Form("%s_s_%s", taskName.Data(), configNames[i].Data()));
    myTask[i]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());
    myTask[i]->SetDebugLevel(0);
    myTask[i]->HOCFASetCentralityBinning(9);
    myTask[i]->HOCFASetCentralityArray("0. 5. 10. 20. 30. 40. 50. 60. 70. 80.");
    myTask[i]->HOCFASetMinMultiplicity(10);
    myTask[i]->HOCFASetParticleWeights(kTRUE);
    myTask[i]->HOCFASetNumberCombi(Ncombi);
    myTask[i]->HOCFASetHarmoArray(Form("%s", combiArray.Data())); // 6 combis.
    mgr->AddTask((AliAnalysisTask *) myTask[i]);
  }

// Connect the input and output.
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *jHist[2*maxNsets];

  for (Int_t i = 0; i < Nsets; i++) {
    mgr->ConnectInput(fJCatalyst[i], 0, cinput);
    mgr->ConnectInput(myTask[i], 0, cinput);

    jHist[i] = new AliAnalysisDataContainer();
    jHist[Nsets+i] = new AliAnalysisDataContainer();    
    jHist[i] = mgr->CreateContainer(Form ("%s", myTask[i]->GetName()),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(myTask[i], 1, jHist[i]);

    jHist[Nsets+i] = mgr->CreateContainer(Form ("%s", fJCatalyst[i]->GetName()),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(fJCatalyst[i], 1, jHist[Nsets+i]);
  }

  return myTask[0];
}
