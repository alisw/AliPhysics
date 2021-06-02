AliAnalysisTask *AddTaskJHOCFAMaster10h (TString taskName = "JHOCFAMaster10h", double ptMin = 0.5, Bool_t removebadarea = kFALSE,
  Bool_t saveCatalystQA = kFALSE)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

// Configuration of the wagons.
  const int Nsets = 7;  // Number of configurations.
  TString configNames[Nsets] = {
    "tpconly",  // 0 = Default.
    "hybrid",   // 1: filtering.
    "V0M",      // 2: centrality.
    "zvtx6",    // 3: PVz < 6cm.
    "zvtx12",   // 4: PVz < 12cm.
    "ntpc50",   // 5: N_TPC > 50.
    "ntpc100"   // 6: N_TPC > 100.
  };

// Loading of the correction map.
  TString MAPfilenames[Nsets];
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
  AliJCatalystTask *fJCatalyst[Nsets];  // One catalyst needed per configuration.
  for (Int_t i = 0; i < Nsets; i++) {
    fJCatalyst[i] = new AliJCatalystTask(Form("JCatalystTask_%s_s_%s", taskName.Data(), configNames[i].Data()));
    cout << "Setting the catalyst: " << fJCatalyst[i]->GetJCatalystTaskName() << endl;

  // TBI: Do we need the flag for outliers here?
    fJCatalyst[i]->SelectCollisionCandidates(configTrigger);
    fJCatalyst[i]->SetSaveAllQA(saveCatalystQA); 
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
    } // Else: do nothing, default is 10.

    if (strcmp(configNames[i].Data(), "ntpc50") == 0) {    
      fJCatalyst[i]->SetNumTPCClusters(50);
    }
    else if (strcmp(configNames[i].Data(), "ntpc100") == 0) {    
      fJCatalyst[i]->SetNumTPCClusters(100);
    } // Else: do nothing, default is 10.

    fJCatalyst[i]->SetPtRange(ptMin, 5.0);
    fJCatalyst[i]->SetEtaRange(-0.8, 0.8);
    fJCatalyst[i]->SetPhiCorrectionIndex (i);

    mgr->AddTask((AliAnalysisTask *)fJCatalyst[i]);
  }

// Configuration of the analysis task itself.
  AliJHOCFATask *myTask[Nsets];
  for (Int_t i = 0; i < Nsets; i++){
    myTask[i] = new AliJHOCFATask(Form("%s_s_%s", taskName.Data(), configNames[i].Data()));
    myTask[i]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());
    myTask[i]->HOCFASetCentralityBinning(9);
    myTask[i]->HOCFASetCentralityArray("0. 5. 10. 20. 30. 40. 50. 60. 70. 80.");
    myTask[i]->HOCFASetMinMultiplicity(6);
    myTask[i]->HOCFASetParticleWeights(kTRUE);
    myTask[i]->HOCFASetNumberCombi(1);
    myTask[i]->HOCFASetHarmoArray("2 3 4");
    mgr->AddTask((AliAnalysisTask *) myTask[i]);
  }

// Connect the input and output.
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *jHist[2*Nsets];

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
