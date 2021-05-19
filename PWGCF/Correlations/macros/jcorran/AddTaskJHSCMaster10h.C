AliAnalysisTask *AddTaskJHSCMaster10h(TString taskName = "JHSCMaster10h", double ptMin = 0.5, Bool_t removebadarea = kFALSE)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  //-------- JFlucWagons -------
  const int Nsets  = 5; // number of configurations
  TString configNames[Nsets] = {
    "tpconly", 
    "hybrid", // 1
    "V0M",    // 2
    "zvtx6",  // 3
    "zvtx12"  // 4
    //"nqq",    //  not yet implemented in JHSC code
    //"subA",   //not yet implemented in JHSC code
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
  Int_t tpconlyCut = 128;
  Int_t hybridCut = 768;

  UInt_t selEvt;
  selEvt = AliVEvent::kMB;// Minimum bias trigger for LHC10h.

  AliJCatalystTask *fJCatalyst[Nsets];  // One catalyst needed per configuration.
  for (Int_t i = 0; i < Nsets; i++) {
    fJCatalyst[i] = new AliJCatalystTask(Form("fJCatalystTask_%s", configNames[i].Data()));
    cout << "Setting the catalyst: " << fJCatalyst[i]->GetJCatalystTaskName() << endl;
    fJCatalyst[i]->SelectCollisionCandidates(selEvt);
    if (i == 2) {    
      fJCatalyst[i]->SetCentDetName("VOM"); // V0M for syst
    } else {
      fJCatalyst[i]->SetCentDetName("CL1"); // SPD clusters in the default analysis.
    }
    if (i == 1) {
      fJCatalyst[i]->SetTestFilterBit(hybridCut);
    } else {
      fJCatalyst[i]->SetTestFilterBit(tpconlyCut); // default
    }
    fJCatalyst[i]->SetPtRange(ptMin, 5.0);
    fJCatalyst[i]->SetEtaRange(-0.8, 0.8);
    fJCatalyst[i]->SetPhiCorrectionIndex (i);
    mgr->AddTask((AliAnalysisTask *)fJCatalyst[i]);
  }

// Configuration of the analysis task itself.
  AliJHSCTask *myTask[Nsets];
  for (Int_t i = 0; i < Nsets; i++){
    myTask[i] = new AliJHSCTask(Form("%s_s_%s", taskName.Data(), configNames[i].Data()));
    myTask[i]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());
    mgr->AddTask((AliAnalysisTask *) myTask[i]);
  }

// Connect the input and output.
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  for (Int_t i = 0; i < Nsets; i++) {
    mgr->ConnectInput(fJCatalyst[i], 0, cinput);
    mgr->ConnectInput(myTask[i], 0, cinput);
    AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form ("%scontainer", myTask[i]->GetName()),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", AliAnalysisManager::GetCommonFileName(),myTask[i]->GetName()));
    mgr->ConnectOutput(myTask[i], 1, jHist);
  }

  return myTask[0];
}
