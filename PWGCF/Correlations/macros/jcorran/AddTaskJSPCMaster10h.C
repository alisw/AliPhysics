AliAnalysisTask *AddTaskJSPCMaster10h(TString taskName = "JSPCMaster10h", double ptMin = 0.5, Bool_t removebadarea = kFALSE)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  //-------- JFlucWagons -------
  const int Nsets  = 1; // number of configurations
  TString configNames[Nsets] = {
    "tpconly"//, 
    //"hybrid", // 1
    //"V0M",    // 2
    //"zvtx6",  // 3
    //"zvtx12"  // 4
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
   
    fJCatalyst[i]->SetSaveAllQA(kTRUE); 
    fJCatalyst[i]->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
    fJCatalyst[i]->SetInitializeCentralityArray();

    if (strcmp(configNames[i].Data(), "V0M") == 0) {    
      fJCatalyst[i]->SetCentDetName("V0M"); // V0M for syst
    } else {
      fJCatalyst[i]->SetCentDetName("CL1"); // SPD clusters in the default analysis.
    }
    if (strcmp(configNames[i].Data(), "hybrid") == 0) {
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
  const int SPCCombination = 3;
  TString SPC[SPCCombination] = { "2SPC", "3SPC", "4SPC" };
  AliJSPCTask *myTask[Nsets][SPCCombination];

  for (Int_t i = 0; i < Nsets; i++){
  
   for(Int_t j = 0; j < SPCCombination; j++){

	myTask[i][j] = new AliJSPCTask(Form("%s_s_%s_%s", taskName.Data(), configNames[i].Data(), SPC[j].Data()));
  myTask[i][j]->SetJCatalystTaskName(fJCatalyst[i]->GetJCatalystTaskName());

  myTask[i][j]->JSPCSetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
  myTask[i][j]->JSPCSetSaveAllQA(kTRUE);
  myTask[i][j]->JSPCSetMinNuPar(14.);
  myTask[i][j]->JSPCSetFisherYates(kFALSE, 1.); 
  myTask[i][j]->JSPCSetUseWeights(kFALSE);

	if(j==0){
	  myTask[i][j]->JSPCSetCorrSet1(4., 6.,-2.,-2.,-2., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet2(3., 6.,-3.,-3., 0., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet3(3., 4.,-2.,-2., 0., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet4(3., 8.,-4.,-4., 0., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet5(5., 3., 3.,-2.,-2.,-2., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet6(5., 8.,-2.,-2.,-2.,-2., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet7(0., 0., 0., 0., 0., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet8(0., 0., 0., 0., 0., 0., 0.,0.);
	}

	if(j==1){
	  myTask[i][j]->JSPCSetCorrSet1(4., 3.,-4.,-4., 5., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet2(3., 2., 4.,-6., 0., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet3(3., 2., 3.,-5., 0., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet4(4., 2.,-3.,-3., 4., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet5(5., 2., 3., 3.,-4.,-4., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet6(6., 2., 2., 2., 2.,-3.,-5.,0.);
	  myTask[i][j]->JSPCSetCorrSet7(5., 3., 3., 3.,-4.,-5., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet8(0., 0., 0., 0., 0., 0., 0.,0.);
	}

	if(j==2){
	  myTask[i][j]->JSPCSetCorrSet1(4., 2.,-3.,-4., 5., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet2(6., 2., 2., 2., 3.,-4.,-5.,0.);
	  myTask[i][j]->JSPCSetCorrSet3(5., 2., 2.,-3., 4.,-5., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet4(6., 2., 2., 3., 3.,-4.,-6.,0.);
	  myTask[i][j]->JSPCSetCorrSet5(0., 0., 0., 0., 0., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet6(0., 0., 0., 0., 0., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet7(0., 0., 0., 0., 0., 0., 0.,0.);
	  myTask[i][j]->JSPCSetCorrSet8(0., 0., 0., 0., 0., 0., 0.,0.);
	}

  myTask[i][j]->JSPCSetMixed(kFALSE,2., kFALSE, kTRUE);


	mgr->AddTask((AliAnalysisTask *) myTask[i][j]);
   } 
  }

// Connect the input and output.
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *jHist[Nsets][SPCCombination+1];

  for (Int_t i = 0; i < Nsets; i++) {
    mgr->ConnectInput(fJCatalyst[i], 0, cinput);

    for(Int_t j = 0; j < SPCCombination; j++){
	    mgr->ConnectInput(myTask[i][j], 0, cinput);
      jHist[i][j] = new AliAnalysisDataContainer();     
      jHist[i][j] = mgr->CreateContainer(Form ("%s", myTask[i][j]->GetName()),
        TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
      mgr->ConnectOutput(myTask[i][j], 1, jHist[i][j]);
    }

      jHist[i][SPCCombination] = new AliAnalysisDataContainer();
    jHist[i][SPCCombination] = mgr->CreateContainer(Form ("%s", fJCatalyst[i]->GetName()),
        TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(fJCatalyst[i], 1, jHist[i][SPCCombination]);
  }

  return myTask[0][0];
}

