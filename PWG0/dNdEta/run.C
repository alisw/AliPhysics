void Load(const char* taskName, Bool_t debug)
{
  TString compileTaskName;
  compileTaskName.Form("%s.cxx+", taskName);
  if (debug)
    compileTaskName += "g";

  if (gProof) {
    gProof->Load(compileTaskName);
  } else
    gROOT->Macro(compileTaskName);

  // Enable debug printouts
  if (debug)
  {
    AliLog::SetClassDebugLevel(taskName, AliLog::kDebug+2);
  }
  else
    AliLog::SetClassDebugLevel(taskName, AliLog::kWarning);
}

void run(Int_t runWhat, const Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Int_t aProof = kFALSE, Int_t requiredData = 1, const char* option = "")
{
  // runWhat options: 0 = AlidNdEtaTask
  //                  1 = AlidNdEtaCorrectionTask
  //                  2 = both
  //
  // aProof option: 0 no proof
  //                1 proof with chain
  //                2 proof with dataset
  //
  // requiredData option: 0 = only ESD
  //                      1 = ESD+MC
  //                      2 = RAW (ESD+check on event type)
  //
  // option is passed to the task(s)
  //   option SAVE is removed and results in moving the output files to maps/<ds name>/<trigger>/<det>
  //
  
  TString taskName;
  if (runWhat == 0 || runWhat == 2)
  {
    Printf("Running AlidNdEtaTask");
  }
  if (runWhat == 1 || runWhat == 2)
  {
    Printf("Running AlidNdEtaCorrectionTask");
    if (requiredData != 1)
    {
      Printf("AlidNdEtaCorrectionTask needs MC. Exiting...");
      return;
    }
  }

  if (nRuns < 0)
    nRuns = 1234567890;

  if (aProof)
  {
    TProof::Open("alicecaf"); 

    Bool_t fullAliroot = kFALSE;
    // Enable the needed package
    if (1)
    {
      gProof->UploadPackage("$ALICE_ROOT/STEERBase");
      gProof->EnablePackage("$ALICE_ROOT/STEERBase");
      gProof->UploadPackage("$ALICE_ROOT/ESD");
      gProof->EnablePackage("$ALICE_ROOT/ESD");
      gProof->UploadPackage("$ALICE_ROOT/AOD");
      gProof->EnablePackage("$ALICE_ROOT/AOD");
      gProof->UploadPackage("$ALICE_ROOT/ANALYSIS");
      gProof->EnablePackage("$ALICE_ROOT/ANALYSIS");
      gProof->UploadPackage("$ALICE_ROOT/ANALYSISalice");
      gProof->EnablePackage("$ALICE_ROOT/ANALYSISalice");
    }
    else if (!fullAliroot)
    {
      gProof->UploadPackage("$ALICE_ROOT/AF-v4-18-12-AN.par");
      gProof->EnablePackage("AF-v4-18-12-AN");
    }
    else
    {
      // needed if ITS recpoints are accessed, see AlidNdEtaTask, FULLALIROOT define statement
      gProof->UploadPackage("$ALICE_ROOT/v4-18-15-AN-all.par");
      gProof->EnablePackage("v4-18-15-AN-all");
    
      gProof->Exec("TGrid::Connect(\"alien://\")", kTRUE);
      
      // TODO add this to loadlibs.C
      gProof->Exec("gSystem->Load(\"libXMLParser\")", kTRUE);
    }

    gProof->UploadPackage("$ALICE_ROOT/PWG0base");
    gProof->EnablePackage("$ALICE_ROOT/PWG0base");
  }
  else
  {
    gSystem->AddIncludePath("-I${ALICE_ROOT}/include/ -I${ALICE_ROOT}/PWG0/ -I${ALICE_ROOT}/PWG0/dNdEta/"); 
    gSystem->Load("libVMC");
    gSystem->Load("libTree");
    gSystem->Load("libProof");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWG0base");
  }

  // Create the analysis manager
  mgr = new AliAnalysisManager;

  // Add ESD handler
  
  if (fullAliroot)
    AliESDInputHandler* esdH = new AliESDInputHandlerRP; // for RecPoints
  else
    AliESDInputHandler* esdH = new AliESDInputHandler;
  
  esdH->SetInactiveBranches("FMD AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks MuonTracks TrdTracks");
  mgr->SetInputEventHandler(esdH);

  AliPWG0Helper::AnalysisMode analysisMode = AliPWG0Helper::kSPD | AliPWG0Helper::kFieldOn;
  //AliPWG0Helper::AnalysisMode analysisMode = AliPWG0Helper::kSPD | AliPWG0Helper::kFieldOn | AliPWG0Helper::kSPDOnlyL0;
  //AliPWG0Helper::AnalysisMode analysisMode = AliPWG0Helper::kTPCITS | AliPWG0Helper::kFieldOn;
  
  //AliTriggerAnalysis::Trigger trigger      = AliTriggerAnalysis::kAcceptAll | AliTriggerAnalysis::kOfflineFlag;
  AliTriggerAnalysis::Trigger trigger      = AliTriggerAnalysis::kAcceptAll | AliTriggerAnalysis::kOfflineFlag | AliTriggerAnalysis::kOneParticle;
  
  //AliTriggerAnalysis::Trigger trigger      = AliTriggerAnalysis::kSPDGFOBits | AliTriggerAnalysis::kOfflineFlag;
  //AliTriggerAnalysis::Trigger trigger      = AliTriggerAnalysis::kSPDGFOBits | AliTriggerAnalysis::kOfflineFlag | AliTriggerAnalysis::kOneParticle;
  
  //AliTriggerAnalysis::Trigger trigger      = AliTriggerAnalysis::kV0AND | AliTriggerAnalysis::kOfflineFlag; 
  
  //AliTriggerAnalysis::Trigger trigger      = AliTriggerAnalysis::kV0OR | AliTriggerAnalysis::kOfflineFlag; 
  //AliTriggerAnalysis::Trigger trigger      = AliTriggerAnalysis::kV0OR | AliTriggerAnalysis::kOfflineFlag | AliTriggerAnalysis::kOneParticle; 

  AliPWG0Helper::DiffTreatment diffTreatment = AliPWG0Helper::kMCFlags;
  //AliPWG0Helper::DiffTreatment diffTreatment = AliPWG0Helper::kE710Cuts;
  
  AliPWG0Helper::PrintConf(analysisMode, trigger, diffTreatment);

  AliESDtrackCuts* esdTrackCuts = 0;
  if (!(analysisMode & AliPWG0Helper::kSPD))
  {
    // selection of esd tracks
    gROOT->ProcessLine(".L ../CreateStandardCuts.C");
    esdTrackCuts = CreateTrackCuts(analysisMode);
    if (!esdTrackCuts)
    {
      printf("ERROR: esdTrackCuts could not be created\n");
      return;
    }
    esdTrackCuts->SetHistogramsOn(kTRUE);
  }

  cInput = mgr->GetCommonInputContainer();
  
  // remove SAVE option if set
  Bool_t save = kFALSE;
  TString optStr(option);
  if (optStr.Contains("SAVE"))
  {
    optStr = optStr(0,optStr.Index("SAVE")) + optStr(optStr.Index("SAVE")+4, optStr.Length());
    save = kTRUE;
  }
  
  // physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  physicsSelectionTask = AddTaskPhysicsSelection((requiredData == 2) ? kFALSE : kTRUE);
  
  // 900 GeV 
  if (0 && requiredData == 2)
  {
    physicsSelectionTask->GetPhysicsSelection()->AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL #769 #3119");
    physicsSelectionTask->GetPhysicsSelection()->AddBGTriggerClass("+CINT1A-ABCE-NOPF-ALL #446 #2554");
    physicsSelectionTask->GetPhysicsSelection()->AddBGTriggerClass("+CINT1C-ABCE-NOPF-ALL #1334 #2228");
    physicsSelectionTask->GetPhysicsSelection()->AddBGTriggerClass("+CINT1-E-NOPF-ALL #790");
  }
  
  // 7 TeV, run 114783
  if (0 && requiredData == 2)
  {
    physicsSelectionTask->GetPhysicsSelection()->AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL #345");
    physicsSelectionTask->GetPhysicsSelection()->AddBGTriggerClass("+CINT1A-ABCE-NOPF-ALL #2130");
    physicsSelectionTask->GetPhysicsSelection()->AddBGTriggerClass("+CINT1C-ABCE-NOPF-ALL #3018");
    physicsSelectionTask->GetPhysicsSelection()->AddBGTriggerClass("+CINT1-E-NOPF-ALL #1238");
  }

  // 7 TeV, run 114786,98
  if (0 && requiredData == 2)
  {
    physicsSelectionTask->GetPhysicsSelection()->AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL #346");
    physicsSelectionTask->GetPhysicsSelection()->AddBGTriggerClass("+CINT1A-ABCE-NOPF-ALL #2131");
    physicsSelectionTask->GetPhysicsSelection()->AddBGTriggerClass("+CINT1C-ABCE-NOPF-ALL #3019");
    physicsSelectionTask->GetPhysicsSelection()->AddBGTriggerClass("+CINT1-E-NOPF-ALL #1238");
    //physicsSelectionTask->GetPhysicsSelection()->Initialize(114786);
  }

  // FO efficiency (for MC)
  if (1 && requiredData != 2)
  {
    //const char* fastORFile = "spdFOEff_run104824_52.root";
    //const char* fastORFile = "spdFOEff_run104867_92.root";
    //const char* fastORFile = "spdFOEff_run105054_7.root";
    const char* fastORFile = "spdFOEff_run114931.root";
  
    Printf("NOTE: Simulating FAST-OR efficiency on the analysis level using file %s", fastORFile);
    TFile::Open(fastORFile);
    spdFOEff = (TH1F*) gFile->Get("spdFOEff");
    physicsSelectionTask->GetPhysicsSelection()->Initialize(114931);
    physicsSelectionTask->GetPhysicsSelection()->GetTriggerAnalysis()->SetSPDGFOEfficiency(spdFOEff);
  }
  
  // V0 syst. study
  if (0)
  {
    Printf("NOTE: Systematic study for VZERO enabled!");
    physicsSelectionTask->GetPhysicsSelection()->Initialize(104867);
    for (Int_t i=0; i<1; i++)
    {
      // for MC and data
      //physicsSelectionTask->GetPhysicsSelection()->GetTriggerAnalysis(i)->SetV0HwPars(15, 61.5, 86.5);
      physicsSelectionTask->GetPhysicsSelection()->GetTriggerAnalysis(i)->SetV0AdcThr(6);
      // only for MC
      //physicsSelectionTask->GetPhysicsSelection()->GetTriggerAnalysis(i)->SetV0HwPars(0, 0, 125);
      //physicsSelectionTask->GetPhysicsSelection()->GetTriggerAnalysis(i)->SetV0AdcThr(0);
    }
  }
  
  // BG study
  //physicsSelectionTask->GetPhysicsSelection()->AddCollisionTriggerClass("+CINT1A-ABCE-NOPF-ALL");
  //physicsSelectionTask->GetPhysicsSelection()->AddCollisionTriggerClass("+CINT1C-ABCE-NOPF-ALL");
  
  // Create, add task
  if (runWhat == 0 || runWhat == 2)
  {
    Load("AlidNdEtaTask", aDebug);
    task = new AlidNdEtaTask(optStr);

    if (requiredData == 1)
      task->SetReadMC();
      
    //physicsSelectionTask->GetPhysicsSelection()->SetBin0Callback("AlidNdEtaTask");

    // syst. error flags
    //task->SetUseMCVertex();
    //task->SetUseMCKine();
    //task->SetOnlyPrimaries();
    //task->SetFillPhi();
    //task->SetSymmetrize();

    // INEL>0 definition
    if (trigger & AliTriggerAnalysis::kOneParticle)
      task->SetMultAxisEta1();

    task->SetTrigger(trigger);
    task->SetAnalysisMode(analysisMode);
    task->SetTrackCuts(esdTrackCuts);
    //task->SetDeltaPhiCut(0.064);
    task->SetDiffTreatment(diffTreatment);

    mgr->AddTask(task);

    // Attach input
    mgr->ConnectInput(task, 0, cInput);

    // Attach output
    cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer);
    mgr->ConnectOutput(task, 1, cOutput);
  }

  if (runWhat == 1 || runWhat == 2)
  {
    Load("AlidNdEtaCorrectionTask", aDebug);
    task2 = new AlidNdEtaCorrectionTask(optStr);

    // syst. error flags
    //task2->SetFillPhi();
    //task2->SetOnlyPrimaries();
    //task2->SetSymmetrize();

    // to account for gaps in real life SPD geometry
    task2->SetSkipParticles();

    // INEL>0 definition
    if (trigger & AliTriggerAnalysis::kOneParticle)
      task2->SetMultAxisEta1();

    task2->SetTrigger(trigger);
    task2->SetAnalysisMode(analysisMode);
    task2->SetTrackCuts(esdTrackCuts);
    //task2->SetDeltaPhiCut(0.064);
    task2->SetDiffTreatment(diffTreatment);

    mgr->AddTask(task2);

    // Attach input
    mgr->ConnectInput(task2, 0, cInput);

    // Attach output
    cOutput = mgr->CreateContainer("cOutput2", TList::Class(), AliAnalysisManager::kOutputContainer);
    mgr->ConnectOutput(task2, 0, cOutput);
  }

  if (requiredData == 1) 
  {
    // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // Enable debug printouts
  if (aDebug)
    mgr->SetDebugLevel(2);

  // Run analysis
  mgr->InitAnalysis();
  mgr->PrintStatus();

  if (aProof == 2)
  {
    // process dataset

    mgr->StartAnalysis("proof", data, nRuns, offset);
    
    if (save)
    {
      TString path("maps/");
      path += TString(data).Tokenize("/")->Last()->GetName();
      
      UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) AliTriggerAnalysis::kStartOfFlags;
      switch (triggerNoFlags)
      {
        case AliTriggerAnalysis::kMB1: path += "/mb1"; break;
        case AliTriggerAnalysis::kMB2: path += "/mb2"; break;
        case AliTriggerAnalysis::kMB3: path += "/mb3"; break;
        case AliTriggerAnalysis::kSPDGFO: path += "/spdgfo"; break;
        case AliTriggerAnalysis::kSPDGFOBits: path += "/spdgfobits"; break;
        case AliTriggerAnalysis::kAcceptAll: path += "/all"; break;
        case AliTriggerAnalysis::kV0AND: path += "/v0and"; break;
        case AliTriggerAnalysis::kV0OR: path += "/v0or"; break;
        case AliTriggerAnalysis::kNSD1: path += "/nsd1"; break;
        case AliTriggerAnalysis::kMB1Prime: path += "/mb1prime"; break;
        default: Printf("ERROR: Trigger undefined for path to files"); return;
      }
      
      if (trigger & AliTriggerAnalysis::kOneParticle)
        path += "-onepart";
      
      if (strlen(requireClass) > 0 && strlen(rejectClass) == 0)
      {
        path += Form("/%s", requireClass);
      }
      else if (strlen(rejectClass) > 0)
        path += Form("/%s--%s", requireClass, rejectClass);
      
      if (analysisMode & AliPWG0Helper::kSPD)
        path += "/spd";
      
      if (analysisMode & AliPWG0Helper::kSPDOnlyL0)
        path += "onlyL0";
      
      if (analysisMode & AliPWG0Helper::kTPC)
        path += "/tpc";
        
      if (analysisMode & AliPWG0Helper::kTPCITS)
        path += "/tpcits";

      gSystem->mkdir(path, kTRUE);
      if (runWhat == 0 || runWhat == 2)
      {
        gSystem->Rename("analysis_esd_raw.root", path + "/analysis_esd_raw.root");
        if (requiredData == 1)
          gSystem->Rename("analysis_mc.root", path + "/analysis_mc.root");
      }
      if (runWhat == 1 || runWhat == 2)
      {
        if (optStr.Contains("process-types"))
          gSystem->Rename("correction_mapprocess-types.root", path + "/correction_mapprocess-types.root");
        else
          gSystem->Rename("correction_map.root", path + "/correction_map.root");
      }
      gSystem->Rename("event_stat.root", path + "/event_stat.root");
      
      Printf(">>>>> Moved files to %s", path.Data());
    }
  }
  else if (aProof == 3)
  {
    gROOT->ProcessLine(".L CreateChainFromDataSet.C");
    ds = gProof->GetDataSet(data)->GetStagedSubset();
    chain = CreateChainFromDataSet(ds, "esdTree", nRuns);
    mgr->StartAnalysis("local", chain, 1234567890, offset);
  }
  else
  {
    // Create chain of input files
    gROOT->LoadMacro("../CreateESDChain.C");

    chain = CreateESDChain(data, nRuns, offset);
    //chain = CreateChain("TE", data, nRuns, offset);

    mgr->StartAnalysis((aProof > 0) ? "proof" : "local", chain);
  }
}

