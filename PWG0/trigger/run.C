void Load(const char* taskName, Bool_t debug)
{
  TString compileTaskName;
  compileTaskName.Form("%s.cxx++", taskName);
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

void GetTimes(UInt_t run, UInt_t* startTime = 0, UInt_t* endTime = 0)
{
  gSystem->Load("libXMLParser");
  gSystem->Load("libXMLIO");
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("alien://folder=/alice/data/2009/OCDB");
  //man->SetRun(run);
  AliCDBPath cdb("GRP", "GRP", "Data");
  obj = man->Get(cdb, run);
  grp = (AliGRPObject*) obj->GetObject();
  
  if (startTime)
    *startTime = grp->GetTimeStart();
  if (endTime)
    *endTime = grp->GetTimeEnd();
  
  Printf("Got start and endtime from OCDB: %d, %d", grp->GetTimeStart(), grp->GetTimeEnd());
}

void run(const Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Int_t aProof = kFALSE, UInt_t startTime = 0, UInt_t endTime = 0, const char* option = "")
{
  // aProof option: 0 no proof
  //                1 proof with chain
  //                2 proof with dataset
  //
  // option is passed to the task(s)

/*
 .x run.C("/PWG0/jgrosseo/ERP_run98097", -1, 0, kFALSE, 2, 1258045012, 1258045458)
 .x run.C("/PWG0/jgrosseo/ERP_run98576", -1, 0, kFALSE, 2, 1258123911, 1258124103)
 .x run.C("/PWG0/jgrosseo/ERP_run98569", -1, 0, kFALSE, 2, 1258122187, 1258122524)
 .x run.C("/PWG0/jgrosseo/run101235", -1, 0, kFALSE, 2, 1258821541, 1258822595) 
 .x run.C("/PWG0/jgrosseo/bgRuns_101424to102689_allITS_2", -1, 0, kFALSE, 2, 1258950519, 1259567587) 
 .x run.C("/PWG0/jgrosseo/run000102556", -1, 0, kFALSE, 2, 1259481473, 1259494364) 
 
  timestamps:
 .x run.C("/ALIREC/aliprod/run101498", -1, 0, kFALSE, 2, 1258990726, 1258993311) 
  orbits:
 .x run.C("/ALIREC/aliprod/run101498", -1, 0, kFALSE, 2, 13587, 16749493) 
 
 .x run.C("/ALIREC/aliprod/run104044", -1, 0, 1, 2, 1260078397, 1260082975) 
 
 .x run.C("/ALIREC/aliprod/run104065", -1, 0, 1, 2, 1260084071, 1260084250)
 .x run.C("/ALIREC/aliprod/run104065", -1, 0, 1, 2, 1353, 1504)

 .x run.C("/ALIREC/aliprod/run104321", -1, 0, 1, 2, 1260237771, 1260249377)

*/
  
  if (nRuns < 0)
    nRuns = 1234567890;

  if (aProof)
  {
    TProof::Mgr("alicecaf")->SetROOTVersion("v5-24-00a"); 
    TProof::Open("alicecaf"); 
    //gProof->SetParallel(1);

    // Enable the needed package
    if (0)
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
    else
    {
      gProof->UploadPackage("$ALICE_ROOT/AF-v4-18-12-AN.par");
      gProof->EnablePackage("AF-v4-18-12-AN");
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
  
  if (startTime == endTime && startTime > 0)
  {
    // get times from OCDB, startTime must be run number

    // WARNING only works if your par files loaded above are compatible with the libraries loaded here...
    GetTimes(startTime, &startTime, &endTime);
  }

  // Create the analysis manager
  mgr = new AliAnalysisManager;

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  esdH->SetInactiveBranches("AliESDACORDE ALIESDTZERO AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks Kinks Cascades AliESDTZERO ALIESDACORDE MuonTracks TrdTracks CaloClusters");
  mgr->SetInputEventHandler(esdH);

  cInput = mgr->GetCommonInputContainer();
  
  Load("AliTriggerTask", aDebug);
  TString optStr(option);
  
  // remove SAVE option if set
  Bool_t save = kFALSE;
  if (optStr.Contains("SAVE"))
  {
    optStr = optStr(0,optStr.Index("SAVE")) + optStr(optStr.Index("SAVE")+4, optStr.Length());
    save = kTRUE;
  }
  
  task = new AliTriggerTask(optStr);
  task->SetTimes(startTime, endTime);
  //task->SetUseOrbits(kTRUE);

  mgr->AddTask(task);

  // Attach input
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer);
  mgr->ConnectOutput(task, 0, cOutput);

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
      
      gSystem->mkdir(path, kTRUE);
      gSystem->Rename("trigger.root", path + "/trigger.root");
      
      Printf(">>>>> Moved files to %s", path.Data());
    }  
  }
  else if (aProof == 3)
  {
    gROOT->ProcessLine(".L CreateChainFromDataSet.C");
    ds = gProof->GetDataSet(data)->GetStagedSubset();
    chain = CreateChainFromDataSet(ds);
    mgr->StartAnalysis("local", chain, nRuns, offset);
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
