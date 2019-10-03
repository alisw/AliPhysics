//-----------------Global variables for the the analysis------------------------
TString aliroot_version = "v4-21-11-AN";
TString analysis_mode = "";
TString output_file = "PWG3hfe.root";
TString root_version = "";
TString api_version = "V1.1x";
TString input_data = "";

//___________________________________________________________
Bool_t FindDependentPackages(){
  // Find packages where aliroot depends on
  // Method is the following: the dependencies are on a webpage on alimonitor
  // we download file and try to find the aliroot expression in the file. The line
  // we read into a string. The we tokenize the string and look for the expression VO_ALICE@ROOT::
  // Then we have all information we need
  const char *wd = gSystem->pwd(); 
  gSystem->cd("/tmp");
  gSystem->Exec("wget http://alimonitor.cern.ch/packages");
  TString setting = gSystem->GetFromPipe(Form("cat index.html | grep VO_ALICE@AliRoot::%s", aliroot_version.Data()));
  gSystem->Exec("rm index.html");
  gSystem->cd(wd);
  if(!setting.Length()) return kFALSE;
  TObjArray *tokens = setting.Tokenize(";");
  TIter tokiter(tokens);
  TObjString *os = NULL;
  while((os = dynamic_cast<TObjString *>(tokiter()))){
    TString &mys = os->String();
    if(!mys.Contains("VO_ALICE@ROOT")) continue;
    mys.ReplaceAll("VO_ALICE@ROOT::", "");
    mys.ReplaceAll("&quot","");
    root_version = mys;
    break;
  }
  // print results:
  printf("Found Packages for Analysis:\n");
  printf("=====================================\n");
  printf("ALIROOT: %s\n", aliroot_version.Data());
  printf("ROOT:    %s\n", root_version.Data());
  printf("\n");
  return kTRUE;
}

//___________________________________________________________
AliAnalysisAlien *CreateAlienHandler(Bool_t isProof){
  //if(!AliAnalysisGrid::CreateToken()) return NULL;
  if(!FindDependentPackages()) return NULL;
  AliAnalysisAlien *alienplugin = new AliAnalysisAlien;
   
  // common settings
  alienplugin->SetAliROOTVersion(aliroot_version);
  alienplugin->SetROOTVersion(root_version);
  alienplugin->SetAPIVersion(api_version);
  alienplugin->SetAdditionalLibs("libPWGHFhfe.so");
  alienplugin->SetDefaultOutputs(kTRUE);
  if(isProof){
    // proof mode
    if(analysis_mode.Contains("test")) alienplugin->SetRunMode("test");
    alienplugin->SetProofCluster("alice_caf.cern.ch");
    alienplugin->SetRootVersionForProof(root_version);
    alienplugin->SetAliRootMode("aliroot");
  } else {
    // grid mode
    alienplugin->SetRunMode(analysis_mode.Data());
    // default setting that need no deeper logic
    alienplugin->SetJDLName("hfeanalysis.jdl");
    alienplugin->SetExecutable("hfeanalysis.sh");
    alienplugin->SetAnalysisMacro("hfeanalysis.C");
    alienplugin->SetValidationScript("hfevalidate.sh");
    alienplugin->SetGridWorkingDir("hfeanalysis");
    /*alienplugin->SetOutputArchive();
    alienplugin->SetDefaultOutputs();*/
    alienplugin->SetOverwriteMode();
    alienplugin->SetFastReadOption();
    alienplugin->SetSplitMaxInputFileNumber(5);
    alienplugin->SetTTL(30000);
    alienplugin->SetInputFormat("xml-single");
    alienplugin->SetPrice(1);
    alienplugin->SetSplitMode("se");
    alienplugin->SetNtestFiles(1);

    // Merging setting, only needed in terminate or full mode
    if(analysis_mode.Contains("full") || analysis_mode.Contains("terminate")){
      alienplugin->SetMergeViaJDL();
      alienplugin->SetMaxMergeFiles(50);
      alienplugin->SetMaxMergeStages(5);
    }
  }
  return alienplugin;
}

//___________________________________________________________
void DecodeDataString(const TString &datastring, TString &sample, TArrayI &listofruns){
  TObjArray *toks = datastring.Tokenize(":");
  sample = (dynamic_cast<TObjString *>(toks->At(0)))->String();
  TString &listrunstring = (dynamic_cast<TObjString *>(toks->At(1)))->String();
  TObjArray *runstrings = listrunstring.Tokenize(",");
  TIter runiter(runstrings);
  listofruns.Set(runstrings->GetEntriesFast());
  TObjString *myrunstring = NULL;
  Int_t counter = 0;
  while((myrunstring = dynamic_cast<TObjString *>(runiter()))) listofruns[counter++] = myrunstring->String().Atoi();  
  // Print summary:
  printf("Selected sample: %s\n", sample.Data());
  printf("========================================\n");
  for(Int_t irun = 0; irun < listofruns.GetSize(); irun++){
    printf("\trun %d\n", listofruns[irun]);
  }
  printf("\n");
  delete toks; delete runstrings;
}

//___________________________________________________________
Int_t GetYear(TString &sample){
  TString yearstring = sample(3,4); 
  Int_t year = yearstring.Atoi();
  return 2000 + year;
}

//___________________________________________________________
void AddInput(AliAnalysisAlien *alienhandler, TString sample, TArrayI &listofruns, Bool_t MC){
  Int_t year = GetYear(sample);
  TString specialpattern, datapattern;
  if(!MC){
    // handle Data
    datapattern = "ESDs/pass2/*ESDs.root";
    specialpattern = Form("data/%d", year);
    //alienhandler->SetRunPrefix("%09d");
  } else {
    // handle MC
    datapattern = "*ESDs.root";
    specialpattern = "sim";
  }
  alienhandler->SetGridDataDir(Form("/alice/%s/%s", specialpattern.Data(), sample.Data()));
  alienhandler->SetDataPattern(datapattern.Data());
  for(Int_t ien = 0; ien < listofruns.GetSize(); ien++){
    if(!MC) alienhandler->AddRunNumber(Form("%09d", listofruns[ien]));
    else alienhandler->AddRunNumber(listofruns[ien]);
  }
  alienhandler->SetGridOutputDir(sample.Data());
  alienhandler->SetOutputToRunNo();
}

//___________________________________________________________
void runAlien(TString data, TString mode = "test", Bool_t MC = kFALSE){  
  if(!gSystem->Getenv("ALICE_ROOT")){
    printf("AliRoot has to be initialized\n");  
    return;
  }
  
  // check for valid modes
  const int kModes = 5;
  TString allowed_modes[kModes] = {"proof", "prooftest", "test", "full", "submit"}; 
  Bool_t isValid = kFALSE;
  mode.ToLower();
  for(int imode = 0; imode < kModes; imode++){
    if(!mode.CompareTo(allowed_modes[imode])) isValid = kTRUE;
  }
  if(!isValid){
    printf("invalid analysis mode selected\n");
    return;
  }
  analysis_mode = mode; 
  Bool_t proofmode = mode.Contains("proof");
  // libraries to be loaded
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGhfe");
  
  // Create Analysis Manager
  AliAnalysisManager *runAnalysis = new AliAnalysisManager("Heavy Flavour Electron Analysis");
  runAnalysis->SetCommonFileName(output_file.Data());
  runAnalysis->SetInputEventHandler(new AliESDInputHandler); 
  if(MC) runAnalysis->SetMCtruthEventHandler(new AliMCEventHandler);
  AliAnalysisAlien *alienhandler = CreateAlienHandler(proofmode);
  printf("alienhandler %p\n", alienhandler);
  runAnalysis->SetGridHandler(alienhandler);
  //return;
  
  // Specify input (runs or dataset)
  if(!proofmode){
    // Query sample ID and runs
    TString sample;
    TArrayI listofruns;
    DecodeDataString(data, sample, listofruns);
    AddInput(alienhandler, sample, listofruns, MC);
  } else {
    alienhandler->SetProofDataSet(data);
  }

  // Add Tasks
  gROOT->LoadMacro(Form("%s/OADB/macros/AddTaskPhysicsSelection.C", gSystem->Getenv("ALICE_ROOT")));
  gROOT->LoadMacro(Form("%s/PWG3/hfe/macros/AddTaskHFE.C", gSystem->Getenv("ALICE_ROOT")));
  AddTaskPhysicsSelection(MC);
  AddTaskHFE();     // @TODO: MC and PbPb flag to be fixed

  // Run Analysis
  TString anamode = proofmode ? "proof" : "grid";
  if(runAnalysis->InitAnalysis()){
    runAnalysis->PrintStatus();
    runAnalysis->StartAnalysis(anamode);
  }
}
