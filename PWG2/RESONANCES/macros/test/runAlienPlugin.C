//
// This is an example steering macro for running RSN analysis task
// with the AliEn plugin to launch a multiple analysis.
//
// Inputs:
//   - runList     = list of runs to be processed
//   - runPath     = path containing the runs
//   - runMode     = AliEn plugin run mode
//   - pluginMacro = macro which loads and initializes the plugin
//   - addTaskName = name of the macro to add the RSN analysis task
//                   (assumed to have inside it a function named like the file)
//   - inputSource = name of the file containing all the inputs
//                   ---> to run on a local collection, the collection file 
//                        must contain on each line the full path 
//                        of one input file and it must have the ".txt" extension
//                   ---> to run on an AliEn collection, the collection file must be an XML
//                        file collection like those built from the "find -x" method in aliensh.
//   - dataLabel   = a label which is used to know what kind of data are being read
//                   (it is propagated to the 'addTask' macro for eventual setting up of something
//   - outName     = name for the file with RSN package outputs (without ROOT extension)
//                   in this case it is fundamental to define the names of all plugin objects
//
// Notes:
//   - in case the source is an ESD, and if inputs are a MC production
//     the MC input handler is created by default
// 
//
// In principle, the user should never modify this macro. 
//
void runAlienPlugin
(
  const char *runList     = "117112-117116-117099-117220-117048-117109-117060-117054-117065",
  const char *runPath     = "/alice/data/2010/LHC10b",
  const char *runMode     = "terminate",
  const char *pluginMacro = "PluginDataByRun.C",
  const char *addTaskName = "AddAnalysisTaskRsnTest.C",
  const char *dataLabel   = "7TeV_pass2_data_ESD",
  const char *outName     = "rsnTest"
)
{
  
  // convert the last argument into a BOOL variable
  TString strDataLabel(dataLabel);
  Bool_t isESD = strDataLabel.Contains("ESD");
  Bool_t isAOD = strDataLabel.Contains("AOD");
  Bool_t isSim = strDataLabel.Contains("sim");   
  
  //AliLog::SetGlobalDebugLevel(AliLog::kDebug+2);

  // load compiled libraries (for aliroot session)
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG2resonances.so");
  
  //
  // === PLUGIN CONFIGURATION =====================================================================
  //
  
  // check token
  if (!AliAnalysisGrid::CreateToken()) return;
  
  // load and execute plugin configuration macro
  // pass to the macro, as FIRST argument, the common name
  // which is used for the output, since it must be equal
  // to the one defined here for the common output (for merging)
  TString splugin(pluginMacro);
  gROOT->LoadMacro(pluginMacro);
  splugin.ReplaceAll(".C", Form("(\"%s\",\"%s\",\"%s\")", outName, runList, runPath));
  AliAnalysisAlien *plugin = (AliAnalysisAlien*)gROOT->ProcessLine(splugin);
  
  // set run mode
  plugin->SetRunMode(runMode);
  
  //
  // === ANALYSIS MANAGER CONFIGURATION ===========================================================
  //

  // create analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("taskRsnTest");
  mgr->SetGridHandler(plugin);
  mgr->SetCommonFileName(Form("%s.root", outName));
  
  // create input handler
  if (isESD)
  {
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdHandler);
    // if possible, create also MC handler
    if (isSim)
    {
      AliMCEventHandler *mcHandler  = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcHandler);
    }
  }
  else if (isAOD)
  {
    AliAODInputHandler *aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);
  }
  else
  {
    ::Error("Required an ESD or AOD input data set");
    return;
  }
  
  //
  // === ANALYSIS TASK CREATION AND INCLUSION =====================================================
  //
  
  // add event selection for data
  gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isSim);
  
  // add task macro
  gROOT->ProcessLine(Form(".x %s(\"%s\")", addTaskName, dataLabel));

  // initialize and start analysis
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");
}

//_________________________________________________________________________________________________
Bool_t LoadPars(const char *parList, const char *path)
{
//
// Load PAR libraries locally
// ---
// Arguments:
//  - parList = list of PARs without extension, separated by ':'
//  - path    = path where PARs are stored
//

  // store position of working directory
  TString ocwd = gSystem->WorkingDirectory();

  // tokenize list
  TString     pars(parList);
  TObjArray  *array = pars.Tokenize(":");

  // loop on list
  TObjString *ostr;
  TString     str;
  Char_t      parName[200], parFile[200];
  for (Int_t i = 0; i < array->GetEntriesFast(); i++)
  {
    ostr = (TObjString*) array->At(i);
    str = ostr->GetString();
    sprintf(parName, "%s", str.Data());
    sprintf(parFile, "%s/%s.par", path, str.Data());

    // check that file exists
    if (!gSystem->AccessPathName(parFile))
    {
      // explode tar-ball and enter it
      gROOT->ProcessLine(Form(".! tar xzf %s", parFile));
      gSystem->ChangeDirectory(Form("%s", parName));
      // checks for BUILD.sh and execute it
      if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh"))
      {
        ::Info("", ">> Building PARs: %s", parName);
        if (gSystem->Exec("PROOF-INF/BUILD.sh"))
        {
          ::Error("LoadPars", Form("BUILD.sh error for '%s'", parFile));
          gSystem->ChangeDirectory(ocwd);
          return kFALSE;
        }
      }
      // check and execute SETUP.C
      if (!gSystem->AccessPathName("PROOF-INF/SETUP.C"))
      {
        ::Info("", ">> Setting up PARs: %s", parName);
        if (gROOT->Macro("PROOF-INF/SETUP.C"))
        {
          Error("LoadPars", Form("SETUP.C error for '%s'", parFile));
          gSystem->ChangeDirectory(ocwd);
          return kFALSE;
        }
      }
    }
    else
    {
      Error("LoadParsLocal", Form("File '%s' not found", parFile));
      return kFALSE;
    }
  }

  gSystem->ChangeDirectory(ocwd);
  return kTRUE;
}
