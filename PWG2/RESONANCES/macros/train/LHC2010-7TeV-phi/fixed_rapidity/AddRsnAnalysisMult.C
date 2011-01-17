//
// This macro serves to add the RSN analysis task to the steering macro.
//
// Inputs:
//   - dataLabel   = a string with informations about the type of data
//                   which could be needed to be ported to the config macro
//                   to set up some cuts
//   - configMacro = macro which configures the analysis; it has *ALWAYS*
//                   defined inside a function named 'RsnConfigTask()',
//                   whatever the name of the macro itself, whose first two
//                   arguments must have to be the task and the 'dataLabel' argument.
//
Bool_t AddRsnAnalysisMult
(
  const char *options,
  const char *configs = "RsnConfigNoSA.C RsnConfigSA.C",
  const char *path    = "$(ALICE_ROOT)/PWG2/RESONANCES/macros/train/LHC2010-7TeV-phi"
)
{
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  // initialize multiplicity bins
  Int_t multMin[] = {0, 0, 6, 10, 15, 23       ,  6, 11, 21, 31, 51,  76, 101, 126, 151, 201      };
  Int_t multMax[] = {0, 5, 9, 14, 22, 100000000, 10, 20, 30, 50, 75, 100, 125, 150, 200, 100000000};
  Int_t nmult = sizeof(multMin) / sizeof(multMin[0]);
  
  // initialize several tasks, each one with different multiplicity cut
  // and all with the same primary vertex + pile-up cut
  for (Int_t i = 0; i < nmult; i++)
  {
    // create the task and connect with physics selection
    AliRsnAnalysisSE *task = new AliRsnAnalysisSE(Form("RsnAnalysis_%d", i));
    task->SetZeroEventPercentWarning(100.0);
    task->SelectCollisionCandidates();

    // add the task to manager
    mgr->AddTask(task);

    // load and execute all required configuration macroes in the string (arg #2)
    TString    sList   = configs;
    TObjArray *list    = sList.Tokenize(" ");
    Int_t      nConfig = list->GetEntries();
    Int_t      iConfig = 0;
    for (iConfig = 0; iConfig < nConfig; iConfig++)
    {
      TObjString *ostr = (TObjString*)list->At(iConfig);
      
      // the config macro is assumed to be stored in the path in argument #3
      // and to have three arguments: task name, a free string of options and the path where it is stored
      // --> all of them is a string, and then it must be passed with the quote marks
      const char *macro     = ostr->GetString().Data();
      const char *argName   = Form("\"%s\"", task->GetName());
      const char *argOption = Form("\"%s\"", options);
      const char *argPath   = Form("\"%s\"", path);
      const char *argMult   = Form("%d, %d", multMin[i], multMax[i]);
      gROOT->ProcessLine(Form(".x %s/%s(%s,%s,%s, %s)", path, macro, argName, argOption, argPath, argMult));
    }

    // connect input container according to source choice
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    // create paths for the output in the common file
    Char_t commonPath[500];
    sprintf(commonPath, "%s", AliAnalysisManager::GetCommonFileName());

    // create containers for output
    AliAnalysisDataContainer *outputInfo = mgr->CreateContainer(Form("RsnInfo_%d", i), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
    AliAnalysisDataContainer *outputHist = mgr->CreateContainer(Form("RsnHist_%d", i), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
    mgr->ConnectOutput(task, 1, outputInfo);
    mgr->ConnectOutput(task, 2, outputHist);
    
  } // end loop on tasks

  return kTRUE;
}
