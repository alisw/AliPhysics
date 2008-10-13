//
// This macro creates the AnalysisManager containing specific tasks for resonances.
// The created object is used in another macro which runs effectively the analysis.
//
// Arguments:
//  - a flag to know if the macro is being executed in PROOF or not (should not be touched)
//  - output filename
//  - a string containing all macroes which define pair manager for the analysis,
//    which are supposed to be defined in the standard ROOT style, that is a macro
//    with a function with the same name; they can be more than one, separated by colons (':')
//
// In this macro, the user should play with the following things:
//  - what pairs to use (second argument)
//  - what kind of source data to use (third argument)
//

AliAnalysisManager* CreateAnalysisManager
(
  Bool_t                             isProof,
  const char                        *fileOut   = "rsn.root",
  const char                        *macroList = "CreatePairsPhi.C",
  AliRsnAnalysisTaskBase::EInputType inputType = AliRsnAnalysisTaskBase::kESDMC,
)
{
  // default path containing all macroes described in 2nd argument:
  // ...for local analysis, it goes to the AliRoot path:
  TString strMacroListPath(Form("%s/PWG2/RESONANCES/macros", getenv("ALICE_ROOT")));
  // ...for PROOF analysis from lxplus, the macro is probably in the same directory
  if (isProof) strMacroListPath = ".";
  // (of course one can play with this variable)
  cout << "Loaded macro path: " << strMacroListPath.Data() << endl;
  
  // initialize analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("RsnTEST");
  if (!mgr) AliError("no Analysis Mgr");

  // initialize analysis task
  // the method "SetInputType" is used to tell the task what kind of input is read;
  // in this case we are reading ESD data with MC info, otherwise one must use
  // kESD, kAOD, kMC or kRSN to read files already in the internal non-standard AOD of this package
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("AliRsnAnalysisSE");
  task->SetInputType(inputType, mgr, kTRUE);
  task->SetNumOfEventsInBuffer(2000);
  
  // Settings for the ESD-to-RSN conversion
  // in this case we reject splitted tracks (worst one)
  AliRsnReader *reader = task->GetReader();
  reader->SetCheckSplit(kTRUE);
  
  // PID settings:
  // define prior probabilities
  AliRsnPID *pid = task->GetPID();
  pid->SetPriorProbability(AliRsnPID::kElectron, 0.02);
  pid->SetPriorProbability(AliRsnPID::kMuon,     0.02);
  pid->SetPriorProbability(AliRsnPID::kPion,     0.83);
  pid->SetPriorProbability(AliRsnPID::kKaon,     0.07);
  pid->SetPriorProbability(AliRsnPID::kProton,   0.06);
  pid->SetMaxPt(10.0);
  pid->SetMinProb(0.5);
  
  // tokenize the argument to retrieve all macroes which define pair managers
  TString strMacros(macroList);
  TObjArray *macros = strMacros.Tokenize(":");
  TObjString *macro = 0x0;
  TObjArrayIter next(macros);
  while ( (macro = (TObjString*)next()) ) {
    TString str = macro->GetString();
    // different methods for PROOF and others
    if (isProof) {
      cout << "PROOF: loading " << str.Data() << endl;
      gProof->Load(Form("%s/%s", strMacroListPath.Data(), str.Data())); 
    }
    else {
      cout << "LOCAL: loading " << str.Data() << endl;
      gROOT->LoadMacro(Form("%s/%s", strMacroListPath.Data(), str.Data()));
    }
    str.ReplaceAll(".C", "()");
    AliRsnPairMgr *pairMgr = (AliRsnPairMgr*)gROOT->ProcessLine(str.Data());
    task->AddPairMgr(pairMgr);
  }

  // define containers:
  // the AOD container of the AliAnalysisTaskSE is created as required
  // but is not used for output, while another output container for histograms is created
  AliAnalysisDataContainer *input  = mgr->CreateContainer("in", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *dummy  = mgr->CreateContainer("dummy", TTree::Class(), AliAnalysisManager::kOutputContainer, "default");
  AliAnalysisDataContainer *output = mgr->CreateContainer("Histograms", TList::Class(), AliAnalysisManager::kOutputContainer, fileOut);

  // connect containers to AnalysisManager and return it
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, input);
  mgr->ConnectOutput(task, 0, dummy);
  mgr->ConnectOutput(task, 1, output);

  return mgr;
}
