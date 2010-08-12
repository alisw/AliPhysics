//=============================================================================
//
// This is an example of steering macro for running RSN analysis task
// in a PROOF cluster, using a dataset stored there.
//
// Usage requirements:
// -- 1: run with root
// -- 2: have available the PARs of: STEERBase, ESD, AOD, ANALYSIS
//                                   ANALYSISalice, CORRFW and PWG2resonances
//
//=============================================================================

// enum to define data types
enum ERsnData
{
  kRsnESD,       // ESD
  kRsnAOD        // AOD
};

//_________________________________________________________________________________________________
//
// Steering macro.
// The user must call this one to run the job.
// ---
// Arguments:
//  - nRead       = how many events to read
//  - nSkip       = how many events to skip (starting from first in collection)
//  - inputSource = dataset name
//  - outName     = output file name
//  - dataType    = kESD or kAOD (see enum above)
//  - isMC        = tells if this is a MonteCarlo (otherwise is data)
//  - pathStd     = path of {STEERBase|ESD|AOD|ANALYSIS|ANALYSISalice|CORRFW}.par
//  - pathRsn     = path of PWG2resonances.par
//
// When 'isMC' is true, the MC handler is created by default.
//
void runProof
(
  Int_t       nRead       = 1e5,
  Int_t       nSkip       = 0,
  const char *addMacro    = "AddAnalysisTaskRsnTest.C",
  //const char *inputSource = "/ALICE/pp000900/MC_LHC09d10_104821",
//   const char *inputSource = "/alice/sim/LHC10a12_104157",
    const char *inputSource = "/alice/data/LHC10b_000117112_p2",
  const char *outName     = "rsn_proof.root",
  ERsnData    dataType    = kRsnESD,
  Bool_t      isMC        = kFALSE,
  const char *pathStd     = "/home/pulvir/ALICE/ALIROOT/head",
  const char *pathRsn     = "/home/pulvir/ALICE/ALIROOT/head",
  const char *datalabel   = "7TeV_pass2_data_ESD"
)
{

  // connect to PROOF
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  TProof::Open("skaf.saske.sk");
  //TProof::Open("pulvir@localhost");

  // this i will do that AAF will load (needed for TOF)
  gProof->Exec("gSystem->Load(\"libXMLParser.so\");");
  
  // needed for tof too
  gProof->Exec("TGrid::Connect(\"alien:\/\/\");");
  
  // setup PARs
//   gProof->ClearPackages();

  // setup aliroot mode in AAF (for now using SIM mode since ALIROOT mode doesn't load correctly (this is tmp))
  TList *listAliroot = new TList();
  listAliroot->Add(new TNamed("ALIROOT_MODE", "SIM"));
  
  Bool_t usePWG2resonancesPAR = kTRUE;
  TString alirootVer = "VO_ALICE@AliRoot::v4-20-03-AN-proof";
  
  if (usePWG2resonancesPAR) {
    listAliroot->Add(new TNamed("ALIROOT_EXTRA_LIBS", "ANALYSIS:ANALYSISalice:CORRFW"));
    listAliroot->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", "TOF"));
    gProof->EnablePackage(alirootVer.Data(),listAliroot);
    gProof->UploadPackage("PWG2resonances.par");
    if (gProof->EnablePackage("PWG2resonances")) {
      Error("runAAF.C","Error in PWG2resonances !!!");
      return;
    }
  }
  else {
    listAliroot->Add(new TNamed("ALIROOT_EXTRA_LIBS", "ANALYSIS:ANALYSISalice:CORRFW:PWG2resonances"));
    listAliroot->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", "PWG2/RESONANCES:TOF"));
    gProof->EnablePackage(alirootVer.Data(),listAliroot);
  }
  
  // create analysis manager and set filename
  AliAnalysisManager *mgr = new AliAnalysisManager("RsnAnalysis");
  mgr->SetCommonFileName(outName);

  // create input handler according to data type
  // if it is ESD and MC is required, add also that
  AliESDInputHandler *hesd = 0x0;
  AliAODInputHandler *haod = 0x0;
  AliMCEventHandler  *hmc  = 0x0;
  switch (dataType)
  {
    case kRsnESD:
      hesd = new AliESDInputHandler();
      mgr->SetInputEventHandler(hesd);
      if (isMC)
      {
        hmc  = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(hmc);
      }
      break;
    case kRsnAOD:
      haod = new AliAODInputHandler();
      mgr->SetInputEventHandler(haod);
      break;
    default:
      ::Error("rsnLocal.C", "Data type not supported");
      return;
  }
  
//   // add event selection for data
//   gROOT->LoadMacro("AddTaskPhysicsSelection.C");
//   AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);

  // scan list of macros
  TString        macro, taskList(addMacro);
  TObjArray     *list = taskList.Tokenize(" ");
  TObjString    *ostr = 0;
  TObjArrayIter *next = (TObjArrayIter*)list->MakeIterator();
  while ( (ostr = (TObjString*)(*next)()) )
  {
    // get tokenized string
    macro = ostr->GetString();
    Info("rsnLocal.C", "Adding macro: %s", macro.Data());
    // load the macro and execute it
    gROOT->LoadMacro(macro.Data());
    macro.ReplaceAll(".C",Form("(\"%s\");",datalabel));
    gROOT->ProcessLine(macro.Data());
  }

  // initialize analysis and run it
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis("proof", inputSource, nRead, nSkip);
}

//_________________________________________________________________________________________________
Bool_t LoadPars(const char *parList, const char *path = ".")
{
//
// Load PAR libraries in a PROOF environment.
// ---
// Arguments:
//  - parList = list of PARs without extension, separated by ':'
//  - path    = path where PARs are stored
//

  // check PROOF initialization
  if (!gProof) {
    Error("CleanPars", "gProof object not initialized");
    return kFALSE;
  }

  TString     str, pars(parList), fileName;
  TObjArray  *array = pars.Tokenize(":");
  TObjString *ostr;

  for (Int_t i = 0; i < array->GetEntriesFast(); i++)
  {
    ostr = (TObjString*) array->At(i);
    str = ostr->GetString();

    Info("", ">> Creating PAR: %s", str.Data());

    // compose filename with path and extension
    fileName = path;
    fileName += '/';
    fileName.Append(str.Data());
    fileName.Append(".par");

    // upload package (uses filename)
    gProof->UploadPackage(fileName.Data());

    // enable package (uses only PAR name)
    gProof->EnablePackage(str.Data());
  }

  return kTRUE;
}


