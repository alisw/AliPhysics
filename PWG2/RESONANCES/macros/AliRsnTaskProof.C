//=========================================================================
// This macro loops on a dataset in CAF and performs a single-step
// resonance analysis, which reads the ESD events (with MC if available)
// and saves directly the invariant mass spectra.
//
// All required configurations (reader, PID, pairs, dataSet) is defined
// inside a macro which must be built in the same way as 
// the example named "PhiExample.C" in the "macros/config" directory.
//
// The output file name is defined from the config file name
// irrespectively of its location.
//=========================================================================

//
// Macro to build an analysis chain from a dataset
//
TChain* CreateChainFromDataSet
(TFileCollection* coll, const char* treeName, const Int_t nMaxFiles)
{
  TIter iter(coll->GetList());
  TChain* target = new TChain(treeName);
	
  TFileInfo* fileInfo = 0;
  Int_t nFilesAdded = 0;
  while ((fileInfo = dynamic_cast<TFileInfo*> (iter())) && nFilesAdded < nMaxFiles)
  {
    if (fileInfo->GetFirstUrl()) {
	  target->Add(fileInfo->GetFirstUrl()->GetUrl());
	  nFilesAdded++;
	}
  }
	
  Printf("Added %d files to chain", target->GetListOfFiles()->GetEntries());
  return target;
}

//
// Core macro which executes analysis
//
void AliRsnTaskProof
(
  const char   *fileOut         = "rsn.root",
  const char   *macro           = "CreateAnalysisManager.C",
  Int_t         nFilesToProcess = 2,
  const Char_t *dataSetName     = "/COMMON/COMMON/LHC08c11_10TeV_0.5T",
  const Char_t *treeName        = "esdTree"
)
{
  // path for the "CreateAnalysisManager.C" macro which creates the AnalysisManager
  // in this case, assuming that one copies this in an area in lxplus, the path
  // is set to the same directory
  //TString strTaskPath(Form("%s/PWG2/RESONANCES/macros", getenv("ALICE_ROOT")));
  TString strTaskPath(".");
  
  // connect to CAF
  TProof::Open("alicecaf");

  // upload packages
  gProof->UploadPackage("STEERBase.par");
  gProof->EnablePackage("STEERBase");
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");
  gProof->UploadPackage("ANALYSISalice.par");
  gProof->EnablePackage("ANALYSISalice");
  gProof->UploadPackage("PWG2resonances.par");
  gProof->EnablePackage("PWG2resonances");
	
  // create chains for processing
  TChain *analysisChain = 0;
  TFileCollection *fc = gProof->GetDataSet(dataSetName)->GetStagedSubset();
  if (nFilesToProcess <= 0) nFilesToProcess = 10000;
  analysisChain = CreateChainFromDataSet(fc, treeName, nFilesToProcess);
  Printf("Found %d entries", analysisChain->GetEntries());
	
  // load and execute the macro which creates the AnalysisManager
  // this macro is expected to define a function with the standard name
  // *** "CreateAnalysisManager()" *** [without arguments],
  // irrespectively of the filename of the macro itself,
  // which returns an AliAnalysisManager object
  cout << Form("%s/%s", strTaskPath.Data(), macro) << endl;
  gROOT->LoadMacro(Form("%s/%s", strTaskPath.Data(), macro));
  AliAnalysisManager *mgr = CreateAnalysisManager(kTRUE);
  if (!mgr) Error("", "no Analysis Mgr");
	
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("proof", analysisChain);
}
