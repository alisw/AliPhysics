// ---------------------------------------------------------------------
// Example of macro to run AliComparisonTask (PWG1 library). 
// The task can be run locally or on the Proof (aProof == kTRUE).
// 
// To run AliComparisonTask on Proof user must set the enviroment variables 
// and load AliRoot libraries using e.g. LoadMyLibs.C macro.
//
// As the input the AliComparisonTask needs a chain of root files with ESDcmpTracks
// tree inside.
// The output file (Output.root) contains MC vs Reconstruction comparison objects
// filled under conditions set in RunAliComparisonTask.C
//
//----------------------------------------------------------------------

/*
 
1. Example (run locally):

//-- Load toolkit
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
AliXRDPROOFtoolkit tool;

//-- Make chain 
TChain * chain = tool.MakeChain("cmpESDTracks.txt","ESDcmpTracks","",1000,0); 
 
//-- Run AliComparisonTask task
gROOT->LoadMacro("RunAliComparisonTask.C");
RunAliComparisonTask(chain, kFALSE);


// ----------------------------------------------------------------------

2. Example (run on Proof/GSI):

// -- Connect to Proof as USER_NAME
TProof::Open("USER_NAME@gsiaf.gsi.de");

// -- Set Proof ROOT version (e.g. 5.18/00a)
gProof->GetManager()->SetROOTVersion("5.18/00a");

// -- Load AliRoot Libraries
gProof->Exec(Form("TString str(gSystem->ExpandPathName(\"%s\")); gSystem->Setenv(\"ALICE_ROOT\", str);", gSystem->ExpandPathName("$ALICE_ROOT")), kTRUE);
gProof->AddDynamicPath(Form("%s/lib/tgt_linuxx8664gcc", gSystem->ExpandPathName("$ALICE_ROOT")));
gProof->Exec(Form("gROOT->Macro(\"%s/PWG1/Macros/LoadMyLibs.C\")",gSystem->ExpandPathName("$ALICE_ROOT")),kTRUE);

// -- Load toolkit
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
AliXRDPROOFtoolkit tool;


// -- Make chain
TChain * chain = tool.MakeChain("cmpESDTracks.txt","ESDcmpTracks","",1000,0); 

// -- Run AliComparisonTask task
gROOT->LoadMacro("RunAliComparisonTask.C");
RunAliComparisonTask(chain, kTRUE);

*/

//-----------------------------------------------------------------------------

void RunGSI(Bool_t aProof=kFALSE) {

  // -- RUN  LOCALLY
  if(aProof == kFALSE) 
  {
    // -- Load toolkit
    gSystem->Load("/usr/local/grid/XRootd/GSI/lib64/libXrdClient.so");
    gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
    gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
    AliXRDPROOFtoolkit tool;
    
    // -- Make chain
    TChain * chain = tool.MakeChain("cmpESDTracks_post_v4-16-Rev-05.txt","ESDcmpTracks","",10,0);
    
    // -- Run AliComparisonTask task
    gROOT->LoadMacro("RunAliComparisonTask.C");
    if(chain) {
      RunAliComparisonTask(chain, kFALSE);
	} else {
      AliDebug(AliLog::kError, "ERROR: Cannot get chain");
	}
  }
  //-- RUN ON PROOF
  else 
  {
    // -- Connect to Proof as USER_NAME
    TProof::Open("jacek@gsiaf.gsi.de");
     
    // -- Set Proof ROOT version (e.g. 5.18/00a)
    gProof->GetManager()->SetROOTVersion("5.18/00a");
    
    // -- Load AliRoot Libraries
    gProof->Exec(Form("TString str(gSystem->ExpandPathName(\"%s\"));gSystem->Setenv(\"ALICE_ROOT\",str);",gSystem->ExpandPathName("$ALICE_ROOT")),kTRUE);
    gProof->AddDynamicPath(Form("%s/lib/tgt_linuxx8664gcc", gSystem->ExpandPathName("$ALICE_ROOT")));
    gProof->Exec(Form("gROOT->Macro(\"%s/PWG1/Macros/LoadMyLibs.C\")",gSystem->ExpandPathName("$ALICE_ROOT")),kTRUE);
    
    // -- Load toolkit
    gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
    gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
    AliXRDPROOFtoolkit tool;
    
    // -- Make chain
    TChain * chain = tool.MakeChain("cmpESDTracks.txt","ESDcmpTracks","",1000,0);
    
    // -- Run AliComparisonTask task
    gROOT->LoadMacro("RunAliComparisonTask.C");
    if(chain) {
      RunAliComparisonTask(chain, kTRUE);
	} else {
      AliDebug(AliLog::kError, "ERROR: Cannot get chain");
	}
  }
}
