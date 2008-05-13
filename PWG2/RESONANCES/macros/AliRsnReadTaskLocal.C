//=========================================================================
// This macro loops on the entries of a TChain (argument) containing ESDs
// and saves a file containing a TTree of AliRsnEvents.
//=========================================================================

void AliRsnReadTaskLocal()
{
    // by default, assume that a file named "AliESDs.root"
    // exists in the working directory and connects to it
    TChain *analysisChain = new TChain("esdTree");
    analysisChain->Add("AliESDs.root");
    
    // load read macro
    TString str(getenv("ALICE_ROOT"));
    str.Append("/PWG2/RESONANCES/macros/AliRsnReadTask.C");
    gROOT->LoadMacro(str.Data());
    
    AliRsnReadTask(analysisChain);
}
