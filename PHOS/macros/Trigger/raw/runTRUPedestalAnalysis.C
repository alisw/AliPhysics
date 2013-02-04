#include <TSystem.h>
#include <TROOT.h>

void runTRUPedestalAnalysis(char* filesList = "files.txt")
{
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PHOS");
  gSystem->AddIncludePath("-I$ALICE_ROOT/RAW");

  gSystem->Load("libPHOSrec");
  
  gROOT->LoadMacro("AliTRUPedestalOutput.cxx+g");
  gROOT->LoadMacro("AliTRUPedestalAnalysis.cxx+g");

  gROOT->LoadMacro("truPedestalAnalysis.C+g");
  
  // for if the files are on alien, use responsibly!
  // TGrid::Connect("alien://");
  // AliCDBManager::Instance()->SetDefaultStorage("raw://");

  truPedestalAnalysis(filesList);
}
