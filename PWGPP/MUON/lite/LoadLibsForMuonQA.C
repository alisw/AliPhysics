#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"

#endif

//_____________________________
void LoadRootAnalysis()
{
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libProof");
}

//_____________________________
void LoadAnalysis(const char* option = "")
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libOADB.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  TString opt(option);
  opt.ToUpper();
  if ( opt.Contains("PWG") ) {
    gSystem->Load("libPWGmuon.so");
  }
  if ( opt.Contains("PWGPP") ) {
    gSystem->Load("libPWGPPMUONlite.so");
  }
}

//_____________________________
void IncludeAliroot()
{
  TString envList[3] = {"ALICE_ROOT","ALICE_INSTALL","ALICE_BUILD"};
  for ( Int_t ienv=0; ienv<3; ienv++ ) {
    if ( ! gSystem->Getenv(envList[ienv].Data()) ) continue;
    if ( gSystem->AccessPathName(gSystem->ExpandPathName(Form("${%s}/include",envList[ienv].Data()))) ) continue;
    gSystem->AddIncludePath(Form("-I${%s}/include",envList[ienv].Data()));
  }
}

//_____________________________
void IncludeMuon()
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/MUON");
}

//_____________________________
void LoadLibsForMuonQA ( const char* option )
{
  TString opt(option);
  opt.ToLower();
  if ( opt.Contains("maketrend") ) {
    IncludeAliroot();
    gSystem->AddIncludePath("-I${ALICE_ROOT}/PWGPP/MUON/lite");
    LoadRootAnalysis();
    LoadAnalysis("PWGPP");
  }
  if ( opt.Contains("trigtrend") ) {
    IncludeAliroot();
    IncludeMuon();
    LoadAnalysis("PWG");
    gSystem->Load("libPWGmuondep.so");
  }
  if (opt.Contains("tracktrend") ) {
    IncludeAliroot();
    LoadAnalysis("PWG");
  }
}
