//____________________________________________________________________
//
// $Id$
//
// Script to compile (using ACLic) and load a script.  It sets the
// include path to contain the relevant directories. 
//
/** @defgroup FMD_script Scripts
    @brief Scripts for the FMD 
    @ingroup FMD
*/
/** @defgroup FMD_simple_script Simple scripts
    @ingroup FMD_script
    @brief Scripts for the FMD 
*/
/** Compile an FMD script using ACLic
    @param script Script to compile
    @param option Compile option 
    @ingroup FMD_script
*/
void
AddInclude(const char* what)
{
  TString path(gSystem->GetIncludePath());
  if (path.Contains(what)) return;
  gSystem->AddIncludePath(what);
}
Bool_t
Compile(const char* script=0, Option_t* option="g")
{
  if (!script || script[0] == '\0') { 
    std::cerr << "No script to compile!" << std::endl;
    return kFALSE;
  }
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libFMDutil");
  TString macroPath(gROOT->GetMacroPath());
  macroPath.Append(":${ALICE_ROOT}/FMD/scripts");
  gROOT->SetMacroPath(macroPath.Data());
  AddInclude("-I`root-config --incdir`");
  AddInclude("-I${ALICE_ROOT}");
  AddInclude("-I${ALICE_ROOT}/include");
  AddInclude("-I${ALICE_ROOT}/FMD");
  AddInclude("-I${ALICE_ROOT}/geant3/TGeant3");
  AddInclude("-I${ALICE_ROOT}/../master-src");
  AddInclude("-I${ALICE_ROOT}/../master-src/FMD");
  AddInclude("-I${ALICE_ROOT}/../master-src/RAW");
  Long_t ret = gROOT->ProcessLine(Form(".L %s+%s", script, option));
  return ret == 0;
}

//____________________________________________________________________
//
// EOF
// 
