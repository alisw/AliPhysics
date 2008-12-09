//____________________________________________________________________
//
// $Id$
//
// Script to compile (using ACLic) and load a script.  It sets the
// include path to contain the relevant directories. 
//
/** @defgroup FMD_script Scripts
    @brief Scripts for the FMD 
*/
/** @defgroup simple_script Simple scripts
    @ingroup FMD_script
    @brief Scripts for the FMD 
*/
/** Compile an FMD script using ACLic
    @param script Script to compile
    @param option Compile option 
    @ingroup FMD_script
*/
Bool_t
Compile(const char* script, Option_t* option="g")
{
  if (!script || script[0] == '\0') { 
    std::cerr << "No script to compile!" << std::endl;
    return kFALSE;
  }
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libFMDanalysis.so");
  gSystem->Load("libFMDutil.so");
  gSystem->Load("libFMDflow.so");
  TString macroPath(gROOT->GetMacroPath());
  macroPath.Append(":${ALICE_ROOT}/FMD/scripts");
  gROOT->SetMacroPath(macroPath.Data());
  gSystem->SetIncludePath("-I`root-config --incdir` "
			  "-I${ALICE_ROOT} " 
			  "-I${ALICE_ROOT}/include " 
			  "-I${ALICE_ROOT}/FMD "
			  "-I${ALICE_ROOT}/geant3/TGeant3");
  Long_t ret = gROOT->ProcessLine(Form(".L %s+%s", script, option));
  return ret == 0;
}

//____________________________________________________________________
//
// EOF
// 
