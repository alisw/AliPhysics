//____________________________________________________________________
//
// $Id: Compile.C 30305 2008-12-09 05:45:53Z cholm $
//
// Script to compile (using ACLic) and load a script.  It sets the
// include path to contain the relevant directories. 
//
/** 
 * Compile an FMD script using ACLic
 *
 * @param script Script to compile
 * @param option Compile option 
 * 
 * @return true on success
 * 
 * @ingroup pwglf_forward_scripts
 */
Bool_t
Compile(const char* script, Option_t* option="g")
{
  if (!script || script[0] == '\0') { 
    std::cerr << "No script to compile!" << std::endl;
    return kFALSE;
  }
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGLFforward2");
  TString macroPath(gROOT->GetMacroPath());
  macroPath.Append(":${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2");
  macroPath.Append(":${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2/scripts");
  gROOT->SetMacroPath(macroPath.Data());
  gSystem->SetIncludePath("-I`root-config --incdir` "
			  "-I${ALICE_PHYSICS} " 
			  "-I${ALICE_PHYSICS}/include " 
			  "-I${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2 "
			  "-I${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2 ");
  Long_t ret = gROOT->ProcessLine(Form(".L %s+%s", script, option));
  return ret == 0;
}

//____________________________________________________________________
//
// EOF
// 
