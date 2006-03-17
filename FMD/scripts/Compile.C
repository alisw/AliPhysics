//____________________________________________________________________
//
// $Id$
//
// Script to compile (using ACLic) and load a script.  It sets the
// include path to contain the relevant directories. 
//
void
Compile(const char* script, Option_t* option="g")
{
  gSystem->Load("libFMDutil.so");
  gSystem->SetIncludePath("-I`root-config --incdir` "
			  "-I${ALICE_ROOT}/include " 
			  "-I${ALICE_ROOT}/FMD "
			  "-I${ALICE_ROOT}/geant3/TGeant3");
  gROOT->ProcessLine(Form(".L %s+%s", script, option));
}

//____________________________________________________________________
//
// EOF
// 
