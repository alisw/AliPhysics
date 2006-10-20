// $Id$
//
// Macro for loading libraries for mapping in AliRoot.
// It sets also the environment variable used to find mapping data.

void mlibs() 
{
  // add include path
  // gSystem->SetIncludePath(" -I$MINSTALL/include");
  gSystem->SetIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/MUON/mapping");

  // load Root libraries
  gSystem->Load("libPhysics"); 
  
  // load mapping library
  gSystem->Load("libMUONmapping"); 

  // set path to mapping data
  if (! gSystem->Getenv("MINSTALL")) {    
    TString dirPath = gSystem->Getenv("ALICE_ROOT");
    dirPath += "/MUON/mapping"; 
    AliMpFiles::SetTopPath(dirPath);
    gSystem->Setenv("MINSTALL", dirPath.Data());
    // cout << "AliMpFiles top path set to " << dirPath << endl;	  
  }  
}  
