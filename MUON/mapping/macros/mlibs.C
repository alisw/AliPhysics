// $Id$
//
// Macro for loading libraries for mapping in AliRoot

void mlibs() 
{
  // add include path
  // gSystem->SetIncludePath(" -I$MINSTALL/include");
  
  // set path to mapping
  if (! gSystem->Getenv("MINSTALL")) {    
    TString dirPath = gSystem->Getenv("ALICE_ROOT");
    dirPath += "/MUON/mapping"; 
    gSystem->Setenv("MINSTALL", dirPath.Data());
    //cout << "AliMpFiles top path set to " << dirPath << endl;	  
  }

  // load Root libraries
  gSystem->Load("libPhysics"); 
  
  // load mapping library
  gSystem->Load("libMUONmapping"); 
}  
