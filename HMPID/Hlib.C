void Hlib()
{
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  
  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  gSystem->Load("libEGPythia6");

  gSystem->Load("libNet");
  gSystem->Load("libTree");
 
  gSystem->Load("libProof");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libCDB");
  
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  
  gSystem->Load("libSTEER");
  
  gSystem->Load("libRAWDatasim");

  gSystem->Load("libHMPIDbase");
  gSystem->Load("libHMPIDsim");
  gSystem->Load("libHMPIDrec");

//  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/HMPID"); 
  gInterpreter->AddIncludePath("$HOME/HMPID");
  gInterpreter->AddIncludePath("$ALICE_ROOT/include");
  gInterpreter->AddIncludePath("$ALICE_ROOT/HMPID");

}





