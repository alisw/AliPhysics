void Hlib()
{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libSTEER.so");
  gSystem->Load("libCDB.so");

  gSystem->Load("libEVGEN.so");
  gSystem->Load("libSTRUCT.so");
  
  gSystem->Load("libANALYSIS.so");

  gSystem->Load("libRAWDatasim.so");
  
  gSystem->Load("libHMPIDbase.so");
  gSystem->Load("libHMPIDsim.so");
  gSystem->Load("libHMPIDrec.so");

  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/HMPID"); 
}
