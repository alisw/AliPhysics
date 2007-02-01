void Hlib()
{
  gSystem->Load("libMinuit.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libSTEER.so");
  gSystem->Load("libCDB.so");

  gSystem->Load("libRAWDatasim.so");
  
  gSystem->Load("libHMPIDbase.so");
  gSystem->Load("libHMPIDsim.so");
  gSystem->Load("libHMPIDrec.so");
  
  gROOT->LoadMacro("Hmenu.C");
}
