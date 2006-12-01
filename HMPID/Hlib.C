void h()
{
  gSystem->Load("libMinuit.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libSTEER.so");
  gSystem->Load("libCDB.so");

  gSystem->Load("libHMPIDbase.so");
  gSystem->Load("libHMPIDsim.so");
  gSystem->Load("libHMPIDrec.so");
}
