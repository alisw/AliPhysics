void loadSteeringLibs () 
{
  gSystem->Load("libmicrocern");
  gSystem->Load("libPhysics");
  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  gSystem->Load("libRAW");
  gSystem->Load("libESD");
  gSystem->Load("libSTEER");
  gSystem->Load("libCONTAINERS");
}
