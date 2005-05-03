void loadlibs () 
{
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libmicrocern");
  gSystem->Load("libpdf");
  gSystem->Load("libpythia6");
  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libEGPythia6");

  gSystem->Load("libRAWData");
  gSystem->Load("libESD");
  gSystem->Load("libSTEER");
  gSystem->Load("libEVGEN");
  gSystem->Load("libFASTSIM");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libSTRUCT");
  gSystem->Load("libMUONmapping");
  gSystem->Load("libMUONgeometry");
  gSystem->Load("libMUONbase");

  gSystem->Load("libMUONsim");
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");

  gSystem->Load("libMUONrec");
  
}
