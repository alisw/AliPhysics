void loadlibraw () 
{
  gSystem->Load("libPhysics");
  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  gSystem->Load("libRAWData");
  gSystem->Load("libESD");
  gSystem->Load("libSTEER");

  gSystem->Load("libMUONmapping");
  gSystem->Load("libMUONgeometry");
  gSystem->Load("libMUONraw");

}
