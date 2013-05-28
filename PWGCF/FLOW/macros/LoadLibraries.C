void LoadLibraries()
{
  //--------------------------------------
  //  // Load the needed flow libraries most of them already loaded by aliroot
  //    //--------------------------------------
  //gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");

  // for AliRoot
  //gSystem->Load("libSTAT");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
} // end of void LoadLibrariesRF(const libModes mode)
