void LoadLibraries()
{
  //--------------------------------------
  //  // Load the needed flow libraries most of them already loaded by aliroot
  //    //--------------------------------------
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libXMLParser");
  gSystem->Load("libProof");
  gSystem->Load("libMinuit");


  // for AliRoot
  //gSystem->Load("libSTAT");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTEER");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGmuon");
  gSystem->Load("libTPCbase");
  gSystem->Load("libTOFbase");
  gSystem->Load("libTOFsim");
  gSystem->Load("libTOFrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libT0base");
  gSystem->Load("libT0rec");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");

  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  
} // end of void LoadLibrariesRF(const libModes mode)
