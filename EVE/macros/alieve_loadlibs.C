// $Header$

void alieve_loadlibs () 
{
  // Macro which loads the libraries needed for simulation and reconstruction

  printf("Loading ALICE libraries ...");
  fflush(stdout);

  gSystem->Load("libPhysics");

  // Uncomment the following line for macosx
  // Waiting for a better solution
  //  gSystem->Load("libg2c_sh");
  // And the following works for gcc-4.
  gSystem->Load("/usr/lib/libg2c");
  // gSystem->Load("/usr/lib/gcc/i386-redhat-linux/4.1.0/libgfortran");
  gSystem->Load("libmicrocern");
  gSystem->Load("libpdf");
  gSystem->Load("libpythia6");

  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  gSystem->Load("libEGPythia6");

  gSystem->Load("libESD");
  gSystem->Load("libSTEER");
  
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");

  // Moved before libRAW
  // gSystem->Load("libESD");
  // gSystem->Load("libSTEER");

  gSystem->Load("libEVGEN");
  gSystem->Load("libFASTSIM");
  gSystem->Load("libAliPythia6");

  gSystem->Load("libhijing");
  gSystem->Load("libTHijing");// AliGenHijingEventHeader needed by libZDCsim.so

  gSystem->Load("libSTRUCT");
  gSystem->Load("libPHOSbase");
  gSystem->Load("libPHOSsim");
  gSystem->Load("libPHOSrec");
  gSystem->Load("libMUONmapping");
  gSystem->Load("libMUONgeometry");
  gSystem->Load("libMUONbase");
  gSystem->Load("libMUONrec");
  gSystem->Load("libMUONraw");
  gSystem->Load("libMUONsim");
  gSystem->Load("libFMDbase");
  gSystem->Load("libFMDsim");
  gSystem->Load("libFMDrec");
  gSystem->Load("libPMDbase");
  gSystem->Load("libPMDsim");
  gSystem->Load("libPMDrec");
  gSystem->Load("libHMPIDbase");
  gSystem->Load("libHMPIDsim");
  gSystem->Load("libHMPIDrec");
  gSystem->Load("libT0base");
  gSystem->Load("libT0sim");
  gSystem->Load("libT0rec");
  gSystem->Load("libZDCbase");
  gSystem->Load("libZDCsim");
  gSystem->Load("libZDCrec");
  gSystem->Load("libACORDE");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROsim");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libEMCALbase");
  gSystem->Load("libEMCALsim");
  gSystem->Load("libEMCALrec");
  gSystem->Load("libEMCALjet");

  // The following lines have to be commented on Darwin
  // for the moment due to cross dependencies
  gSystem->Load("libTPCbase");
  gSystem->Load("libTPCrec");
  gSystem->Load("libTPCsim");
  gSystem->Load("libTPCfast");
  gSystem->Load("libITSbase");
  gSystem->Load("libITSsim");
  gSystem->Load("libITSrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libTRDsim");
  gSystem->Load("libTRDrec");
  gSystem->Load("libTRDfast");
  gSystem->Load("libTOFbase");
  gSystem->Load("libTOFsim");
  gSystem->Load("libTOFrec");

  gSystem->Load("libAliL3ITS");
  gSystem->Load("libAliL3Src");
  gSystem->Load("libAliL3Misc");
  gSystem->Load("libAliL3Comp");
  gSystem->Load("libThread");
  gSystem->Load("libAliL3Hough");
  gSystem->Load("libANALYSIS");

  gSystem->Load("libAlieve.so");

  printf(" done.\n");
}
