void loadlibs () 
{
  // Macro which loads the libraries needed for simulation and reconstruction
  // Possible usage: In a Root session (no AliRoot) one does
  // root [0] .x loadlibs.C
  // root [1] gAlice = new AliRun("gAlice","test")           
  // root [2] AliSimulation sim
  // root [3] sim.Run()
  // root [4] AliReconstruction rec
  // root [5] rec.Run()

  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");

  // Uncomment the following line for macosx
  // Waiting for a better solution
  //  gSystem->Load("libg2c_sh");
  gSystem->Load("libmicrocern");
  gSystem->Load("libpythia6");
  gSystem->Load("libpdf");

  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  gSystem->Load("libEGPythia6");

  gSystem->Load("libESD");
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libRAWDatasim");
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
  gSystem->Load("libMUONraw");
  gSystem->Load("libMUONsim");
  gSystem->Load("libMUONrec");
  gSystem->Load("libFMDbase");
  gSystem->Load("libFMDsim");
  gSystem->Load("libFMDrec");
  gSystem->Load("libPMDbase");
  gSystem->Load("libPMDsim");
  gSystem->Load("libPMDrec");
  gSystem->Load("libRICHbase");
  gSystem->Load("libRICHsim");
  gSystem->Load("libRICHrec");
  gSystem->Load("libSTARTbase");
  gSystem->Load("libSTARTsim");
  gSystem->Load("libSTARTrec");
  gSystem->Load("libZDCbase");
  gSystem->Load("libZDCsim");
  gSystem->Load("libZDCrec");
  gSystem->Load("libCRT");
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
}
