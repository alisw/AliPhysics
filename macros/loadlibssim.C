void loadlibssim () 
{
  // Macro which loads the libraries needed for simulation
  // Possible usage: In a Root session (no AliRoot) one does
  // root [0] .x loadlibssim.C
  // root [1] gAlice = new AliRun("gAlice","test")           
  // root [2] AliSimulation sim
  // root [3] sim.Run()

  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libProof");

  // Uncomment the following line for macosx
  // Waiting for a better solution
  //  gSystem->Load("libg2c_sh");
  gSystem->Load("libmicrocern");
  gSystem->Load("liblhapdf");
  gSystem->Load("libpythia6");

  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  gSystem->Load("libEGPythia6");

  gSystem->Load("libNet");
  gSystem->Load("libTree");
  gSystem->Load("libGui");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libAOD");
  gSystem->Load("libSTEER");
  gSystem->Load("libRAWDatasim");
  gSystem->Load("libFASTSIM");
  gSystem->Load("libEVGEN");
  gSystem->Load("libAliPythia6");

  gSystem->Load("libhijing");
  gSystem->Load("libTHijing");// AliGenHijingEventHeader needed by libZDCsim.so

  gSystem->Load("libSTRUCT");
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libPHOSbase");
  gSystem->Load("libPHOSsim");
  gSystem->Load("libMUONcore");
  gSystem->Load("libMUONmapping");
  gSystem->Load("libMUONgeometry");
  gSystem->Load("libMUONcalib");
  gSystem->Load("libMUONraw");
  gSystem->Load("libMUONtrigger");
  gSystem->Load("libMUONbase");
  gSystem->Load("libMUONsim");
  gSystem->Load("libMUONrec"); // Needed by libAliHLTMUON
  gSystem->Load("libFMDbase");
  gSystem->Load("libFMDsim");
  gSystem->Load("libPMDbase");
  gSystem->Load("libPMDsim");
  gSystem->Load("libHMPIDbase");
  gSystem->Load("libHMPIDsim");
  gSystem->Load("libT0base");
  gSystem->Load("libT0sim");
  gSystem->Load("libT0rec");
  gSystem->Load("libZDCbase");
  gSystem->Load("libZDCsim");
  gSystem->Load("libACORDEbase");
  gSystem->Load("libACORDErec");
  gSystem->Load("libACORDEsim");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libVZEROsim");
  gSystem->Load("libEMCALbase");
  gSystem->Load("libEMCALsim");

  // The following lines have to be commented on Darwin
  // for the moment due to cross dependencies
  gSystem->Load("libTPCbase");
  gSystem->Load("libTPCsim");
  gSystem->Load("libTPCrec");
  gSystem->Load("libITSbase");
  gSystem->Load("libITSsim");
  gSystem->Load("libITSrec"); // Needed by libAliHLTITS
  gSystem->Load("libTRDbase");
  gSystem->Load("libTRDsim");
  gSystem->Load("libTOFbase");
  gSystem->Load("libTOFrec");
  gSystem->Load("libTOFsim");

  gSystem->Load("libHLTbase");
  gSystem->Load("libHLTinterface");
  gSystem->Load("libHLTsim");
  gSystem->Load("libHLTrec");
}
