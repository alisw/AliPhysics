#include <THtml.h>
#include <TROOT.h>
#include <TSystem.h>

void mkhtml (char *macro=0, Int_t force=0) {
// to run this macro, you must have the correct .rootrc file
// in your galice directory.
// The gAlice classes summary documentation go to directory html
// The gAlice classes source  documentation go to directory html/src
// The example macros documentation go to directory html/examples
   
  // gROOT->LoadMacro("loadlibs.C");
  // loadlibs();
  THtml html;
  html.SetProductName("AliRoot");
  if(macro) {
    gROOT->LoadMacro(macro);
    html.Convert(macro,"Example Macro");
  } else {
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations

    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libANALYSIScalib.so");
    gSystem->Load("libEventMixing.so");

    gSystem->Load("libPWG0base.so");
    gSystem->Load("libPWG0dep.so");
    gSystem->Load("libPWG0selectors.so");

    gSystem->Load("libTENDER.so");
    gSystem->Load("libPWG1.so");

    gSystem->Load("libCORRFW.so");
    gSystem->Load("libPWG2AOD.so");
    gSystem->Load("libPWG2ebye.so");
    gSystem->Load("libPWG2evchar.so");
    gSystem->Load("libPWG2femtoscopy.so");
    gSystem->Load("libPWG2femtoscopyUser.so");
    gSystem->Load("libPWG2flowCommon.so");
    gSystem->Load("libPWG2flowTasks.so");
    gSystem->Load("libPWG2forward.so");
    gSystem->Load("libPWG2kink.so");
    gSystem->Load("libPWG2resonances.so");
    gSystem->Load("libPWG2spectra.so");
    gSystem->Load("libPWG2unicor.so");

    gSystem->Load("libPWG3base.so");
    gSystem->Load("libPWG3hfe.so");
    gSystem->Load("libPWG3muondep.so");
    gSystem->Load("libPWG3muon.so");
    gSystem->Load("libPWG3vertexingHF.so");

    gSystem->Load("libJETAN.so");
    gSystem->Load("libPWG4CaloCalib.so");
    gSystem->Load("libPWG4GammaConv.so");
    gSystem->Load("libPWG4JetTasks.so");
    gSystem->Load("libPWG4omega3pi.so");
    gSystem->Load("libPWG4PartCorrBase.so");
    gSystem->Load("libPWG4PartCorrDep.so");

    // EVE
    gSystem->Load("libEve.so");
    gSystem->Load("libEveBase.so");
    gSystem->Load("libEveDet.so");
    gSystem->Load("libEveHLT.so");

    html.MakeAll(force,"[A-Z]*");
  }
}
