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

    // ANALYSIS
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libANALYSIScalib.so");
    gSystem->Load("libEventMixing.so");
    gSystem->Load("libTENDER.so");
    gSystem->Load("libTENDERSupplies.so");

    // CORRFW
    gSystem->Load("libCORRFW.so"); // 

    // PWG
    gSystem->Load("libPWGTools.so");
    gSystem->Load("libPWGGlauber.so");
    gSystem->Load("libPWGflowBase.so");
    gSystem->Load("libPWGflowTasks.so");
    gSystem->Load("libPWGmuon.so");
    gSystem->Load("libPWGmuondep.so");
    gSystem->Load("libPWGEMCAL.so");
    gSystem->Load("libPWGCaloTrackCorrBase.so");

    // PWGCF
    gSystem->Load("libPWGCFCorrelationsBase.so");
    gSystem->Load("libPWGCFCorrelationsDPhi.so");
    gSystem->Load("libPWGCFCorrelationsJCORRAN.so");
    gSystem->Load("libPWGCFChaoticity.so");
    gSystem->Load("libPWGCFFEMTOSCOPYAOD.so");
    gSystem->Load("libPWGCFfemtoscopy.so");
    gSystem->Load("libPWGCFfemtoscopyUser.so");
    gSystem->Load("libPWGCFunicor.so");
    gSystem->Load("libPWGCFebye.so");
    //PH    gSystem->Load("libPWGCFK0Analysis.so");

    // PWGDQ
    gSystem->Load("libPWGDQbase.so");
    gSystem->Load("libPWGDQdielectron.so");

    // PWGGA
    gSystem->Load("libPWGGACaloTasks.so");
    gSystem->Load("libPWGGACaloTrackCorrelations.so");
    gSystem->Load("libPWGGAEMCALTasks.so");
    gSystem->Load("libPWGGAGammaConv.so");
    gSystem->Load("libPWGGAPHOSTasks.so");

    // PWGHF
    gSystem->Load("libPWGHFbase.so");
    gSystem->Load("libPWGHFhfe.so");
    gSystem->Load("libPWGHFcorrelationHF.so");
    gSystem->Load("libPWGHFvertexingHF.so");

    // PWGJE
    gSystem->Load("libJETAN.so");
    // gSystem->Load("libPWGJE.so");
    // gSystem->Load("libPWGJEEMCALJetTasks.so");

    // PWGLF
    gSystem->Load("libPWGLFSTRANGENESS.so");
    gSystem->Load("libPWGLFforward.so");
    gSystem->Load("libPWGLFforward2.so");
    gSystem->Load("libPWGLFresonances.so");
    gSystem->Load("libPWGLFrsnextra.so");
    gSystem->Load("libPWGLFspectra.so");
    // gSystem->Load("libPWGLFtotEt.so")

    // PWGPP
    gSystem->Load("libPWGPP.so");
    gSystem->Load("libPWGPPMUONdep.so");
    gSystem->Load("libPWGPPMUONlite.so");
    gSystem->Load("libPWGPPevchar.so");

    // PWGUD

    gSystem->Load("libPWGUDbase.so");
    gSystem->Load("libPWGUDFP.so");
    gSystem->Load("libPWGUDdiffractive.so");
    gSystem->Load("libPWGUDselectors.so");

    // ITS/UPGRADE
    gSystem->Load("libITSUpgradeBase.so");
    gSystem->Load("libITSUpgradeSim.so");
    gSystem->Load("libITSUpgradeRec.so");

    // EVE
    gSystem->Load("libEve.so");
    gSystem->Load("libEveBase.so");
    gSystem->Load("libEveDet.so");
    gSystem->Load("libEveHLT.so");

    html.MakeAll(force,"[A-Z]*");
  }
}
