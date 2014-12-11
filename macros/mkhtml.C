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
    gSystem->Load("liblhapdf");      // Parton density functions
    gSystem->Load("libEGPythia6");   // TGenerator interface
    gSystem->Load("libpythia6");     // Pythia
    gSystem->Load("libAliPythia6");  // ALICE specific implementations

    // ANALYSIS
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libANALYSIScalib");
    gSystem->Load("libESDfilter");
    gSystem->Load("libEventMixing");
    gSystem->Load("libTender");
    gSystem->Load("libTenderSupplies");

    // CORRFW
    gSystem->Load("libCORRFW"); // 

    // PWG
    gSystem->Load("libPWGTools");
    gSystem->Load("libPWGGlauber");
    gSystem->Load("libPWGflowBase");
    gSystem->Load("libPWGflowTasks");
    gSystem->Load("libPWGmuon");
    gSystem->Load("libPWGmuondep");
    gSystem->Load("libPWGEMCAL");
    gSystem->Load("libPWGCaloTrackCorrBase");

    // PWGCF
    gSystem->Load("libPWGCFCorrelationsBase");
    gSystem->Load("libPWGCFCorrelationsDPhi");
    gSystem->Load("libPWGCFCorrelationsJCORRAN");
    gSystem->Load("libPWGCFChaoticity");
    gSystem->Load("libPWGCFFEMTOSCOPYAOD");
    gSystem->Load("libPWGCFfemtoscopy");
    gSystem->Load("libPWGCFfemtoscopyUser");
    gSystem->Load("libPWGCFunicor");
    gSystem->Load("libPWGCFebye");
    //PH    gSystem->Load("libPWGCFK0Analysis");

    // PWGDQ
    gSystem->Load("libPWGDQbase");
    gSystem->Load("libPWGDQdielectron");

    // PWGGA
    gSystem->Load("libPWGGACaloTasks");
    gSystem->Load("libPWGGACaloTrackCorrelations");
    gSystem->Load("libPWGGAEMCALTasks");
    gSystem->Load("libPWGGAGammaConv");
    gSystem->Load("libPWGGAPHOSTasks");

    // PWGHF
    gSystem->Load("libPWGHFbase");
    gSystem->Load("libPWGHFhfe");
    gSystem->Load("libPWGHFcorrelationHF");
    gSystem->Load("libPWGHFvertexingHF");

    // PWGJE
    gSystem->Load("libJETAN");
    // gSystem->Load("libPWGJE");
    // gSystem->Load("libPWGJEEMCALJetTasks");

    // PWGLF
    gSystem->Load("libPWGLFSTRANGENESS");
    gSystem->Load("libPWGLFforward2");
    gSystem->Load("libPWGLFresonances");
    gSystem->Load("libPWGLFrsnextra");
    gSystem->Load("libPWGLFspectra");
    // gSystem->Load("libPWGLFtotEt")

    // PWGPP
    gSystem->Load("libPWGPP");
    gSystem->Load("libPWGPPMUONdep");
    gSystem->Load("libPWGPPMUONlite");
    gSystem->Load("libPWGPPevchar");

    // PWGUD

    gSystem->Load("libPWGUDbase");
    gSystem->Load("libPWGUDFP");
    gSystem->Load("libPWGUDdiffractive");
    gSystem->Load("libPWGUDselectors");

    // ITS/UPGRADE
    gSystem->Load("libITSUpgradeBase");
    gSystem->Load("libITSUpgradeSim");
    gSystem->Load("libITSUpgradeRec");

    // EVE
    gSystem->Load("libEve");
    gSystem->Load("libEveBase");
    gSystem->Load("libEveDet");
    gSystem->Load("libEveHLT");

    html.MakeAll(force,"[A-Z]*");
  }
}
