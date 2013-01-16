//Useful macro when you just want to open files
void loadLibraries(){
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libXMLIO.so");

    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");

    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWGUDbase.so");

    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->AddIncludePath("-I$ALICE_ROOT/PWGUD");
  gROOT->ProcessLine(".L AliAnalysisEtCuts.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisHadEtCorrections.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtCommon.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtSelector.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtSelectorPhos.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtSelectorEmcal.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtTrackMatchCorrections.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtRecEffCorrection.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEt.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtMonteCarlo.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtMonteCarloPhos.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtMonteCarloEmcal.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtReconstructed.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtReconstructedPhos.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtReconstructedEmcal.cxx+g");  
  gROOT->ProcessLine(".L AliAnalysisTaskTransverseEnergy.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEmEtMonteCarlo.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEmEtReconstructed.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisTaskTotEt.cxx+g");
}


