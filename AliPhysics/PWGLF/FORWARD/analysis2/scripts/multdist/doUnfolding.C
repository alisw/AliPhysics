{
  gROOT->Macro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");
  gSystem->Load("libPWGUDbase");
  //gSystem->Load("libPWGmuondep");
  gSystem->Load("libPWGUDselectors");
  gSystem->Load("$HOME/Desktop/RooUnfold-1.0.3/libRooUnfold");
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include -I${ALICE_ROOT}/PWGPP/");

  gROOT->LoadMacro("unfoldBase.C++g");
  
  unfoldBase("900GeV_materialCorrection_unfolded.root",kBayes, "/home/caz/ALICE/AliRoot.old/PWG2/FORWARD/analysis2/finalResponse_900GeV.root", "/home/caz/ALICE/thesis/analysisFiles/900GeV/900GeV_materialCorrection.root");
  
  
}

