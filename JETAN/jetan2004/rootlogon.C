// $Id$

{
  printf("\nWELCOME to the magic world of jets\n\n"); 

  gSystem->SetIncludePath("-Wno-deprecated -DALICEINTERFACE -I$ALICE_ROOT/include -I$ALICE_ROOT/JETAN -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/HLT/ITS");

  Int_t saveErrIgLevel=gErrorIgnoreLevel;
  gErrorIgnoreLevel=kFatal;

  gSystem->Load("libAliHLTITS.so");

  gSystem->Load("libTkJetFinder.so");
  gSystem->Load("libJetFinder.so");
  gSystem->Load("libJETAN.so");

  gSystem->Load("testJets_C.so");
  gSystem->Load("anaAliJets_C.so");
  gSystem->Load("anaTracks_C.so");
  gSystem->Load("createEvents_C.so");
  gSystem->Load("findJets_C.so");
  gSystem->Load("findJetsMixed_C.so");

  gErrorIgnoreLevel=saveErrIgLevel;
  gErrorIgnoreLevel=kWarning;
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.036,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetTitleOffset(1.3,"XYZ");
  gStyle->SetLabelSize(0.036,"XYZ");
  gStyle->SetTitleColor(1,"XYZ");
  gStyle->SetHistLineWidth(3);
  //gStyle->SetLabelSize(0.06,"XYZ");
  //gStyle->SetLabelOffset(0.009,"XYZ");

  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat("NEMR");
  gROOT->ForceStyle();
}
