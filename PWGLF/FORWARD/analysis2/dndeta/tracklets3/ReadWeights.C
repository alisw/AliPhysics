void ReadWeights()
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include -I${ALICE_PHYSICS}/include");
  gROOT->LoadMacro("AliAODTracklet.C+g");
  gROOT->LoadMacro("AliTrackletWeights.C+g");

  TFile* file = TFile::Open("weights.root","READ");
  file->ls();
  AliTrackletWeights* weights =
    static_cast<AliTrackletWeights*>(file->Get("weights"));
  TCanvas* c = new TCanvas("c","c");
  weights->Draw();
}

