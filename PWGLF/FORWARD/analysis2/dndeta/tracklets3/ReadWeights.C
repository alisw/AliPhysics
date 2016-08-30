void ReadWeights(const char* filename="weights.root")
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include -I${ALICE_PHYSICS}/include");
  gROOT->LoadMacro("AliAODTracklet.C+g");
  gROOT->LoadMacro("AliTrackletWeights.C+g");

  TFile* file = TFile::Open(filename,"READ");
  file->ls();
  AliTrackletBaseWeights* weights =
    static_cast<AliTrackletBaseWeights*>(file->Get("weights"));
  TCanvas* c = new TCanvas("c","c");
  weights->Draw();
}

