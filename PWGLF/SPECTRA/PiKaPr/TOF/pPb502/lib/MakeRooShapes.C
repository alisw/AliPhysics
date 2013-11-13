{
  gSystem->Load("libRooFit");
  using namespace RooFit;
  gROOT->ProcessLine(".L RooFermiCutoff.cxx++");
  gROOT->ProcessLine(".L RooGaussianTail.cxx++");
  gROOT->ProcessLine(".L RooInverseGaussianTail.cxx++");
}
