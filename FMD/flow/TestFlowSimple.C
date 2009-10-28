Double_t GeneratePhi(Double_t psiR, Double_t v1, Double_t v2)
{
  Double_t phi  =  gRandom->Uniform(0, 2 * TMath::Pi());
  Double_t rel  =  phi - psiR;
  Double_t dphi =  -2 * TMath::Sin(rel) * v1;
  dphi          -= TMath::Sin(2 * rel) * v2;
  phi           += dphi;
  return phi;
}

void TestFlowSimple()
{
  gSystem->Load("libFMDflow.so");
  AliFMDFlowAxis      axis(10, -5, 5);
  AliFMDFlowBinned1D  flow("flow", "analysed", 2, 1, axis);
  TArrayD             phis(20000);
  TArrayD             xs(20000);

  for (Int_t i = 0; i < 100; i++) { 
    Double_t psiR = gRandom->Uniform(0, 2*TMath::Pi());
    Int_t    nObs = gRandom->Integer(20000);
    
    for (Int_t j = 0; j < nObs; j++) { 
      xs[j]   = gRandom->Uniform(-5, 5);
      phis[j] = GeneratePhi(psiR, 0, 0.05);
    }
    flow.Event(nObs, phis.fArray, xs.fArray);
  }
  flow.Draw("ht");
}

