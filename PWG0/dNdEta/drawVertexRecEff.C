/* $Id$ */

// draws the result of AlidNdEtaVertexRecEffSelector

void drawVertexRecEff()
{
  TFile* fout = TFile::Open("vertexRecEff.root");

  TH1F* dNGen = dynamic_cast<TH1F*> fout->Get("dNGen");
  TH1F* dNRec = dynamic_cast<TH1F*> fout->Get("dNRec");

  TH1F* vtxGen = dynamic_cast<TH1F*> fout->Get("VtxGen");
  TH1F* vtxRec = dynamic_cast<TH1F*> fout->Get("VtxRec");

  // calculate ratios
  TH1F* dNRatiodN = dynamic_cast<TH1F*> (dNRec->Clone("dNRatiodN"));
  dNRatiodN->SetTitle("Ratio");
  dNRatiodN->Divide(dNGen);

  //vtxGen->Rebin(4);
  //vtxRec->Rebin(4);

  // calculate correction ratio number
  Float_t sumGen = 0;
  Float_t sumRec = 0;
  for (Int_t i=1; i<=dNGen->GetNbinsX(); ++i)
  {
    sumGen += dNGen->GetBinCenter(i) * dNGen->GetBinContent(i);
    sumRec += dNRec->GetBinCenter(i) * dNRec->GetBinContent(i);
  }
  Float_t ratio = sumRec / dNRec->Integral(1, dNGen->GetNbinsX()) * dNGen->Integral(1, dNGen->GetNbinsX()) / sumGen;

  cout << "Ratio: " << ratio << endl;

  TH1F* dNRatioVtx = dynamic_cast<TH1F*> (vtxRec->Clone("dNRatioVtx"));
  dNRatioVtx->SetTitle("Ratio");
  dNRatioVtx->Divide(vtxGen);

  TCanvas* canvas = new TCanvas("dN", "dN", 1000, 1000);
  canvas->Divide(2, 2);

  canvas->cd(1);
  dNGen->Draw();
  dNRec->SetLineColor(kRed);
  dNRec->Draw("SAME");

  canvas->cd(2);
  dNRatiodN->Draw();

  canvas->cd(3);
  vtxGen->Draw();
  vtxRec->SetLineColor(kRed);
  vtxRec->Draw("SAME");

  canvas->cd(4);
  dNRatioVtx->Draw();
}
