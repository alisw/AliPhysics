void extrapolateXSection()
{
  TFile* file = TFile::Open("crosssection.root");
  TH1* xSection2 = dynamic_cast<TH1*> (gFile->Get("xSection2"));
  TH1* xSection15 = dynamic_cast<TH1*> (gFile->Get("xSection15"));

  gROOT->cd();

  TH1F* xSection2Ex = new TH1F("xSection2Ex", ";Npart", 1001, -0.5, 1000.5);
  TH1F* xSection15Ex = new TH1F("xSection15Ex", ";Npart", 1001, -0.5, 1000.5);

  new TCanvas;
  xSection2->Draw();
  xSection2->Fit("expo", "", "", 200, 250);
  gPad->SetLogy();

  for (Int_t i=1; i<=1000; ++i)
  {
    if (i < 250)
    {
      xSection2Ex->SetBinContent(i, xSection2->GetBinContent(i));
      xSection2Ex->SetBinError(i, xSection2->GetBinError(i));
    }
    else
      xSection2Ex->SetBinContent(i, xSection2->GetFunction("expo")->Eval(i));
  }

  new TCanvas;
  xSection2Ex->Draw();
  xSection2->SetLineColor(2);
  xSection2->Draw("SAME");
  gPad->SetLogy();

  new TCanvas;
  xSection15->Draw();
  xSection15->Fit("expo", "", "", 145, 250);
  gPad->SetLogy();

  for (Int_t i=1; i<=1000; ++i)
  {
    if (i < 145)
    {
      xSection15Ex->SetBinContent(i,xSection15->GetBinContent(i));
      xSection15Ex->SetBinError(i, xSection15->GetBinError(i));
    }
    else
      xSection15Ex->SetBinContent(i, xSection15->GetFunction("expo")->Eval(i));
  }

  new TCanvas;
  xSection15Ex->Draw();
  xSection15->SetLineColor(2);
  xSection15->Draw("SAME");
  gPad->SetLogy();

  TFile* file2 = TFile::Open("crosssectionEx.root", "RECREATE");

  xSection2Ex->Write();
  xSection15Ex->Write();

  file2->Close();
}



