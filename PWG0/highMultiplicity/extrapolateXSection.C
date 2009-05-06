void extrapolateXSection()
{
  TFile* file2 = TFile::Open("crosssectionEx_10TeV.root", "RECREATE");
    
  c = new TCanvas;
  TH1* base = 0;
  
  for (Int_t fileId = 0; fileId < 2; fileId++)
  {
    if (fileId == 1)
    {
      TFile* file = TFile::Open("out_phojet.root");
    }
    else
      TFile* file = TFile::Open("out_pythia.root");
    
    TH1* xSection2 = dynamic_cast<TH1*> (gFile->Get("fMult"));
    if (!base)
      base = xSection2;
    xSection2->Sumw2();
    xSection2->Scale(1.0 / xSection2->Integral());
    //TH1* xSection15 = dynamic_cast<TH1*> (gFile->Get("xSection15"));
  
    //TH1F* xSection15Ex = new TH1F("xSection15Ex", ";Npart", 1001, -0.5, 1000.5);
  
    TF1* func[3];
    Float_t lowLimit[] = { 150, 175, 200 };
    if (fileId == 1)
    {
      lowLimit[0] = 50;
      lowLimit[1] = 75;
      lowLimit[2] = 100;
    }
    
    c->cd();
    xSection2->Draw((fileId == 0) ? "" : "SAME");
    
    for (Int_t i=0; i<3; i++)
    {
      func[i] = new TF1("func", "[0]*exp([1]*x)", 0, 1000);
      func[i]->SetParameters(1, -1e-4);
    
      xSection2->Fit(func[i], "0", "", lowLimit[i], 250);
      func[i]->SetRange(lowLimit[i], 500);
      func[i]->SetLineColor(i+1);
      func[i]->Draw("SAME");
      gPad->SetLogy();
    }
    
    base->GetXaxis()->SetRangeUser(0, 500);
    base->GetYaxis()->SetRangeUser(func[2]->Eval(500), xSection2->GetMaximum());
    
    for (Int_t j=0; j<3; j++)
    {
      gROOT->cd();
      TH1F* xSection2Ex = new TH1F(Form("xSection2Ex_%d_%d", fileId, j), ";Npart", 1001, -0.5, 1000.5);
      
      for (Int_t i=1; i<=1000; ++i)
      {
        if (i < lowLimit[j])
        {
          xSection2Ex->SetBinContent(i, xSection2->GetBinContent(i));
          xSection2Ex->SetBinError(i, xSection2->GetBinError(i));
        }
        else
          xSection2Ex->SetBinContent(i, func[j]->Eval(i));
      }
      
      new TCanvas;
      xSection2Ex->Draw();
      gPad->SetLogy();
      
      file2->cd();
      xSection2Ex->Write();
    }
  }
  
  file2->Close();
}



