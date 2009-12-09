void show(const char* fileName = "trigger.root")
{
  TH1* fStats[100];

  TFile::Open(fileName);
  
  Int_t count = 0;
  while (1)
  {
    hist = (TH1*) gFile->Get(Form("fStats_%d", count));
    if (!hist)
      break;
    
    c = new TCanvas;
    hist->Rebin(10);
    hist->Draw();
    c->SaveAs(Form("trigger_%d.png", count));
    Printf("%s: %d", hist->GetTitle(), (Int_t) hist->GetEntries());
    
    fStats[count] = hist;
    count++;
  }

  Int_t base = 1;
  if (fStats[0])
    base = (Int_t) fStats[0]->Integral(0, fStats[0]->GetNbinsX()+1);

  for (Int_t i=0; i<count; i++)
    if (fStats[i])
    {
      c->cd(i+1);
      fStats[i]->Draw();
      Printf("%s: %d triggers | %f %% of all triggered", fStats[i]->GetTitle(), (UInt_t) fStats[i]->Integral(0, fStats[i]->GetNbinsX()+1), 100.0 * fStats[i]->Integral(0, fStats[i]->GetNbinsX()+1) / base);
    }

}
