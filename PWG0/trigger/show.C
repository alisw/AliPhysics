void show(const char* fileName = "trigger.root")
{
  TFile::Open(fileName);
  
  Int_t count = 0;
  while (1)
  {
    hist = (TH1*) gFile->Get(Form("fStats_%d", count));
    if (!hist)
      break;
    
    c = new TCanvas;
    hist->Draw();
    c->SaveAs(Form("trigger_%d.png", count));
    Printf("%s: %d", hist->GetTitle(), (Int_t) hist->GetEntries());
    
    count++;
  }
}
