//
// Script to make ROOT files with LEGO information histograms
//
TH1*
process(TFile* file, const char* name, const char* opt) 
{
  TH2F*  h2d = static_cast<TH2F*>(file->Get(name));
  if (!h2d) {
    Error("process", "Couldn't get %s from %s", name, file->GetName());
    return 0;
  }
  TH1D*  h1d = h2d->ProjectionY();
  h1d->SetTitle(Form("%s", h1d->GetTitle()));
  h1d->Scale(1. / 360.);

  return h1d;
  
  return heta;
}

void
MakeLego(const Char_t* what) 
{
  TString config("FMD/scripts/ConfigInner.C");
  TString opt(what);
  if      (opt == "ITS") 
    config = "FMD/scripts/ConfigItsOnly.C";
  else if (opt == "PIPE") 
    config = "FMD/scripts/ConfigPipeOnly.C";
  else if (opt == "FMD") 
    config = "FMD/scripts/ConfigFmdOnly.C";
  else if (opt == "Nothing") 
    config = "FMD/scripts/ConfigNothing.C";
  else 
    opt = "Inner";
  
      
  cout << "Running AliRun::RunLego(" << config
       << ",180,0,180,360,0,360,0,10000,0,10000); " << endl;
  
  gAlice->RunLego(config.Data(), 180, 0, 180, 360, 0, 360, 0, 100000,
		  100000000, 0);
  
  TFile* galice   = TFile::Open("galice.root", "READ");
  TFile* output   = TFile::Open(Form("Lego_%s.root", opt.Data()),"RECREATE");
  TH1F* habso_eta = process(galice, "habso", opt.Data());
  TH1F* hradl_eta = process(galice, "hradl", opt.Data());
  TH1F* hgcm2_eta = process(galice, "hgcm2", opt.Data());
  hgcm2_eta->Draw();
  
  output->Write();
  output->Close();
  galice->Close();  
}
