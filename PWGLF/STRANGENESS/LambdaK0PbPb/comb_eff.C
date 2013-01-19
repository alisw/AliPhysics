TH1D *combine(TH1D *pure, TH1D *injected) {
   TH1D *eff=new TH1D(*pure);
   eff->SetLineColor(2);
   eff->Add(injected);
   eff->Scale(0.5);
   for (Int_t i=1; i<16; i++) {
      Double_t c=pure->GetBinContent(i);
      Double_t e=pure->GetBinError(i);
      eff->SetBinContent(i,c);
      eff->SetBinError(i,e);       
   }
   for (Int_t i=25; i<=eff->GetNbinsX(); i++) {
      Double_t c=injected->GetBinContent(i);
      Double_t e=injected->GetBinError(i);
      eff->SetBinContent(i,c);
      eff->SetBinError(i,e);       
   }
   return eff;
}

void comb_eff(const Char_t *centr) {
  TFile *non=TFile::Open("SpectraV0CutVariations_nonInj.root");
  TFile *all=TFile::Open("SpectraV0CutVariations_all.root");

  TFile *file=TFile::Open("comb_eff.root","update");

  TString name;
  TH1D *hnon, *hall, *comb, *r;

  TCanvas *c=new TCanvas();
  c->Divide(3,2);


  non->cd();
  name="eff_K0s_nonInj_"; name+=centr;
  hnon=(TH1D*)gDirectory->Get(name.Data());
     all->cd();
     name="eff_K0s_"; name+=centr;
     hall=(TH1D*)gDirectory->Get(name.Data());
  comb=combine(hnon,hall);
  name="eff_K0s_comb_"; name+=centr;
  comb->SetName(name.Data());
  c->cd(1);
  hall->Draw(); comb->Draw("same");
  c->cd(4);  gPad->SetGridx(); gPad->SetGridy();
  r=new TH1D(*hall); r->Divide(comb); 
  r->SetMinimum(0.8); r->SetMaximum(1.2); r->Draw();
  file->cd(); comb->Write();

  non->cd(); 
  name="eff_Lambda_nonInj_"; name+=centr;
  hnon=(TH1D*)gDirectory->Get(name.Data());
     all->cd();
     name="eff_Lambda_"; name+=centr;
     hall=(TH1D*)gDirectory->Get(name.Data());
  comb=combine(hnon,hall);
  name="eff_Lambda_comb_"; name+=centr;
  comb->SetName(name.Data());
  c->cd(2);
  hall->Draw(); comb->Draw("same");
  c->cd(5);  gPad->SetGridx(); gPad->SetGridy();
  r=new TH1D(*hall); r->Divide(comb); 
  r->SetMinimum(0.8); r->SetMaximum(1.2); r->Draw();
  file->cd(); comb->Write();


  non->cd(); 
  name="eff_LambdaBar_nonInj_"; name+=centr;
  hnon=(TH1D*)gDirectory->Get(name.Data());
     all->cd();
     name="eff_LambdaBar_"; name+=centr;
     hall=(TH1D*)gDirectory->Get(name.Data());
  comb=combine(hnon,hall);
  name="eff_LambdaBar_comb_"; name+=centr;
  comb->SetName(name.Data());
  c->cd(3);
  hall->Draw(); comb->Draw("same");
  c->cd(6); gPad->SetGridx(); gPad->SetGridy();
  r=new TH1D(*hall); r->Divide(comb); 
  r->SetMinimum(0.8); r->SetMaximum(1.2); r->Draw();
  file->cd(); comb->Write();

  file->Close();
}

void Loop() {
  const Char_t *cent[]={"0005","0010","1020","2040","4060","6080","8090"};
  Int_t nCent=sizeof(cent)/sizeof(const Char_t *);
  for (Int_t i=0; i<nCent; i++) {
    cout<<cent[i]<<endl;
    comb_eff(cent[i]);
  }
}
