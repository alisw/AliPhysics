void DrawSliceProton(Double_t p, Bool_t singlePiGauss){
  //TVirtualFitter::Set
  TFile *in = TFile::Open("HFEtask.root");
  TList *qa = (TList *)in->Get("HFE_QA");

  // Make Plots for TPC
  TList *pidqa = (TList *)qa->FindObject("HFEpidQA");
  AliHFEtpcPIDqa *tpcqa = (AliHFEtpcPIDqa *)pidqa->FindObject("TPCQA");
  TH2 *hTPCsig = tpcqa->MakeSpectrumNSigma(AliHFEdetPIDqa::kBeforePID);

  // define cut model
  TF1 cutmodel("cutmodel", "[0] * TMath::Exp([1]*x) + [2]", 0, 20);
  cutmodel.SetParameter(0, -2.75);
  cutmodel.SetParameter(1, -0.8757);
  cutmodel.SetParameter(2, -0.9);

  Int_t nBinsP = 2;
  Int_t pbin = hTPCsig->GetXaxis()->FindBin(p);
  TH1 *hslice = hTPCsig->ProjectionY("hslice", pbin - nBinsP, pbin + nBinsP);

  // Functions needed
  TF1 *pi1 = new TF1("pi1", "[0] * TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5]) + [6]*TMath::Gaus(x, [7], [8])", -20, 20),
      *pi2 = new TF1("pi2", "[0] * TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5]) + [6]*TMath::Gaus(x, [7], [8])", -20, 20),
      *el1 = new TF1("el1", "[0] * TMath::Gaus(x, [1], [2])", -20, 20),
      *el2 = new TF1("el1", "[0] * TMath::Gaus(x, [1], [2])", -20, 20),
      *combined = new TF1("combined",  "[0] * TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5]) +[6] * TMath::Gaus(x, [7], [8]) + [9]*TMath::Gaus(x, [10], [11])" , -20, 20);

  // Constraints
  Double_t minele = 1e0;
  Double_t maxele = 1e3;
  Double_t minpr = 1e1; 
  Double_t maxpr = 1e4;
  Double_t minpi = 1e2;
  Double_t maxpi = 1e5;
  pi1->SetParLimits(0, minpi, maxpi); pi1->SetParLimits(3, minpi, maxpi);
  pi1->SetParLimits(2, 0.4, 1.5); pi1->SetParLimits(5, 0.4, 1.5);
  pi1->SetParLimits(1, -8,-4); pi1->SetParLimits(4, -5.5,-3);
  if(singlePiGauss){
    pi1->FixParameter(3,0);
    pi1->FixParameter(4,0);
    pi1->FixParameter(5,0);
  }
  // Set Protons
  pi1->SetParLimits(6, minpr, maxpr);
  pi1->SetParLimits(7, -7., -5.8);
  pi1->SetParLimits(8, 0.1, 1.3);
  // Set Electrons
  el1->SetParLimits(0, minele, maxele); 
  el1->SetParLimits(1, -0.9., 0.);
  el1->SetParLimits(2, 0.9, 1.5); 
  el1->SetParameter(1, -0.1);
  //el1->FixParameter(1, -1.34726e-01);
  //el1->FixParameter(2, 9.77243e-01);
  hslice->Fit(pi1, "N", "", -8, -1);
  hslice->Fit(el1, "NL", "", -.8, 3);

  // Use single fits to constrain combined fit
  combined->FixParameter(1, pi1->GetParameter(1));
  combined->FixParameter(2, pi1->GetParameter(2));
  combined->FixParameter(4, pi1->GetParameter(4));
  combined->FixParameter(5, pi1->GetParameter(5));
  if(singlePiGauss) 
    combined->FixParameter(3, 0);
  else
    combined->SetParLimits(3, minpi, maxpi);
  // Protons
  combined->FixParameter(7, pi1->GetParameter(7));
  combined->FixParameter(8, pi1->GetParameter(8));
  // Electrons
  combined->FixParameter(10, el1->GetParameter(1));
  combined->FixParameter(11, el1->GetParameter(2));
  combined->SetParLimits(0, minpi, maxpi);
  combined->SetParLimits(6, minpr, maxpr);
  combined->SetParLimits(9, minele, maxele);

  hslice->Fit(combined, "N", "", -8., 3.);

  for(Int_t ipar = 0; ipar < 9; ipar++) pi2->SetParameter(ipar, combined->GetParameter(ipar));
  for(Int_t ipar = 0; ipar < 3; ipar++) el2->SetParameter(ipar, combined->GetParameter(ipar+9));

  combined->SetLineColor(kBlack);
  el2->SetLineColor(kGreen);
  pi2->SetLineColor(kBlue);
  combined->SetLineWidth(2);
  el2->SetLineWidth(2);
  pi2->SetLineWidth(2);
  el2->SetLineStyle(2);
  pi2->SetLineStyle(2);

  hslice->SetStats(kFALSE);
  hslice->SetTitle("TPC Signal");
  hslice->GetXaxis()->SetTitle("TPC dE/dx - <dE/dx>|_{electrons} [#sigma]");

  TCanvas *output = new TCanvas("output", "Slice output");
  output->SetLogy();
  output->cd();
  hslice->Draw();
  combined->Draw("same");
  el2->Draw("same");
  pi2->Draw("same");
  
  TLegend *leg = new TLegend(0.5, 0.7, 0.89, 0.89);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hslice, "Measurement", "l");
  leg->AddEntry(combined, "Electrons + Pions", "l");
  leg->AddEntry(el2, "Electrons", "l");
  leg->AddEntry(pi2, "Pions", "l");
  leg->Draw();

  TPaveText *momentum = new TPaveText(0.6, 0.45, 0.8, 0.55, "NDC");
  momentum->SetBorderSize(0);
  momentum->SetFillStyle(0);
  momentum->AddText(Form("p = %.2fGeV/c", p));
  momentum->Draw();

  // Draw also Selection Bands
  TF1 *fBandPions = new TF1(*pi2), *fBandElectrons = new TF1(*el2);
  fBandPions->SetRange(cutmodel.Eval(p), 5);
  fBandElectrons->SetRange(cutmodel.Eval(p), 5);
  
  fBandPions->SetFillColor(kBlue);
  fBandPions->SetFillStyle(3004);
  fBandElectrons->SetFillColor(kGreen);
  fBandElectrons->SetFillStyle(3005);
  fBandPions->SetLineWidth(0);
  fBandElectrons->SetLineWidth(0);
  
  fBandElectrons->Draw("same");
  fBandPions->Draw("same");

  // Now calculate the contamination
  Double_t nPions = pi2->Integral(cutmodel.Eval(p), 5);
  Double_t nCandidates = combined->Integral(cutmodel.Eval(p), 5);
  Double_t contamination = nPions/nCandidates;
  Double_t nCandidatesHisto = hslice->Integral(hslice->GetXaxis()->FindBin(cutmodel.Eval(p)), hslice->GetXaxis()->FindBin(5));
  Double_t nPionsBefore = pi1->Integral(cutmodel.Eval(p), 5);
  Double_t contaminationHisto = nPionsBefore/nCandidatesHisto;
  printf("Contamination @ %.3fGeV/c: %f\n", p, contamination);
  printf("Contamination @ %.3fGeV/c: %f (determined from double gauss and histogram)\n", p, contaminationHisto);

  TPaveText *tcont = new TPaveText(0.45, 0.35, 0.9, 0.45, "NDC");
  tcont->SetBorderSize(0);
  tcont->SetFillStyle(0);
  tcont->AddText(Form("contamination: %f", contamination));
  tcont->Draw();

  // Save Canvas
  TFile *out = new TFile(Form("slice%dpr.root", ((Int_t)(100*p))), "RECREATE");
  out->cd();
  output->Write();
  out->Close(); delete out;
}
