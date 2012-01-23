void DrawContamination(const char *fname = "HFEtask.root"){
  //TVirtualFitter::Set
  TFile *in = TFile::Open(fname);
  TList *qa = (TList *)in->Get("HFE_QA");

  // Make Plots for TPC
  TList *pidqa = (TList *)qa->FindObject("PIDqa");
  TList *tpc = (TList *)pidqa->FindObject("list_fQAhistosTPC");
  // Create histograms by projecting the THnSparse
  THnSparseF *hTPC = (THnSparseF *)tpc->FindObject("fHistTPCsignal");
  hTPC->GetAxis(4)->SetRange(1,1);
  TH2 *hTPCsig = hTPC->Projection(3,1);
  hTPCsig->SetName("hTPCsig");

  // define cut model
  TF1 cutmodel("cutmodel", "[0] * TMath::Exp([1]*x) + [2]", 0, 20);
  cutmodel.SetParameter(0, -2.75);
  cutmodel.SetParameter(1, -0.8757);
  cutmodel.SetParameter(2, -0.9);

  Int_t offset = hTPCsig->GetXaxis()->FindBin(0.35);
  Int_t nPoints = hTPCsig->GetXaxis()->FindBin(4.) - offset;
  Int_t points = 0;
  TArrayD ap(nPoints), apErr(nPoints), acont(nPoints), aeff(nPoints), aemean(nPoints), aesigma(nPoints), aemeanErr(nPoints), aeSigmaErr(nPoints);
  Double_t p = 0., pmin = 0, pmax = 0, lowerbound = 0.;
  TH1 *hde = NULL;
  TCanvas outfit("outfit","Output of the fit", 640, 480);
  TCanvas singlefit("singlefit","Output of the fit", 640, 480);
  for(Int_t ipbin = hTPCsig->GetXaxis()->FindBin(0.5); ipbin < hTPCsig->GetXaxis()->FindBin(4.);ipbin+=12){
    p = hTPCsig->GetXaxis()->GetBinCenter(ipbin);
    lowerbound = cutmodel.Eval(p);
    printf("ipbin = %d, p = %f, Lower limit %f\n", ipbin, p, lowerbound);
    hde = hTPCsig->ProjectionY("py", ipbin-6, ipbin+6);
    pmin =  hTPCsig->GetXaxis()->GetBinLowEdge(ipbin-6);
    pmax =  hTPCsig->GetXaxis()->GetBinUpEdge(ipbin+6);

    // First fit gaus for electrons and pions 
    TF1 gpsingle("gpsingle", "[0]*TMath::Gaus(x, [1], [2])", -7, -1.5),
        gesingle("gpsingle", "[0]*TMath::Gaus(x, [1], [2])", -1.5, 3);

    gesingle.SetParLimits(0, 1, 1e9); gpsingle.SetParLimits(0, 1, 1e9);
    gesingle.SetParLimits(1, -0.7, 1.2); gpsingle.SetParLimits(1, -6, -4);
    gesingle.SetParLimits(2, 0.95, 1.35); gpsingle.SetParLimits(2, 0.75, 1.35);

    Double_t lowerlimitpi = -6, upperlimitpi = -2;
    if(p > 2.7) {lowerlimitpi = -5; upperlimitpi = -1; gpsingle.SetParLimits(-5, -2)}
    hde->Fit(&gpsingle, "", "N", lowerlimitpi, upperlimitpi);
    hde->Fit(&gesingle, "", "N", -1, 2);

    singlefit.cd();
    singlefit.SetLogy();
    hde->Draw();
    gesingle.SetLineColor(kBlue);
    gesingle.Draw("lsame");
    gpsingle.SetLineColor(kGreen);
    gpsingle.Draw("lsame");
    singlefit.SaveAs(Form("singlegaus_input%dGeV.png", static_cast<Int_t>(p*100)));

    printf("Model restricted, fit with double gauss\n");
    // Define 2 gaussian model: 1st gauss electron, 2nd gauss pion
    TF1 fitfun("2gausmodel", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x,[4],[5])", -20, 20);
    Double_t errmp = 0.1 * TMath::Abs(gpsingle.GetParameter(1)),
             errme = 0.2 * TMath::Abs(gesingle.GetParameter(1)),
             errsp = 0.1 * gpsingle.GetParameter(2),
             errse = 0.1 * gesingle.GetParameter(2);
    printf("Limits:\n=============================================\n");
    printf("Electrons: Mean:\t\t%.3f, Sigma:\t\t%.3f\n", errme, errse);
    printf("Mean       Lower:\t\t%.3f, Upper:\t\t%.3f\n", gesingle.GetParameter(1) - errme, gesingle.GetParameter(1) + errme);
    printf("Sigma      Lower:\t\t%.3f, Upper:\t\t%.3f\n", gesingle.GetParameter(2) - errse, gesingle.GetParameter(2) + errse);
    printf("Pions:     Mean:\t\t%.3f, Sigma:\t\t%.3f\n", errmp, errsp);
    printf("Mean       Lower:\t\t%.3f, Upper:\t\t%.3f\n", gpsingle.GetParameter(1) - errmp, gpsingle.GetParameter(1) + errmp);
    printf("Mean       Lower:\t\t%.3f, Upper:\t\t%.3f\n", gpsingle.GetParameter(2) - errsp, gpsingle.GetParameter(2) + errsp);

    fitfun.SetParLimits(0, 0.1, 1e9);
    fitfun.SetParLimits(3, 0.1, 1e9);
    fitfun.SetParLimits(1, gesingle.GetParameter(1) - errme, gesingle.GetParameter(1) + errme);
    fitfun.SetParLimits(4, gpsingle.GetParameter(1) - errmp, gpsingle.GetParameter(1) + errmp);
    fitfun.SetParLimits(2, gesingle.GetParameter(2) - errse, gesingle.GetParameter(2) + errse);
    fitfun.SetParLimits(5, gpsingle.GetParameter(2) - errsp, gpsingle.GetParameter(2) + errsp);
    fitfun.SetParameter(1, gesingle.GetParameter(1));
    fitfun.SetParameter(4, gpsingle.GetParameter(1));
    fitfun.SetParameter(2, gesingle.GetParameter(2));
    fitfun.SetParameter(5, gpsingle.GetParameter(2));
    // Fit the histogram
    hde->Fit(&fitfun, "","", lowerlimitpi, 3);
    // Save monitoring plot
    outfit.cd();
    outfit.SetLogy();
    outfit.Clear();
    hde->SetTitle();
    hde->GetYaxis()->SetTitle("Yield");
    hde->SetStats(kFALSE);
    hde->Draw();
    TPaveText plabel(0.1, 0.7, 0.3, 0.89, "NDC");
    plabel.AddText(Form("p = %0.2fGeV/c", p));
    plabel.SetBorderSize(0);
    plabel.SetFillStyle(0);
    plabel.Draw();
    TPaveText fitpar(0.6, 0.5, 0.89, 0.89, "NDC"); fitpar.AddText("Fit Patameters");
    fitpar.SetBorderSize(0); fitpar.SetFillStyle(0);
    for(Int_t ipar = 0; ipar < 6; ipar++) fitpar.AddText(Form("p%d: %0.3e #pm %0.3e", ipar, fitfun.GetParameter(ipar), fitfun.GetParError(ipar)));
    fitpar.Draw();
    outfit.SaveAs(Form("contaminationFit%dGeV.png", static_cast<Int_t>(p*100)));
    
    // Make single gaus functions, integrate them within the electron selection band
    TF1 gausele("gausele", "[0]*TMath::Gaus(x, [1], [2])", -20, 20), 
        gauspi("gausele", "[0]*TMath::Gaus(x, [1], [2])", -20, 20);
    for(Int_t ipar = 0; ipar < 3; ipar++){
      gausele.SetParameter(ipar, fitfun.GetParameter(ipar));
      gauspi.SetParameter(ipar, fitfun.GetParameter(ipar+3));
    }
    Double_t cele = gausele.Integral(lowerbound, 2.);
    Double_t celt = gausele.Integral(-10, 10);
    Double_t cpi  = gauspi.Integral(lowerbound, 2.);
    Double_t ct   = fitfun.Integral(lowerbound, 2);
    Double_t cont = cpi / ct;
    Double_t eff = celt ? cele/celt : 0;
    printf("Contamination level at p = %0.2f: %f\n", p, cont);
    printf("Efficiency level at p = %0.2f: %f\n", p, eff);

    ap[points] = p;
    apErr[points] = (pmax - pmin)/2;
    acont[points] = cont;
    aeff[points] = eff;
    aemean[points] = fitfun.GetParameter(1); 
    aesigma[points] = fitfun.GetParameter(2);
    aemeanErr[points] = fitfun.GetParError(1);
    aeSigmaErr[points] = fitfun.GetParError(2);
    points++;
    delete hde;
  }
  TGraph *contamination = new TGraph(points);
  TGraph *efficiency = new TGraph(points);
  TGraphErrors *electronMean = new TGraphErrors(points);
  TGraphErrors *electronSigma = new TGraphErrors(points);
  for(Int_t ip = 0; ip < points; ip++){
    contamination->SetPoint(ip, ap[ip], acont[ip]);
    efficiency->SetPoint(ip, ap[ip], aeff[ip]);
    electronMean->SetPoint(ip, ap[ip], aemean[ip]);
    electronMean->SetPointError(ip, apErr[ip], aemeanErr[ip]);
    electronSigma->SetPoint(ip, ap[ip], aesigma[ip]);
    electronSigma->SetPointError(ip, apErr[ip], aeSigmaErr[ip]);
  }

  contamination->SetName("contamination");
  efficiency->SetName("efficiency");
 
  TFile *outfile = new TFile("Fitresults.root", "recreate");
  outfile->cd();
  efficiency->Write();
  contamination->Write();
  gROOT->cd();

  // Draw Plots
  TCanvas *cCont = new TCanvas("cContamination", "Contamination", 800, 600);
  cCont->cd();
  contamination->SetMarkerColor(kBlue);
  contamination->SetMarkerStyle(22);
  contamination->GetXaxis()->SetTitle("p_{T} / GeV/c");
  contamination->GetYaxis()->SetTitle("contamination");
  contamination->Draw("ap");
  TCanvas *cEff = new TCanvas("cEfficiency", "Efficiency", 800, 600);
  cEff->cd();
  efficiency->SetMarkerColor(kBlue);
  efficiency->SetMarkerStyle(22);
  efficiency->GetXaxis()->SetTitle("p_{T} / GeV/c");
  efficiency->GetYaxis()->SetTitle("efficiency");
  efficiency->Draw("ap");
  outfile->Close();
  TCanvas *cElectrons = new TCanvas("cElectrons", "Electron Mean and Sigma");
  cElectrons->cd();
  electronMean->SetMarkerStyle(22);
  electronMean->SetMarkerColor(kBlue);
  electronMean->SetLineColor(kBlue);
  electronMean->SetLineWidth(1);
  electronMean->GetXaxis()->SetTitle("p [GeV/c]");
  electronMean->GetYaxis()->SetTitle("Electron Mean/Sigma [n#sigma]");
  electronMean->GetYaxis()->SetRangeUser(-1., 2);
  electronSigma->SetMarkerStyle(22);
  electronSigma->SetMarkerColor(kRed);
  electronSigma->SetLineColor(kRed);
  electronSigma->SetLineWidth(1);
  electronSigma->GetXaxis()->SetTitle("p [GeV/c]");
  electronSigma->GetYaxis()->SetTitle("Electron Mean/Sigma [n#sigma]");
  electronSigma->GetYaxis()->SetRangeUser(-1., 2.);
  electronMean->Draw("ape");
  electronSigma->Draw("epsame");
  TLegend *leg = new TLegend(0.7, 0.75, 0.89, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(electronMean, "Mean", "p");
  leg->AddEntry(electronSigma, "Sigma", "p");
  leg->Draw();
  delete outfile;
}
