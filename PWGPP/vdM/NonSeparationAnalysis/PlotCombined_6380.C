// -*- C++ -*-

void PlotCombined_6380()
{
  PlotCombinedSingle(0);
  gPad->SaveAs("pdf/6380/6380_b-by-b.pdf(");
  PlotCombinedSingle(1);
  gPad->SaveAs("pdf/6380/6380_b-by-b.pdf");
  PlotCombinedSingle(2);
  gPad->SaveAs("pdf/6380/6380_b-by-b.pdf)");
}
void PlotCombinedSingle(Int_t iType=0)
{
  AliTriggerBCMask mask("B", "344H1L39H1L39H1L39H1L362H1L39H1L39H1L39H1L39H1L39H1L39H1L39H1L39H1L39H1L370H1L39H1L39H1L39H1L1498H1L39H1L39H1L39H1L266H");
  Int_t n = mask.GetNUnmaskedBCs();

  const char* scan="Scan1_Scan2";
  const char* types[3] = {
    "",
    "_fix_alpha_XZ",
    "_fix_alpha_XZ_fix_alpha_YZ"
  };
  const char* type= types[iType];
  TH2 *hR = new TH2D("hR", Form("fill 6380 %s%s;BC;R", scan, type), n,0,n, 40, 0.99, 1.01);
  TH1 *h1R = new TH1D("h1R", Form("fill 6380 %s%s;BC;R", scan, type), n,0,n);
  TH1 *hChi2NDF = new TH1D("hChi2NDF", Form("fill 6380 %s%s;BC;#chi^{2}/ndf", scan, type), n,0,n);

  TVectorD *R = NULL;
  Double_t *r= NULL;
  TLine *line = new TLine;
  TVectorD r0(n);
  TParameter<Double_t> *chi2=NULL, *ndf=NULL;
  for (Int_t bc=0,counter=0; bc<3564; ++bc) {
    if (mask.GetMask(bc))
      continue;

    hR->GetXaxis()->SetBinLabel(1+counter, Form("%d", bc));
    h1R->GetXaxis()->SetBinLabel(1+counter, Form("%d", bc));
    hChi2NDF->GetXaxis()->SetBinLabel(1+counter, Form("%d", bc));
    TFile::Open(Form("pdf/6380/fit_%s%s_bcid%d.pdf_canvas.root", scan, type, bc));
    gFile->GetObject("R", R);
    if (R) {
      Printf("%p %d", R, bc);
      r = R->GetMatrixArray();
      Printf("%p", r);
      for (Int_t j=1; j<R->GetNoElements(); ++j) {
        hR->Fill(counter, r[j]);
      }
      r0(counter) = r[0];
      if (R->GetNoElements() > 1) {
        h1R->SetBinContent(1+counter, r[0]);
        h1R->SetBinError(1+counter, TMath::RMS(R->GetNoElements()-1, r+1));
      } else {
        h1R->SetMinimum(0.99);
      }
      delete R;
    }
    gFile->Close();


    TFile::Open(Form("root/6380/fit_%s%s_bcid%d.root", scan, type, bc));
    gFile->GetObject("chi2", chi2);
    gFile->GetObject("ndf", ndf);
    if (chi2 && ndf) {
      Printf("%4d: %.2f", bc, chi2->GetVal()/ndf->GetVal());
      hChi2NDF->SetBinContent(1+counter, chi2->GetVal()/ndf->GetVal());
    }
    gFile->Close();
    ++counter;
  }

  gStyle->SetOptFit(111);
  gStyle->SetFitFormat(".5f");
  TCanvas *c1 = new TCanvas;
  c1->Divide(2,2);
  c1->cd(1)->SetGrid(1,1);
  hR->SetLabelSize(0.05, "Y");
  hR->SetLabelSize(0.07, "X");
//  hR->SetStats(0);
  hR->Draw("COLZ");
  line->SetLineWidth(3);
  line->SetLineColor(kMagenta);
  for (Int_t i=0; i<n; ++i) {
    line->DrawLine(i, r0(i), i+1, r0(i));
  }
  hR->Fit("pol0");
  c1->cd(3)->SetGrid();
  hChi2NDF->SetLabelSize(0.05, "Y");
  hChi2NDF->SetLabelSize(0.07, "X");
  hChi2NDF->SetStats(0);
  hChi2NDF->SetMinimum(0);
  hChi2NDF->Draw();
  c1->cd(2);
  hR->ProjectionY()->Draw();
  c1->cd(4);
  h1R->SetLabelSize(0.07, "X");
  h1R->Draw();
  h1R->Fit("pol0", "");
  c1->cd();
}
