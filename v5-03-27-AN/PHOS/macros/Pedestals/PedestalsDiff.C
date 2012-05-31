PedestalsDiff(const Int_t runNum1=0, const Int_t runNum2=0)
{
  // This macro read two root files with pedestal histograms produced
  // by Pedestals.C and compare them channel-by-channel. Difference of
  // the mean pedestals is filled into histograms
  //-
  // Arguments: runNum1, runNum2 are the run numbers which define the
  // filenames of the pedestal histogram
  //-
  // Yuri Kharlov. 15.03.2011

  TString pedhisto1 = Form("ped%d.root",runNum1);
  TString pedhisto2 = Form("ped%d.root",runNum2);
  TFile *f1 = new TFile(pedhisto1,"readonly");
  TFile *f2 = new TFile(pedhisto2,"readonly");

  TH2F  *h1m2 = (TH2F*)f1->Get("hPedHiMeanm2");
  TH2F  *h2m2 = (TH2F*)f2->Get("hPedHiMeanm2");
  TH2F  *h1m3 = (TH2F*)f1->Get("hPedHiMeanm3");
  TH2F  *h2m3 = (TH2F*)f2->Get("hPedHiMeanm3");
  TH2F  *h1m4 = (TH2F*)f1->Get("hPedHiMeanm4");
  TH2F  *h2m4 = (TH2F*)f2->Get("hPedHiMeanm4");
  h2m2->Add(h1m2,-1.);
  h2m3->Add(h1m3,-1.);
  h2m4->Add(h1m4,-1.);

  TString nameM2  =Form("Pedestal difference in runs %d and %2, module 2",runNum1,runNum2);
  h2m2->SetTitle(nameM2);
  TString nameM3  =Form("Pedestal difference in runs %d and %2, module 3",runNum1,runNum2);
  h2m3->SetTitle(nameM3);
  TString nameM4  =Form("Pedestal difference in runs %d and %2, module 4",runNum1,runNum2);
  h2m4->SetTitle(nameM4);
  
  h2m2->SetXTitle("x, cells");
  h2m3->SetXTitle("x, cells");
  h2m4->SetXTitle("x, cells");
  h2m2->SetYTitle("z, cells");
  h2m3->SetYTitle("z, cells");
  h2m4->SetYTitle("z, cells");

  h2m2->SetAxisRange(-5.,5.,"Z");
  h2m3->SetAxisRange(-5.,5.,"Z");
  h2m4->SetAxisRange(-5.,5.,"Z");

  h2m2->SetStats(0);
  h2m3->SetStats(0);
  h2m4->SetStats(0);

  h2m2->SetTitleOffset(1.0,"X");
  h2m3->SetTitleOffset(1.0,"X");
  h2m4->SetTitleOffset(1.0,"X");

  TString nameDPm2 = Form("#Delta ped (runs %d and %d) in module 2",runNum1,runNum2);
  TH1F *hDPm2 = new TH1F("hDPm2",nameDPm2,100.,-5.,5.);
  TString nameDPm3 = Form("#Delta ped (runs %d and %d) in module 3",runNum1,runNum2);
  TH1F *hDPm3 = new TH1F("hDPm3",nameDPm3,100.,-5.,5.);
  TString nameDPm4 = Form("#Delta ped (runs %d and %d) in module 4",runNum1,runNum2);
  TH1F *hDPm4 = new TH1F("hDPm4",nameDPm4,100.,-5.,5.);
  hDPm2->Sumw2();
  hDPm3->Sumw2();
  hDPm4->Sumw2();

  for (Int_t ix=0; ix<64; ix++) {
    for (Int_t iz=0; iz<56; iz++) {
      Float_t dPedm2 = h2m2->GetBinContent(ix+1,iz+1);
      Float_t dPedm3 = h2m3->GetBinContent(ix+1,iz+1);
      Float_t dPedm4 = h2m4->GetBinContent(ix+1,iz+1);
      hDPm2->Fill(dPedm2);
      hDPm3->Fill(dPedm3);
      hDPm4->Fill(dPedm4);
      if (TMath::Abs(dPedm2) > 1.)
	printf("Bad channel: m=%d, x=%2d, z=%2d, dPed=%.1f\n",2,ix,iz,dPedm2);
      if (TMath::Abs(dPedm3) > 1.)
	printf("Bad channel: m=%d, x=%2d, z=%2d, dPed=%.1f\n",3,ix,iz,dPedm3);
      if (TMath::Abs(dPedm4) > 1.)
	printf("Bad channel: m=%d, x=%2d, z=%2d, dPed=%.1f\n",4,ix,iz,dPedm4);
    }
  }

  TCanvas *c2 = new TCanvas("m2","m2",0,0,600,400);
  h2m2->Draw("colz");
  TCanvas *c3 = new TCanvas("m3","m3",0,0,600,400);
  h2m3->Draw("colz");
  TCanvas *c4 = new TCanvas("m4","m4",0,0,600,400);
  h2m4->Draw("colz");

  gStyle->SetOptStat(1110);
  TCanvas *c5 = new TCanvas("d","d",0,0,1000,400);
  c5->Divide(3,1);
  c5->cd(1);
  gPad->SetLogy();
  hDPm2->Draw("hist");
  c5->cd(2);
  gPad->SetLogy();
  hDPm3->Draw("hist");
  c5->cd(3);
  gPad->SetLogy();
  hDPm4->Draw("hist");
}
