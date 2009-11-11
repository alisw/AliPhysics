


void makeSysErr() {

  TFile* f0 = new TFile("data/scaled25Oct09/TOTALhistosscaled-LHC09b2-0.root");
  TFile* f1 = new TFile("data/scaled25Oct09/TOTALhistosscaled-LHC09b2-1.root");
  TFile* f2 = new TFile("data/scaled25Oct09/TOTALhistosscaled-LHC09b2-2.root");
  TFile* f3 = new TFile("data/scaled25Oct09/TOTALhistosscaled-LHC09b2-3.root");

  double pyscale = (1.E6)*0.5*208*208*100/360; //seconds*lumi*Pb*Pb*acceptance

  TH2F* tte0 = (TH2F*)f0->Get("AnaElectron_hPtNPEleTTEScaled");
  TH2F* emc0 = (TH2F*)f0->Get("AnaElectron_hPtNPEleEMCALScaled");

  TH2F* tte1 = (TH2F*)f1->Get("AnaElectron_hPtNPEleTTEScaled");
  TH2F* emc1 = (TH2F*)f1->Get("AnaElectron_hPtNPEleEMCALScaled");

  TH2F* tte2 = (TH2F*)f2->Get("AnaElectron_hPtNPEleTTEScaled");
  TH2F* emc2 = (TH2F*)f2->Get("AnaElectron_hPtNPEleEMCALScaled");

  TH2F* tte3 = (TH2F*)f3->Get("AnaElectron_hPtNPEleTTEScaled");
  TH2F* emc3 = (TH2F*)f3->Get("AnaElectron_hPtNPEleEMCALScaled");

  tte0->Scale(pyscale);
  tte1->Scale(pyscale);
  tte2->Scale(pyscale);
  tte3->Scale(pyscale);
  emc0->Scale(pyscale);
  emc1->Scale(pyscale);
  emc2->Scale(pyscale);
  emc3->Scale(pyscale);

  alltte0 = (TH1F*)tte0->ProjectionX("alltte0",1,1);
  allemc0 = (TH1F*)emc0->ProjectionX("allemc0",1,1);
  alltte0->Rebin(5);
  allemc0->Rebin(5);

  alltte1 = (TH1F*)tte1->ProjectionX("alltte1",1,1);
  allemc1 = (TH1F*)emc1->ProjectionX("allemc1",1,1);
  alltte1->Rebin(5);
  allemc1->Rebin(5);

  alltte2 = (TH1F*)tte2->ProjectionX("alltte2",1,1);
  allemc2 = (TH1F*)emc2->ProjectionX("allemc2",1,1);
  alltte2->Rebin(5);
  allemc2->Rebin(5);

  alltte3 = (TH1F*)tte3->ProjectionX("alltte3",1,1);
  allemc3 = (TH1F*)emc3->ProjectionX("allemc3",1,1);
  alltte3->Rebin(5);
  allemc3->Rebin(5);
  
  Double_t tsum = 0.;
  Double_t esum = 0.;
  for(Int_t i = 1; i <= alltte0->GetNbinsX(); i++) {
    Double_t t0 = alltte0->GetBinContent(i);
    Double_t t1 = alltte1->GetBinContent(i);
    Double_t t2 = alltte2->GetBinContent(i);
    Double_t t3 = alltte3->GetBinContent(i);
    Double_t e0 = allemc0->GetBinContent(i);
    Double_t e1 = allemc1->GetBinContent(i);
    Double_t e2 = allemc2->GetBinContent(i);
    Double_t e3 = allemc3->GetBinContent(i);

    Double_t td01 = (t0 - t1)/(t0+t1+t2+t3);
    Double_t td02 = (t0 - t2)/(t0+t1+t2+t3);
    Double_t td03 = (t0 - t3)/(t0+t1+t2+t3);

    Double_t tave = (TMath::Abs(td01) + TMath::Abs(td02) + TMath::Abs(td03))/2.;
    tsum += tave;

    Double_t ed01 = (e0 - e1)/(e0+e1+e2+e3);
    Double_t ed02 = (e0 - e2)/(e0+e1+e2+e3);
    Double_t ed03 = (e0 - e3)/(e0+e1+e2+e3);

    Double_t eave = (TMath::Abs(ed01) + TMath::Abs(ed02) + TMath::Abs(ed03))/2.;
    printf("%d Average unc = %2.2f\n",i,eave);
    esum += eave;
  }
  Double_t tfinal = tsum/alltte0->GetNbinsX();
  Double_t efinal = esum/allemc0->GetNbinsX();
  printf("tfinal %f, efinal %f\n",tfinal,efinal);
  
  TCanvas *c1 = new TCanvas();
  gPad->SetLogy();
  alltte0->SetLineWidth(2);
  alltte0->Draw();

  TCanvas *c2 = new TCanvas();
  gPad->SetLogy();
  allemc0->SetLineWidth(2);
  allemc0->Draw();

  TGraphErrors* terr = new TGraphErrors();
  terr->SetName("tteErr");
  TGraphErrors* eerr = new TGraphErrors();
  eerr->SetName("emcErr");
  for(Int_t i = 1; i <= alltte0->GetNbinsX(); i++) {
    terr->SetPoint(i-1,alltte0->GetBinCenter(i),alltte0->GetBinContent(i));
    terr->SetPointError(i-1,0.,(tfinal)*alltte0->GetBinContent(i));

    eerr->SetPoint(i-1,allemc0->GetBinCenter(i),allemc0->GetBinContent(i));
    eerr->SetPointError(i-1,0.,(efinal)*allemc0->GetBinContent(i));
  }
  c1->cd();
  terr->SetFillColor(kRed-8);
  terr->Draw("3same");
  alltte0->Draw("same");

  c2->cd();
  eerr->SetFillColor(kRed-8);
  eerr->Draw("3same");
  allemc0->Draw("same");
}
