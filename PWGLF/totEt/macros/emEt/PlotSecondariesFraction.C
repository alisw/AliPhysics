void SetStyles(TH1 *histo,int marker, int color){
  histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  //histo->GetXaxis()->SetTitle(xtitle);
  //histo->GetYaxis()->SetTitle(ytitle);
}

void PlotSecondariesFraction(Bool_t isPhos = kFALSE){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString filename, detname;
  if(isPhos){
    detname = "PHOS";
    filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.Run139465.root";
  }
  else{
    filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root";
    detname = "EMCal";
  }
  
  TFile *f = TFile::Open(filename, "READ");
  TList *l = dynamic_cast<TList*>(f->Get("out1"));
  TH1F *fHistSecondaryEnergy = l->FindObject("fHistSecondaryEnergy");
  TH1F *fHistSecondaryChargedEnergy = l->FindObject("fHistSecondaryChargedEnergy");
  TH1F *fHistSecondaryNeutronEnergy = l->FindObject("fHistSecondaryNeutronEnergy");
  TH1F *fHistSecondaryGammaEnergy = l->FindObject("fHistSecondaryGammaEnergy");
  TH1F *fHistSecondaryElectronEnergy = l->FindObject("fHistSecondaryElectronEnergy");
  TH1F *fHistSecondaryOtherEnergy = l->FindObject("fHistSecondaryOtherEnergy");

  TH1F *fHistSecondaryChargedEnergyFraction = fHistSecondaryChargedEnergy->Clone("fHistSecondaryChargedEnergyFraction");
  TH1F *fHistSecondaryNeutronEnergyFraction = fHistSecondaryNeutronEnergy->Clone("fHistSecondaryNeutronEnergyFraction");
  TH1F *fHistSecondaryGammaEnergyFraction = fHistSecondaryGammaEnergy->Clone("fHistSecondaryGammaEnergyFraction");
  TH1F *fHistSecondaryElectronEnergyFraction = fHistSecondaryElectronEnergy->Clone("fHistSecondaryElectronEnergyFraction");
  TH1F *fHistSecondaryOtherEnergyFraction = fHistSecondaryOtherEnergy->Clone("fHistSecondaryOtherEnergyFraction");
  SetStyles(fHistSecondaryChargedEnergyFraction,20,1);
  SetStyles(fHistSecondaryNeutronEnergyFraction,21,TColor::kBlue);
  SetStyles(fHistSecondaryGammaEnergyFraction,22,TColor::kYellow);
  SetStyles(fHistSecondaryElectronEnergyFraction,23,TColor::kGreen);
  SetStyles(fHistSecondaryOtherEnergyFraction,24,TColor::kRed);


  fHistSecondaryChargedEnergyFraction->Divide(fHistSecondaryEnergy);
  fHistSecondaryNeutronEnergyFraction->Divide(fHistSecondaryEnergy);
  fHistSecondaryGammaEnergyFraction->Divide(fHistSecondaryEnergy);
  fHistSecondaryElectronEnergyFraction->Divide(fHistSecondaryEnergy);
  fHistSecondaryOtherEnergyFraction->Divide(fHistSecondaryEnergy);

  TLegend *leg = new TLegend(0.169463,0.687166,0.290268,0.911765);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->SetTextSize(0.038682);
    leg->AddEntry(fHistSecondaryChargedEnergyFraction,"Charged");
    leg->AddEntry(fHistSecondaryNeutronEnergyFraction,"Neutron");
    leg->AddEntry(fHistSecondaryGammaEnergyFraction,"Gamma");
    leg->AddEntry(fHistSecondaryElectronEnergyFraction,"Electron");
    leg->AddEntry(fHistSecondaryOtherEnergyFraction,"Other");

    TCanvas *c1 = new TCanvas("c1","Simulation",600,400);
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.03);
    c1->SetLeftMargin(0.11745);
    c1->SetBottomMargin(0.11745);
    c1->SetBorderSize(0);
    c1->SetFillColor(0);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    fHistSecondaryChargedEnergyFraction->SetMaximum(0.6);
    fHistSecondaryChargedEnergyFraction->SetMinimum(0.0);
    fHistSecondaryChargedEnergyFraction->GetXaxis()->SetTitle("Centrality bin");
    fHistSecondaryChargedEnergyFraction->GetYaxis()->SetTitle("fraction");
    fHistSecondaryChargedEnergyFraction->GetYaxis()->SetLabelSize(0.06);
    fHistSecondaryChargedEnergyFraction->GetXaxis()->SetLabelSize(0.06);
    fHistSecondaryChargedEnergyFraction->GetYaxis()->SetTitleSize(0.06);
    fHistSecondaryChargedEnergyFraction->GetXaxis()->SetTitleSize(0.06);
    fHistSecondaryChargedEnergyFraction->Draw();
    fHistSecondaryNeutronEnergyFraction->Draw("same");
    fHistSecondaryGammaEnergyFraction->Draw("same");
    fHistSecondaryElectronEnergyFraction->Draw("same");
    fHistSecondaryOtherEnergyFraction->Draw("same");
    leg->Draw();

    TString outfile = "/tmp/SecondaryFraction"+detname+".png";
    c1->SaveAs(outfile.Data());
}
