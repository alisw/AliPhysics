void SetStyles(TH1 *histo,int marker, int color){
  histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  //histo->GetXaxis()->SetTitle(xtitle);
  //histo->GetYaxis()->SetTitle(ytitle);
}
void PlotEmEtDistribution(Bool_t isPhos = kFALSE){
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
  TH1F *fHistAllEnergy = l->FindObject("fHistAllEnergy");
  TH1F *fHistSignalEnergy = l->FindObject("fHistSignalEnergy");
  TH1F *fHistNeutronEnergy = l->FindObject("fHistNeutronEnergy");
  TH1F *fHistKaonEnergy = l->FindObject("fHistKaonEnergy");
  TH1F *fHistHadronEnergy = l->FindObject("fHistHadronEnergy");
  TH1F *fHistSecondaryEnergy = l->FindObject("fHistSecondaryEnergy");
  TH1F *fHistSignalEnergyFraction = fHistSignalEnergy->Clone("fHistSignalEnergyFraction");
  TH1F *fHistNeutronEnergyFraction = fHistNeutronEnergy->Clone("fHistNeutronEnergyFraction");
  TH1F *fHistKaonEnergyFraction = fHistKaonEnergy->Clone("fHistKaonEnergyFraction");
  TH1F *fHistHadronEnergyFraction = fHistHadronEnergy->Clone("fHistHadronEnergyFraction");
  TH1F *fHistSecondaryEnergyFraction = fHistSecondaryEnergy->Clone("fHistSecondaryEnergyFraction");
  SetStyles(fHistSignalEnergyFraction,20,TColor::kRed);
  SetStyles(fHistNeutronEnergyFraction,21,TColor::kYellow);
  SetStyles(fHistSecondaryEnergyFraction,22,TColor::kGreen);
  SetStyles(fHistHadronEnergyFraction,23,1);
  SetStyles(fHistKaonEnergyFraction,24,TColor::kBlue);
 fHistSignalEnergyFraction->Divide(fHistAllEnergy);
 fHistNeutronEnergyFraction->Divide(fHistAllEnergy);
 fHistSecondaryEnergyFraction->Divide(fHistAllEnergy);
 fHistHadronEnergyFraction->Divide(fHistAllEnergy);
 fHistKaonEnergyFraction->Divide(fHistAllEnergy);
 TLegend *leg = new TLegend(0.67953,0.36631,0.800336,0.590909);
 leg->SetFillStyle(0);
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->SetTextSize(0.03);
 leg->SetTextSize(0.038682);
 leg->AddEntry(fHistSignalEnergyFraction,"Signal");
 leg->AddEntry(fHistHadronEnergyFraction,"Hadron");
 leg->AddEntry(fHistSecondaryEnergyFraction,"Secondary");
 leg->AddEntry(fHistKaonEnergyFraction,"Kaon");
 leg->AddEntry(fHistNeutronEnergyFraction,"Neutron");

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
    fHistSignalEnergyFraction->SetMaximum(0.6);
    fHistSignalEnergyFraction->SetMinimum(0.0);
    fHistSignalEnergyFraction->GetXaxis()->SetTitle("Centrality bin");
    fHistSignalEnergyFraction->GetYaxis()->SetTitle("fraction");
    fHistSignalEnergyFraction->GetYaxis()->SetLabelSize(0.06);
    fHistSignalEnergyFraction->GetXaxis()->SetLabelSize(0.06);
    fHistSignalEnergyFraction->GetYaxis()->SetTitleSize(0.06);
    fHistSignalEnergyFraction->GetXaxis()->SetTitleSize(0.06);
    fHistSignalEnergyFraction->Draw();
    fHistSignalEnergyFraction->Draw("same");
    fHistNeutronEnergyFraction->Draw("same");
    fHistSecondaryEnergyFraction->Draw("same");
    fHistHadronEnergyFraction->Draw("same");
    fHistKaonEnergyFraction->Draw("same");

    leg->Draw();

    TString outfile = "/tmp/EmEtDistribution"+detname+".png";
    c1->SaveAs(outfile.Data());
}
