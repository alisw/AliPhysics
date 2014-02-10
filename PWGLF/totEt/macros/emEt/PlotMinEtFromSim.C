void SetStyles(TH1 *histo,int marker, int color){
  histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  //histo->GetXaxis()->SetTitle(xtitle);
  //histo->GetYaxis()->SetTitle(ytitle);
}
void PlotMinEtFromSim(Bool_t isPhos = kFALSE){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  Float_t min = 0;
  float max = 1;
  TString filename, detname;
  if(isPhos){
    min = 0.655;
    max = 0.785;
    detname = "PHOS";
    filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.Run139465.root";
  }
  else{
    min = 0.58;
    max = 0.725;
    filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root";
    detname = "EMCal";
  }
  
  TFile *f = TFile::Open(filename, "READ");
  TList *l = dynamic_cast<TList*>(f->Get("out1"));
  TH1F *fHistSimulatedGammaEnergyAboveThreshold = l->FindObject("fHistSimulatedGammaEnergyAboveThreshold");
  TH1F *fHistSimulatedGammaEnergy = l->FindObject("fHistSimulatedGammaEnergy");
  SetStyles(fHistSimulatedGammaEnergyAboveThreshold,20,TColor::kRed);
  fHistSimulatedGammaEnergyAboveThreshold->Divide(fHistSimulatedGammaEnergy);

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
    fHistSimulatedGammaEnergyAboveThreshold->SetMaximum(max +0.1);
    fHistSimulatedGammaEnergyAboveThreshold->SetMinimum(min-0.1);
    fHistSimulatedGammaEnergyAboveThreshold->GetXaxis()->SetTitle("Centrality bin");
    fHistSimulatedGammaEnergyAboveThreshold->GetYaxis()->SetTitle("f_{minEt}");
    fHistSimulatedGammaEnergyAboveThreshold->GetYaxis()->SetLabelSize(0.06);
    fHistSimulatedGammaEnergyAboveThreshold->GetXaxis()->SetLabelSize(0.06);
    fHistSimulatedGammaEnergyAboveThreshold->GetYaxis()->SetTitleSize(0.06);
    fHistSimulatedGammaEnergyAboveThreshold->GetXaxis()->SetTitleSize(0.06);
    fHistSimulatedGammaEnergyAboveThreshold->Draw();
    TLine *lineMin = new TLine(-0.5,min,19.5,min);
    lineMin->Draw();
    TLine *lineMax = new TLine(-0.5,max,19.5,max);
    lineMax->Draw();
    lineMin->SetLineColor(TColor::kBlue);
    lineMax->SetLineColor(TColor::kBlue);
    lineMin->SetLineStyle(2);
    lineMax->SetLineStyle(2);

    TString outfile = "/tmp/MinEtFromSim"+detname+".png";
    c1->SaveAs(outfile.Data());
}
