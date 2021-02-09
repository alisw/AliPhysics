void DrawTrends(const char* file="ProductionQA.hist.root")
{
  // Draw trend plots using an output of the automatic QA procedure for PHOS
  // (see https://aliqapho.web.cern.ch/aliqapho).
  
  // Author: Boris Polishchuk

  TGaxis::SetMaxDigits(3);
  gStyle->SetPadRightMargin(0.003);
  gStyle->SetPadLeftMargin(0.03);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TCanvas* c1 = new TCanvas("c1","c1",1600,900);
  c1->Divide(1,3);
  
  TFile *_file0 = TFile::Open(file);
  gStyle->SetOptStat(0);

  c1->cd(1);
  TH1* havCluEnergySM1 = (TH1*)_file0->Get("havCluEnergySM1");
  havCluEnergySM1->GetYaxis()->SetRangeUser(0.31,0.99);
  havCluEnergySM1->SetTitle("Average cluster energy");
  havCluEnergySM1->SetLineColor(kRed);
  havCluEnergySM1->SetLineWidth(2);
  havCluEnergySM1->GetYaxis()->SetTitle("Energy, GeV");
  havCluEnergySM1->GetYaxis()->SetTitleSize(0.06);
  havCluEnergySM1->GetYaxis()->SetTitleOffset(0.25);
  havCluEnergySM1->GetYaxis()->SetLabelSize(0.06);
  havCluEnergySM1->SetLabelSize(0.05,"X");
  havCluEnergySM1->SetLabelOffset(0.002,"Y");
  havCluEnergySM1->Draw();

  TH1* havCluEnergySM2 = (TH1*)_file0->Get("havCluEnergySM2");
  havCluEnergySM2->SetLineColor(kGreen);
  havCluEnergySM2->SetLineWidth(2);
  havCluEnergySM2->Draw("same");

  TH1* havCluEnergySM3 = (TH1*)_file0->Get("havCluEnergySM3");
  havCluEnergySM3->SetLineColor(kBlue);
  havCluEnergySM3->SetLineWidth(2);
  havCluEnergySM3->Draw("same");

  TH1* havCluEnergySM4 = (TH1*)_file0->Get("havCluEnergySM4");
  havCluEnergySM4->SetLineColor(50);
  havCluEnergySM4->SetLineWidth(2);
  havCluEnergySM4->Draw("same");

  TLegend *myLegend = new TLegend(0.6,0.7,0.9,0.8);
  myLegendSetUp(myLegend,0.04);

  myLegend->AddEntry(havCluEnergySM1,"Module 4","l");
  myLegend->AddEntry(havCluEnergySM2,"Module 3","l");
  myLegend->AddEntry(havCluEnergySM3,"Module 2","l");
  myLegend->AddEntry(havCluEnergySM4,"Module 1","l");
  myLegend->Draw();

  TVirtualPad* pad2 = c1->cd(2);
  pad2->SetLogy();
  
  TH1* havCluMultSM1 = (TH1*)_file0->Get("havCluMultSM1");
  havCluMultSM1->GetYaxis()->SetRangeUser(0.01,1.99);
  havCluMultSM1->SetTitle("Average number of clusters per event");
  havCluMultSM1->SetLineColor(kRed);
  havCluMultSM1->SetLineWidth(2);
  havCluMultSM1->GetYaxis()->SetTitle("Number of clusters");
  havCluMultSM1->GetYaxis()->SetTitleSize(0.06);
  havCluMultSM1->GetYaxis()->SetTitleOffset(0.25);
  havCluMultSM1->GetYaxis()->SetLabelSize(0.06);
  havCluMultSM1->GetYaxis()->SetNoExponent(kTRUE);
  havCluMultSM1->SetLabelSize(0.05,"X");
  havCluMultSM1->SetLabelOffset(0.001,"Y");
  havCluMultSM1->Draw();

  TH1* havCluMultSM2 = (TH1*)_file0->Get("havCluMultSM2");
  havCluMultSM2->SetLineColor(kGreen);
  havCluMultSM2->SetLineWidth(2);
  havCluMultSM2->Draw("same");

  TH1* havCluMultSM3 = (TH1*)_file0->Get("havCluMultSM3");
  havCluMultSM3->SetLineColor(kBlue);
  havCluMultSM3->SetLineWidth(2);
  havCluMultSM3->Draw("same");

  TH1* havCluMultSM4 = (TH1*)_file0->Get("havCluMultSM4");
  havCluMultSM4->SetLineColor(50);
  havCluMultSM4->SetLineWidth(2);
  havCluMultSM4->Draw("same");

  myLegend->Draw();
  
  TVirtualPad* pad3 = c1->cd(3);
  
  TH1* havNcellPerCluSM1 = (TH1*)_file0->Get("havNcellPerCluSM1");
  havNcellPerCluSM1->GetYaxis()->SetRangeUser(0.,7.);
  havNcellPerCluSM1->SetTitle("Average number of cells in cluster");
  havNcellPerCluSM1->SetLineColor(kRed);
  havNcellPerCluSM1->SetLineWidth(2);
  havNcellPerCluSM1->GetYaxis()->SetTitle("Number of cells");
  havNcellPerCluSM1->GetYaxis()->SetTitleSize(0.06);
  havNcellPerCluSM1->GetYaxis()->SetTitleOffset(0.25);
  havNcellPerCluSM1->GetYaxis()->SetLabelSize(0.06);
  havNcellPerCluSM1->SetLabelSize(0.05,"X");
  havNcellPerCluSM1->SetLabelOffset(0.002,"Y");
  havNcellPerCluSM1->Draw();

  TH1* havNcellPerCluSM2 = (TH1*)_file0->Get("havNcellPerCluSM2");
  havNcellPerCluSM2->SetLineColor(kGreen);
  havNcellPerCluSM2->SetLineWidth(2);
  havNcellPerCluSM2->Draw("same");

  TH1* havNcellPerCluSM3 = (TH1*)_file0->Get("havNcellPerCluSM3");
  havNcellPerCluSM3->SetLineColor(kBlue);
  havNcellPerCluSM3->SetLineWidth(2);
  havNcellPerCluSM3->Draw("same");

  TH1* havNcellPerCluSM4 = (TH1*)_file0->Get("havNcellPerCluSM4");
  havNcellPerCluSM4->SetLineColor(50);
  havNcellPerCluSM4->SetLineWidth(2);
  havNcellPerCluSM4->Draw("same");

  myLegend->Draw();
  
}
//=============================================================================
void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}
