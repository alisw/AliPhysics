DrawTriggerQA(const char* filename="TriggerQA_253436-253591.root",
	      const char* trigger="L0")
{
  //Author: Yuri Kharlov (Yuri.Kharlov@cern.ch)
  
  file = new TFile(filename);
  hList = (TList*)file->Get(Form("PHOSTriggerQAResults%s",trigger));
  
  h4x4SM1 = (TH2F*)hList->FindObject("h4x4SM1");
  h4x4SM2 = (TH2F*)hList->FindObject("h4x4SM2");
  h4x4SM3 = (TH2F*)hList->FindObject("h4x4SM3");
  h4x4SM4 = (TH2F*)hList->FindObject("h4x4SM4");

  h4x4CluSM1 = (TH2F*)hList->FindObject("h4x4CluSM1");
  h4x4CluSM2 = (TH2F*)hList->FindObject("h4x4CluSM2");
  h4x4CluSM3 = (TH2F*)hList->FindObject("h4x4CluSM3");
  h4x4CluSM4 = (TH2F*)hList->FindObject("h4x4CluSM4");

  hPhotTrigSM1i = (TH1F*)hList->FindObject("hPhotTrigSM1");
  hPhotTrigSM2i = (TH1F*)hList->FindObject("hPhotTrigSM2");
  hPhotTrigSM3i = (TH1F*)hList->FindObject("hPhotTrigSM3");
  hPhotTrigSM4i = (TH1F*)hList->FindObject("hPhotTrigSM4");

  hPhotAllSM1i = (TH1F*)hList->FindObject("hPhotAllSM1");
  hPhotAllSM2i = (TH1F*)hList->FindObject("hPhotAllSM2");
  hPhotAllSM3i = (TH1F*)hList->FindObject("hPhotAllSM3");
  hPhotAllSM4i = (TH1F*)hList->FindObject("hPhotAllSM4");

  h4x4SM1->Rebin2D(2,2);
  h4x4SM2->Rebin2D(2,2);
  h4x4SM3->Rebin2D(2,2);
  h4x4SM4->Rebin2D(2,2);

  h4x4CluSM1->Rebin2D(2,2);  
  h4x4CluSM2->Rebin2D(2,2);  
  h4x4CluSM3->Rebin2D(2,2);  
  h4x4CluSM4->Rebin2D(2,2);  

  Double_t ptBins[] = {0.0,0.2,0.4,0.6,0.8, 1.0,1.2,1.4,1.6,1.8,
		       2.0,2.2,2.4,2.6,2.8, 3.0,3.2,3.4,3.6,3.8,
		       4.0,4.2,4.4,4.6,4.8, 5.0,5.5,6.0,6.5,7.0,
		       7.5,8.0,9.0,10.,12., 14.,20.,30.,40.};
  Int_t nPt = 38;

  TH1F *hPhotTrigSM1 = hPhotTrigSM1i->Rebin(nPt,"hPhotTrigSM1",ptBins);
  TH1F *hPhotTrigSM2 = hPhotTrigSM2i->Rebin(nPt,"hPhotTrigSM2",ptBins);
  TH1F *hPhotTrigSM3 = hPhotTrigSM3i->Rebin(nPt,"hPhotTrigSM3",ptBins);
  TH1F *hPhotTrigSM4 = hPhotTrigSM4i->Rebin(nPt,"hPhotTrigSM4",ptBins);

  hPhotTrigSM1->Sumw2();
  hPhotTrigSM2->Sumw2();
  hPhotTrigSM3->Sumw2();
  hPhotTrigSM4->Sumw2();

  for (Int_t iPt=1; iPt<=nPt; iPt++) {
    Double_t dPt = ptBins[iPt] - ptBins[iPt-1];
    hPhotTrigSM1->SetBinContent(iPt,hPhotTrigSM1->GetBinContent(iPt)/dPt);
    hPhotTrigSM2->SetBinContent(iPt,hPhotTrigSM2->GetBinContent(iPt)/dPt);
    hPhotTrigSM3->SetBinContent(iPt,hPhotTrigSM3->GetBinContent(iPt)/dPt);
    hPhotTrigSM4->SetBinContent(iPt,hPhotTrigSM4->GetBinContent(iPt)/dPt);
    // printf("pt bin %d is [%.1f - %.1f], dN/dpt=%g, dPt=%.1f\n",
    // 	   iPt,ptBins[iPt-1],ptBins[iPt],hPhotTrigSM1->GetBinContent(iPt),dPt);
  }

  SetMyStyle(hPhotTrigSM1,kRed);
  SetMyStyle(hPhotTrigSM2,kRed);
  SetMyStyle(hPhotTrigSM3,kRed);
  SetMyStyle(hPhotTrigSM4,kRed);
  
  TH1F *hPhotAllSM1 = hPhotAllSM1i->Rebin(nPt,"hPhotAllSM1",ptBins);
  TH1F *hPhotAllSM2 = hPhotAllSM2i->Rebin(nPt,"hPhotAllSM2",ptBins);
  TH1F *hPhotAllSM3 = hPhotAllSM3i->Rebin(nPt,"hPhotAllSM3",ptBins);
  TH1F *hPhotAllSM4 = hPhotAllSM4i->Rebin(nPt,"hPhotAllSM4",ptBins);

  hPhotAllSM1->Sumw2();
  hPhotAllSM2->Sumw2();
  hPhotAllSM3->Sumw2();
  hPhotAllSM4->Sumw2();

  for (Int_t iPt=1; iPt<nPt; iPt++) {
    Double_t dPt = ptBins[iPt] - ptBins[iPt-1];
    hPhotAllSM1->SetBinContent(iPt,hPhotAllSM1->GetBinContent(iPt)/dPt);
    hPhotAllSM2->SetBinContent(iPt,hPhotAllSM2->GetBinContent(iPt)/dPt);
    hPhotAllSM3->SetBinContent(iPt,hPhotAllSM3->GetBinContent(iPt)/dPt);
    hPhotAllSM4->SetBinContent(iPt,hPhotAllSM4->GetBinContent(iPt)/dPt);
  }
  
  SetMyStyle(hPhotAllSM1,kBlue);
  SetMyStyle(hPhotAllSM2,kBlue);
  SetMyStyle(hPhotAllSM3,kBlue);
  SetMyStyle(hPhotAllSM4,kBlue);
  
  rTrigRatioSM1 = (TH1F*)hPhotTrigSM1->Clone("rTrigRatioSM1");
  rTrigRatioSM2 = (TH1F*)hPhotTrigSM2->Clone("rTrigRatioSM2");
  rTrigRatioSM3 = (TH1F*)hPhotTrigSM3->Clone("rTrigRatioSM3");
  rTrigRatioSM4 = (TH1F*)hPhotTrigSM4->Clone("rTrigRatioSM4");

  rTrigRatioSM1->Divide(rTrigRatioSM1,hPhotAllSM1,1,1,"B");
  rTrigRatioSM2->Divide(rTrigRatioSM2,hPhotAllSM2,1,1,"B");
  rTrigRatioSM3->Divide(rTrigRatioSM3,hPhotAllSM3,1,1,"B");
  rTrigRatioSM4->Divide(rTrigRatioSM4,hPhotAllSM4,1,1,"B");

  SetMyStyle(rTrigRatioSM1,kBlack);
  SetMyStyle(rTrigRatioSM2,kBlack);
  SetMyStyle(rTrigRatioSM3,kBlack);
  SetMyStyle(rTrigRatioSM4,kBlack);

  rTrigRatioSM1->SetTitle("Triggered clusters / All clusters in SM1");
  rTrigRatioSM2->SetTitle("Triggered clusters / All clusters in SM2");
  rTrigRatioSM3->SetTitle("Triggered clusters / All clusters in SM3");
  rTrigRatioSM4->SetTitle("Triggered clusters / All clusters in SM4");
  //----------------------------------------------------------------------
  
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadLeftMargin(0.05);
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  
  c1 = new TCanvas("c1","c1",0,0,1000,600);
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetLogz();
  h4x4SM1->Draw("colz");
  c1->cd(2);
  gPad->SetLogz();
  h4x4SM2->Draw("colz");
  c1->cd(3);
  gPad->SetLogz();
  h4x4SM3->Draw("colz");
  c1->cd(4);
  gPad->SetLogz();
  h4x4SM4->Draw("colz");
  c1->Print(Form("TrigOccupancy_%s.eps",trigger));
  
  c2 = new TCanvas("c2","c2",0,0,1000,600);
  c2->Divide(2,2);
  c2->cd(1);
  gPad->SetLogz();
  h4x4CluSM1->Draw("colz");
  c2->cd(2);
  gPad->SetLogz();
  h4x4CluSM2->Draw("colz");
  c2->cd(3);
  gPad->SetLogz();
  h4x4CluSM3->Draw("colz");
  c2->cd(4);
  gPad->SetLogz();
  h4x4CluSM4->Draw("colz");
  c2->Print(Form("TrigOccupancyCluster_%s.eps",trigger));

  c3 = new TCanvas("c3","c3",0,0,1000,600);
  c3->Divide(2,2);
  c3->cd(1);
  gPad->SetLogy();
  hPhotTrigSM1->Draw();
  c3->cd(2);
  gPad->SetLogy();
  hPhotTrigSM2->Draw();
  c3->cd(3);
  gPad->SetLogy();
  hPhotTrigSM3->Draw();
  c3->cd(4);
  gPad->SetLogy();
  hPhotTrigSM4->Draw();
  c3->Print(Form("TrigPhotE_triggered_%s.eps",trigger));

  TLegend *leg = new TLegend(0.6,0.78,0.89,0.9);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->AddEntry(hPhotAllSM1,"All clusters","l");
  leg->AddEntry(hPhotTrigSM1,"Triggered clusters","l");
		
  c4 = new TCanvas("c4","c4",0,0,1000,600);
  c4->Divide(2,2);
  c4->cd(1);
  gPad->SetLogy();
  hPhotAllSM1->Draw();
  hPhotTrigSM1->Draw("same");
  leg->Draw();
  c4->cd(2);
  gPad->SetLogy();
  hPhotAllSM2->Draw();
  hPhotTrigSM2->Draw("same");
  c4->cd(3);
  gPad->SetLogy();
  hPhotAllSM3->Draw();
  hPhotTrigSM3->Draw("same");
  c4->cd(4);
  gPad->SetLogy();
  hPhotAllSM4->Draw();
  hPhotTrigSM4->Draw("same");
  c4->Print(Form("TrigPhotE_all_%s.eps",trigger));

  c5 = new TCanvas("c5","c5",0,0,1000,600);
  c5->Divide(2,2);
  c5->cd(1);
  gPad->SetLogy();
  rTrigRatioSM1->Draw();
  c5->cd(2);
  gPad->SetLogy();
  rTrigRatioSM2->Draw();
  c5->cd(3);
  gPad->SetLogy();
  rTrigRatioSM3->Draw();
  c5->cd(4);
  gPad->SetLogy();
  rTrigRatioSM4->Draw();
  c5->Print(Form("TrigPhotE_ratio_%s.eps",trigger));
}
//-----------------------------------------------------------------------------
DrawTriggerQATRU(const char* filename="TriggerQA_253436-253591.root",
		 const char* trigger ="L0",
		 const char* module  ="M2")
{
  file = new TFile(filename);
  hList = (TList*)file->Get(Form("PHOSTriggerQAResults%s",trigger));
  
  hPhotAllTRU1i = (TH1F*)hList->FindObject(Form("hPhotAllS%sTRU1",module));
  hPhotAllTRU2i = (TH1F*)hList->FindObject(Form("hPhotAllS%sTRU2",module));
  hPhotAllTRU3i = (TH1F*)hList->FindObject(Form("hPhotAllS%sTRU3",module));
  hPhotAllTRU4i = (TH1F*)hList->FindObject(Form("hPhotAllS%sTRU4",module));
  hPhotAllTRU5i = (TH1F*)hList->FindObject(Form("hPhotAllS%sTRU5",module));
  hPhotAllTRU6i = (TH1F*)hList->FindObject(Form("hPhotAllS%sTRU6",module));
  hPhotAllTRU7i = (TH1F*)hList->FindObject(Form("hPhotAllS%sTRU7",module));
  hPhotAllTRU8i = (TH1F*)hList->FindObject(Form("hPhotAllS%sTRU8",module));

  hPhotTrigTRU1i = (TH1F*)hList->FindObject(Form("hPhotTrigS%sTRU1",module));
  hPhotTrigTRU2i = (TH1F*)hList->FindObject(Form("hPhotTrigS%sTRU2",module));
  hPhotTrigTRU3i = (TH1F*)hList->FindObject(Form("hPhotTrigS%sTRU3",module));
  hPhotTrigTRU4i = (TH1F*)hList->FindObject(Form("hPhotTrigS%sTRU4",module));
  hPhotTrigTRU5i = (TH1F*)hList->FindObject(Form("hPhotTrigS%sTRU5",module));
  hPhotTrigTRU6i = (TH1F*)hList->FindObject(Form("hPhotTrigS%sTRU6",module));
  hPhotTrigTRU7i = (TH1F*)hList->FindObject(Form("hPhotTrigS%sTRU7",module));
  hPhotTrigTRU8i = (TH1F*)hList->FindObject(Form("hPhotTrigS%sTRU8",module));

  Double_t ptBins[] = {0.0,0.2,0.4,0.6,0.8, 1.0,1.2,1.4,1.6,1.8,
		       2.0,2.2,2.4,2.6,2.8, 3.0,3.2,3.4,3.6,3.8,
		       4.0,4.2,4.4,4.6,4.8, 5.0,5.5,6.0,6.5,7.0,
		       7.5,8.0,9.0,10.,12., 14.,20.,30.,40.};
  Int_t nPt = 38;

  TH1F *hPhotTrigTRU1 = hPhotTrigTRU1i->Rebin(nPt,"hPhotTrigTRU1",ptBins);
  TH1F *hPhotTrigTRU2 = hPhotTrigTRU2i->Rebin(nPt,"hPhotTrigTRU2",ptBins);
  TH1F *hPhotTrigTRU3 = hPhotTrigTRU3i->Rebin(nPt,"hPhotTrigTRU3",ptBins);
  TH1F *hPhotTrigTRU4 = hPhotTrigTRU4i->Rebin(nPt,"hPhotTrigTRU4",ptBins);
  TH1F *hPhotTrigTRU5 = hPhotTrigTRU5i->Rebin(nPt,"hPhotTrigTRU5",ptBins);
  TH1F *hPhotTrigTRU6 = hPhotTrigTRU6i->Rebin(nPt,"hPhotTrigTRU6",ptBins);
  TH1F *hPhotTrigTRU7 = hPhotTrigTRU7i->Rebin(nPt,"hPhotTrigTRU7",ptBins);
  TH1F *hPhotTrigTRU8 = hPhotTrigTRU8i->Rebin(nPt,"hPhotTrigTRU8",ptBins);

  TH1F *hPhotAllTRU1 = hPhotAllTRU1i->Rebin(nPt,"hPhotAllTRU1",ptBins);
  TH1F *hPhotAllTRU2 = hPhotAllTRU2i->Rebin(nPt,"hPhotAllTRU2",ptBins);
  TH1F *hPhotAllTRU3 = hPhotAllTRU3i->Rebin(nPt,"hPhotAllTRU3",ptBins);
  TH1F *hPhotAllTRU4 = hPhotAllTRU4i->Rebin(nPt,"hPhotAllTRU4",ptBins);
  TH1F *hPhotAllTRU5 = hPhotAllTRU5i->Rebin(nPt,"hPhotAllTRU5",ptBins);
  TH1F *hPhotAllTRU6 = hPhotAllTRU6i->Rebin(nPt,"hPhotAllTRU6",ptBins);
  TH1F *hPhotAllTRU7 = hPhotAllTRU7i->Rebin(nPt,"hPhotAllTRU7",ptBins);
  TH1F *hPhotAllTRU8 = hPhotAllTRU8i->Rebin(nPt,"hPhotAllTRU8",ptBins);

  hPhotAllTRU1->Sumw2();
  hPhotAllTRU2->Sumw2();
  hPhotAllTRU3->Sumw2();
  hPhotAllTRU4->Sumw2();
  hPhotAllTRU5->Sumw2();
  hPhotAllTRU6->Sumw2();
  hPhotAllTRU7->Sumw2();
  hPhotAllTRU8->Sumw2();

  hPhotTrigTRU1->Sumw2();
  hPhotTrigTRU2->Sumw2();
  hPhotTrigTRU3->Sumw2();
  hPhotTrigTRU4->Sumw2();
  hPhotTrigTRU5->Sumw2();
  hPhotTrigTRU6->Sumw2();
  hPhotTrigTRU7->Sumw2();
  hPhotTrigTRU8->Sumw2();

  for (Int_t iPt=1; iPt<=nPt; iPt++) {
    Double_t dPt = ptBins[iPt] - ptBins[iPt-1];
    hPhotTrigTRU1->SetBinContent(iPt,hPhotTrigTRU1->GetBinContent(iPt)/dPt);
    hPhotTrigTRU2->SetBinContent(iPt,hPhotTrigTRU2->GetBinContent(iPt)/dPt);
    hPhotTrigTRU3->SetBinContent(iPt,hPhotTrigTRU3->GetBinContent(iPt)/dPt);
    hPhotTrigTRU4->SetBinContent(iPt,hPhotTrigTRU4->GetBinContent(iPt)/dPt);
    hPhotTrigTRU5->SetBinContent(iPt,hPhotTrigTRU5->GetBinContent(iPt)/dPt);
    hPhotTrigTRU6->SetBinContent(iPt,hPhotTrigTRU6->GetBinContent(iPt)/dPt);
    hPhotTrigTRU7->SetBinContent(iPt,hPhotTrigTRU7->GetBinContent(iPt)/dPt);
    hPhotTrigTRU8->SetBinContent(iPt,hPhotTrigTRU8->GetBinContent(iPt)/dPt);

    hPhotAllTRU1->SetBinContent(iPt,hPhotAllTRU1->GetBinContent(iPt)/dPt);
    hPhotAllTRU2->SetBinContent(iPt,hPhotAllTRU2->GetBinContent(iPt)/dPt);
    hPhotAllTRU3->SetBinContent(iPt,hPhotAllTRU3->GetBinContent(iPt)/dPt);
    hPhotAllTRU4->SetBinContent(iPt,hPhotAllTRU4->GetBinContent(iPt)/dPt);
    hPhotAllTRU5->SetBinContent(iPt,hPhotAllTRU5->GetBinContent(iPt)/dPt);
    hPhotAllTRU6->SetBinContent(iPt,hPhotAllTRU6->GetBinContent(iPt)/dPt);
    hPhotAllTRU7->SetBinContent(iPt,hPhotAllTRU7->GetBinContent(iPt)/dPt);
    hPhotAllTRU8->SetBinContent(iPt,hPhotAllTRU8->GetBinContent(iPt)/dPt);
  }

  SetMyStyle(hPhotTrigTRU1,kRed);
  SetMyStyle(hPhotTrigTRU2,kRed);
  SetMyStyle(hPhotTrigTRU3,kRed);
  SetMyStyle(hPhotTrigTRU4,kRed);
  SetMyStyle(hPhotTrigTRU5,kRed);
  SetMyStyle(hPhotTrigTRU6,kRed);
  SetMyStyle(hPhotTrigTRU7,kRed);
  SetMyStyle(hPhotTrigTRU8,kRed);

  SetMyStyle(hPhotAllTRU1,kBlue);
  SetMyStyle(hPhotAllTRU2,kBlue);
  SetMyStyle(hPhotAllTRU3,kBlue);
  SetMyStyle(hPhotAllTRU4,kBlue);
  SetMyStyle(hPhotAllTRU5,kBlue);
  SetMyStyle(hPhotAllTRU6,kBlue);
  SetMyStyle(hPhotAllTRU7,kBlue);
  SetMyStyle(hPhotAllTRU8,kBlue);

  TLegend *leg = new TLegend(0.6,0.78,0.89,0.9);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->AddEntry(hPhotAllTRU1,"All clusters","l");
  leg->AddEntry(hPhotTrigTRU1,"Triggered clusters","l");
  
  TH1F *rTrigRatioTRU1 = (TH1F*)hPhotTrigTRU1->Clone("rTrigRatioTRU1");
  TH1F *rTrigRatioTRU2 = (TH1F*)hPhotTrigTRU2->Clone("rTrigRatioTRU1");
  TH1F *rTrigRatioTRU3 = (TH1F*)hPhotTrigTRU3->Clone("rTrigRatioTRU1");
  TH1F *rTrigRatioTRU4 = (TH1F*)hPhotTrigTRU4->Clone("rTrigRatioTRU1");
  TH1F *rTrigRatioTRU5 = (TH1F*)hPhotTrigTRU5->Clone("rTrigRatioTRU1");
  TH1F *rTrigRatioTRU6 = (TH1F*)hPhotTrigTRU6->Clone("rTrigRatioTRU1");
  TH1F *rTrigRatioTRU7 = (TH1F*)hPhotTrigTRU7->Clone("rTrigRatioTRU1");
  TH1F *rTrigRatioTRU8 = (TH1F*)hPhotTrigTRU8->Clone("rTrigRatioTRU1");

  rTrigRatioTRU1->Divide(rTrigRatioTRU1,hPhotAllTRU1,1,1,"B");
  rTrigRatioTRU2->Divide(rTrigRatioTRU2,hPhotAllTRU2,1,1,"B");
  rTrigRatioTRU3->Divide(rTrigRatioTRU3,hPhotAllTRU3,1,1,"B");
  rTrigRatioTRU4->Divide(rTrigRatioTRU4,hPhotAllTRU4,1,1,"B");
  rTrigRatioTRU5->Divide(rTrigRatioTRU5,hPhotAllTRU5,1,1,"B");
  rTrigRatioTRU6->Divide(rTrigRatioTRU6,hPhotAllTRU6,1,1,"B");
  rTrigRatioTRU7->Divide(rTrigRatioTRU7,hPhotAllTRU7,1,1,"B");
  rTrigRatioTRU8->Divide(rTrigRatioTRU8,hPhotAllTRU8,1,1,"B");

  SetMyStyle(rTrigRatioTRU1,kBlack,Form("Trigger efficiency in TRU1, %s",module),"Triggered clusters / All clusters");
  SetMyStyle(rTrigRatioTRU2,kBlack,Form("Trigger efficiency in TRU2, %s",module),"Triggered clusters / All clusters");
  SetMyStyle(rTrigRatioTRU3,kBlack,Form("Trigger efficiency in TRU3, %s",module),"Triggered clusters / All clusters");
  SetMyStyle(rTrigRatioTRU4,kBlack,Form("Trigger efficiency in TRU4, %s",module),"Triggered clusters / All clusters");
  SetMyStyle(rTrigRatioTRU5,kBlack,Form("Trigger efficiency in TRU5, %s",module),"Triggered clusters / All clusters");
  SetMyStyle(rTrigRatioTRU6,kBlack,Form("Trigger efficiency in TRU6, %s",module),"Triggered clusters / All clusters");
  SetMyStyle(rTrigRatioTRU7,kBlack,Form("Trigger efficiency in TRU7, %s",module),"Triggered clusters / All clusters");
  SetMyStyle(rTrigRatioTRU8,kBlack,Form("Trigger efficiency in TRU8, %s",module),"Triggered clusters / All clusters");
    
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadLeftMargin(0.09);
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  
  c21 = new TCanvas("c21","c21",0,0,1000,600);
  c21->Divide(4,2);
  c21->cd(1);
  gPad->SetLogy();
  hPhotAllTRU1->Draw();
  hPhotTrigTRU1->Draw("same");
  leg->Draw();
  c21->cd(5);
  gPad->SetLogy();
  hPhotAllTRU2->Draw();
  hPhotTrigTRU2->Draw("same");
  c21->cd(2);
  gPad->SetLogy();
  hPhotAllTRU3->Draw();
  hPhotTrigTRU3->Draw("same");
  c21->cd(6);
  gPad->SetLogy();
  hPhotAllTRU4->Draw();
  hPhotTrigTRU4->Draw("same");
  c21->cd(3);
  gPad->SetLogy();
  hPhotAllTRU5->Draw();
  hPhotTrigTRU5->Draw("same");
  c21->cd(7);
  gPad->SetLogy();
  hPhotAllTRU6->Draw();
  hPhotTrigTRU6->Draw("same");
  c21->cd(4);
  gPad->SetLogy();
  hPhotAllTRU7->Draw();
  hPhotTrigTRU7->Draw("same");
  c21->cd(8);
  gPad->SetLogy();
  hPhotAllTRU8->Draw();
  hPhotTrigTRU8->Draw("same");

  c21->Print(Form("TrigPhotEperTRU_%s_%s.eps",trigger,module));

  c22 = new TCanvas("c22","c22",0,0,1000,600);
  c22->Divide(4,2);
  c22->cd(1);
  rTrigRatioTRU1->Draw();
  c22->cd(5);
  rTrigRatioTRU2->Draw();
  c22->cd(2);
  rTrigRatioTRU3->Draw();
  c22->cd(6);
  rTrigRatioTRU4->Draw();
  c22->cd(3);
  rTrigRatioTRU5->Draw();
  c22->cd(7);
  rTrigRatioTRU6->Draw();
  c22->cd(4);
  rTrigRatioTRU7->Draw();
  c22->cd(8);
  rTrigRatioTRU8->Draw();
  
  c22->Print(Form("TrigEffiperTRU_%s_%s.eps",trigger,module));

}
//-----------------------------------------------------------------------------
DrawTriggerEff(const char* filename="TriggerQA_253436-253591.root")
{
  
  file = new TFile(filename);
  hListL0 = (TList*)file->Get("PHOSTriggerQAResultsL0");
  hListL1 = (TList*)file->Get("PHOSTriggerQAResultsL1Low");

  hPhotTrigL0SM1i = (TH1F*)hListL0->FindObject("hPhotTrigSM1");
  hPhotTrigL0SM2i = (TH1F*)hListL0->FindObject("hPhotTrigSM2");
  hPhotTrigL0SM3i = (TH1F*)hListL0->FindObject("hPhotTrigSM3");
  hPhotTrigL0SM4i = (TH1F*)hListL0->FindObject("hPhotTrigSM4");

  hPhotTrigL1SM1i = (TH1F*)hListL1->FindObject("hPhotTrigSM1");
  hPhotTrigL1SM2i = (TH1F*)hListL1->FindObject("hPhotTrigSM2");
  hPhotTrigL1SM3i = (TH1F*)hListL1->FindObject("hPhotTrigSM3");
  hPhotTrigL1SM4i = (TH1F*)hListL1->FindObject("hPhotTrigSM4");

  Double_t ptBins[] = {0.0,0.2,0.4,0.6,0.8, 1.0,1.2,1.4,1.6,1.8,
		       2.0,2.2,2.4,2.6,2.8, 3.0,3.2,3.4,3.6,3.8,
		       4.0,4.2,4.4,4.6,4.8, 5.0,5.5,6.0,6.5,7.0,
		       7.5,8.0,9.0,10.,12., 14.,20.,30.,40.};
  Int_t nPt = 38;

  TH1F *hPhotTrigL0SM1 = hPhotTrigL0SM1i->Rebin(nPt,"hPhotTrigL0SM1",ptBins);
  TH1F *hPhotTrigL0SM2 = hPhotTrigL0SM2i->Rebin(nPt,"hPhotTrigL0SM2",ptBins);
  TH1F *hPhotTrigL0SM3 = hPhotTrigL0SM3i->Rebin(nPt,"hPhotTrigL0SM3",ptBins);
  TH1F *hPhotTrigL0SM4 = hPhotTrigL0SM4i->Rebin(nPt,"hPhotTrigL0SM4",ptBins);

  hPhotTrigL0SM1->Sumw2();
  hPhotTrigL0SM2->Sumw2();
  hPhotTrigL0SM3->Sumw2();
  hPhotTrigL0SM4->Sumw2();

  for (Int_t iPt=1; iPt<=nPt; iPt++) {
    Double_t dPt = ptBins[iPt] - ptBins[iPt-1];
    hPhotTrigL0SM1->SetBinContent(iPt,hPhotTrigL0SM1->GetBinContent(iPt)/dPt);
    hPhotTrigL0SM2->SetBinContent(iPt,hPhotTrigL0SM2->GetBinContent(iPt)/dPt);
    hPhotTrigL0SM3->SetBinContent(iPt,hPhotTrigL0SM3->GetBinContent(iPt)/dPt);
    hPhotTrigL0SM4->SetBinContent(iPt,hPhotTrigL0SM4->GetBinContent(iPt)/dPt);
  }
  
  TH1F *hPhotTrigL1SM1 = hPhotTrigL1SM1i->Rebin(nPt,"hPhotTrigL1SM1",ptBins);
  TH1F *hPhotTrigL1SM2 = hPhotTrigL1SM2i->Rebin(nPt,"hPhotTrigL1SM2",ptBins);
  TH1F *hPhotTrigL1SM3 = hPhotTrigL1SM3i->Rebin(nPt,"hPhotTrigL1SM3",ptBins);
  TH1F *hPhotTrigL1SM4 = hPhotTrigL1SM4i->Rebin(nPt,"hPhotTrigL1SM4",ptBins);

  hPhotTrigL1SM1->Sumw2();
  hPhotTrigL1SM2->Sumw2();
  hPhotTrigL1SM3->Sumw2();
  hPhotTrigL1SM4->Sumw2();

  for (Int_t iPt=1; iPt<=nPt; iPt++) {
    Double_t dPt = ptBins[iPt] - ptBins[iPt-1];
    hPhotTrigL1SM1->SetBinContent(iPt,hPhotTrigL1SM1->GetBinContent(iPt)/dPt);
    hPhotTrigL1SM2->SetBinContent(iPt,hPhotTrigL1SM2->GetBinContent(iPt)/dPt);
    hPhotTrigL1SM3->SetBinContent(iPt,hPhotTrigL1SM3->GetBinContent(iPt)/dPt);
    hPhotTrigL1SM4->SetBinContent(iPt,hPhotTrigL1SM4->GetBinContent(iPt)/dPt);
  }

  hRatioL1L0SM1 = (TH1F*)hPhotTrigL1SM1->Clone("hRatioL1L0SM1");
  hRatioL1L0SM2 = (TH1F*)hPhotTrigL1SM2->Clone("hRatioL1L0SM2");
  hRatioL1L0SM3 = (TH1F*)hPhotTrigL1SM3->Clone("hRatioL1L0SM3");
  hRatioL1L0SM4 = (TH1F*)hPhotTrigL1SM4->Clone("hRatioL1L0SM4");

  hRatioL1L0SM1->SetAxisRange(0.,19.);
  hRatioL1L0SM2->SetAxisRange(0.,19.);
  hRatioL1L0SM3->SetAxisRange(0.,19.);
  hRatioL1L0SM4->SetAxisRange(0.,19.);

  hRatioL1L0SM1->Divide(hRatioL1L0SM1,hPhotTrigL0SM1,1,1,"B");
  hRatioL1L0SM2->Divide(hRatioL1L0SM2,hPhotTrigL0SM2,1,1,"B");
  hRatioL1L0SM3->Divide(hRatioL1L0SM3,hPhotTrigL0SM3,1,1,"B");
  hRatioL1L0SM4->Divide(hRatioL1L0SM4,hPhotTrigL0SM4,1,1,"B");

  SetMyStyle(hRatioL1L0SM1,kBlack);
  SetMyStyle(hRatioL1L0SM2,kBlack);
  SetMyStyle(hRatioL1L0SM3,kBlack);
  SetMyStyle(hRatioL1L0SM4,kBlack);

  hRatioL1L0SM1->SetYTitle("L1/L0");
  hRatioL1L0SM2->SetYTitle("L1/L0");
  hRatioL1L0SM3->SetYTitle("L1/L0");
  hRatioL1L0SM4->SetYTitle("L1/L0");

  hRatioL1L0SM1->SetTitle("Triggered spectra ratio L1L / L0 in M1");
  hRatioL1L0SM2->SetTitle("Triggered spectra ratio L1L / L0 in M2");
  hRatioL1L0SM3->SetTitle("Triggered spectra ratio L1L / L0 in M3");
  hRatioL1L0SM4->SetTitle("Triggered spectra ratio L1L / L0 in M4");
  
  hRatioL1L0SM1->SetAxisRange(0.025,3.9,"Y");
  hRatioL1L0SM2->SetAxisRange(0.025,3.9,"Y");
  hRatioL1L0SM3->SetAxisRange(0.025,3.9,"Y");
  hRatioL1L0SM4->SetAxisRange(0.025,3.9,"Y");

  hRatioL1L0SM1->GetYaxis()->SetMoreLogLabels();
  hRatioL1L0SM2->GetYaxis()->SetMoreLogLabels();
  hRatioL1L0SM3->GetYaxis()->SetMoreLogLabels();
  hRatioL1L0SM4->GetYaxis()->SetMoreLogLabels();

  hRatioL1L0SM1->GetYaxis()->SetNoExponent();
  hRatioL1L0SM2->GetYaxis()->SetNoExponent();
  hRatioL1L0SM3->GetYaxis()->SetNoExponent();
  hRatioL1L0SM4->GetYaxis()->SetNoExponent();
  
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadLeftMargin(0.06);
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  
  c14 = new TCanvas("c4","c4",0,0,1000,600);
  c14->Divide(2,2);
  c14->cd(1);
  hRatioL1L0SM1->Draw();
  c14->cd(2);
  hRatioL1L0SM2->Draw();
  c14->cd(3);
  hRatioL1L0SM3->Draw();
  c14->cd(4);
  hRatioL1L0SM4->Draw();
  c14->Print("TrigEfficiency.eps");
}  
//-----------------------------------------------------------------------------
void SetMyStyle(TH1F * histo, const Int_t color = kRed)
{
  histo->SetLineWidth(2);
  histo->SetLineColor(color);
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(0.5);
  histo->SetMarkerColor(color);
  histo->SetXTitle("E, GeV");
  histo->SetNdivisions(530);
  histo->SetAxisRange(0.,39.);
}
//-----------------------------------------------------------------------------
void SetMyStyle(TH1F * histo, const Int_t color = kRed, char *title, const char *ytitle)
{
  histo->SetLineWidth(2);
  histo->SetLineColor(color);
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(0.5);
  histo->SetMarkerColor(color);
  histo->SetXTitle("E, GeV");
  histo->SetNdivisions(520);
  histo->SetTitle(title);
  histo->SetYTitle(ytitle);
  histo->SetTitleOffset(1.3,"Y");
  histo->SetAxisRange(0.,19.,"X");
  histo->SetAxisRange(0.,1.29,"Y");
}
