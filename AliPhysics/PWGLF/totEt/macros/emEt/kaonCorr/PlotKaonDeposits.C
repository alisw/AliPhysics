void ScaleByBinWidth(TH1 *histo);
void ScaleByBinWidth(TH2 *histo);
void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle, Bool_t scale = kFALSE);
void PlotKaonDeposits(TString filename="rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run137366.root"){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *f = TFile::Open(filename, "READ");
  TList *l = (TList*)f->Get("out1");
  TH1F *fHistSimKaonsInAcceptance = l->FindObject("fHistSimKaonsInAcceptance");
  TH1F *fHistSimK0SInAcceptance = l->FindObject("fHistSimK0SInAcceptance");
  TH1F *fHistSimK0LInAcceptance = l->FindObject("fHistSimK0LInAcceptance");
  TH1F *fHistSimKPlusInAcceptance = l->FindObject("fHistSimKPlusInAcceptance");
  TH1F *fHistSimKMinusInAcceptance = l->FindObject("fHistSimKMinusInAcceptance");
  TH1F *fHistSimKaonsOutOfAcceptance = l->FindObject("fHistSimKaonsOutOfAcceptance");
  TH1F *fHistSimKaonsInAcceptanceWithDepositsPrimaries = l->FindObject("fHistSimKaonsInAcceptanceWithDepositsPrimaries");
  TH1F *fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries = l->FindObject("fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries");
  TH1F *fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries = l->FindObject("fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries");
  
  SetStyles(fHistSimKaonsInAcceptance,20,2,"p_{T}","Number of kaons",kTRUE);
  SetStyles(fHistSimK0SInAcceptance,20,3,"p_{T}","Number of kaons",kTRUE);
  SetStyles(fHistSimK0LInAcceptance,20,4,"p_{T}","Number of kaons",kTRUE);
  SetStyles(fHistSimKMinusInAcceptance,20,5,"p_{T}","Number of kaons",kTRUE);
  SetStyles(fHistSimKPlusInAcceptance,20,6,"p_{T}","Number of kaons",kTRUE);
  SetStyles(fHistSimKaonsOutOfAcceptance,20,7,"p_{T}","Number of kaons",kTRUE);
  SetStyles(fHistSimKaonsInAcceptanceWithDepositsPrimaries,20,2,"p_{T}","Number of kaons",kTRUE);
  SetStyles(fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries,21,4,"p_{T}","Number of kaons",kTRUE);
  SetStyles(fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries,22,3,"p_{T}","Number of kaons",kTRUE);


  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLogy();
  fHistSimKaonsInAcceptance->Draw();
  fHistSimKaonsInAcceptance->GetXaxis()->SetRange(1,fHistSimKaonsInAcceptance->FindBin(3.0));
  fHistSimK0SInAcceptance->Draw("same");
  fHistSimK0LInAcceptance->Draw("same");
  fHistSimKMinusInAcceptance->Draw("same");
  fHistSimKPlusInAcceptance->Draw("same");

  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  c2->SetTopMargin(0.02);
  c2->SetRightMargin(0.02);
  c2->SetBorderSize(0);
  c2->SetFillColor(0);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderMode(0);
  c2->SetLogy();
  fHistSimKaonsInAcceptanceWithDepositsPrimaries->Draw();
  fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries->Draw("same");
  fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries->Draw("same");
  TLegend *leg1 = new TLegend(0.223154,0.685484,0.387584,0.965054);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(fHistSimKaonsInAcceptanceWithDepositsPrimaries,"Primary kaons with |y|<0.5 which deposited energy");
  leg1->AddEntry(fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries,"Primary kaons with |y|>0.5 which deposited energy");
  leg1->AddEntry(fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries,"Secondary kaons with |y|<0.5 which deposited energy");
  leg1->SetTextSize(0.0456989);
  leg1->Draw();


  TCanvas *c3 = new TCanvas("c3","c3",600,400);
  c3->SetTopMargin(0.02);
  c3->SetRightMargin(0.02);
  c3->SetBorderSize(0);
  c3->SetFillColor(0);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetFrameFillColor(0);
  c3->SetFrameBorderMode(0);
  TH1F *hTotal = fHistSimKaonsInAcceptanceWithDepositsPrimaries->Clone("hTotal");
  hTotal->Add(fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries);
  hTotal->Add(fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries);
  TH1F *hKaonsInAccPrimaries = fHistSimKaonsInAcceptanceWithDepositsPrimaries->Clone("hKaonsInAccPrimaries");
  TH1F *hKaonsOutAccPrimaries = fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries->Clone("hKaonsOutAccPrimaries");
  TH1F *hKaonsOutAccSecondary = fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries->Clone("hKaonsOutAccSecondaries");
  hKaonsInAccPrimaries->Divide(hTotal);
  hKaonsOutAccPrimaries->Divide(hTotal);
  hKaonsOutAccSecondaries->Divide(hTotal);
  hKaonsInAccPrimaries->SetMaximum(1.0);
  hKaonsInAccPrimaries->SetMinimum(0.0);
  hKaonsInAccPrimaries->GetXaxis()->SetRange(1,fHistSimKaonsInAcceptance->FindBin(3.0));
  hKaonsInAccPrimaries->Draw();
  hKaonsOutAccPrimaries->Draw("same");
  hKaonsOutAccSecondaries->Draw("same");
  TLegend *leg2 = new TLegend(0.223154,0.685484,0.387584,0.965054);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(hKaonsInAccPrimaries,"Primary kaons with |y|<0.5 which deposited energy");
  leg2->AddEntry(hKaonsOutAccPrimaries,"Primary kaons with |y|>0.5 which deposited energy");
  leg2->AddEntry(hKaonsOutAccSecondaries,"Secondary kaons with |y|<0.5 which deposited energy");
  leg2->SetTextSize(0.0456989);
  leg2->Draw();

  TCanvas *c4 = new TCanvas("c4","c4",600,400);
  c4->SetTopMargin(0.02);
  c4->SetRightMargin(0.02);
  c4->SetBorderSize(0);
  c4->SetFillColor(0);
  c4->SetFillColor(0);
  c4->SetBorderMode(0);
  c4->SetFrameFillColor(0);
  c4->SetFrameBorderMode(0);
  TH1F *hInAcc = fHistSimKaonsInAcceptance->Clone("hInAcc");
  TH1F *hOutAcc = fHistSimKaonsOutOfAcceptance->Clone("hOutAcc");
  TH1F *hTotal2 = fHistSimKaonsInAcceptance->Clone("hTotal2");
  hTotal2->Add(hOutAcc);
  hInAcc->Divide(hTotal2);
  hOutAcc->Divide(hTotal2);
  hInAcc->SetMaximum(0.13);
  hInAcc->SetMinimum(0.0);
  hOutAcc->Scale(0.1);
  hInAcc->GetXaxis()->SetRange(1,hInAcc->FindBin(3.0));
  hInAcc->Draw();
  hOutAcc->Draw("same");
  TLegend *leg3 = new TLegend(0.223154,0.685484,0.387584,0.965054);
  leg3->SetFillStyle(0);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.03);
  leg3->AddEntry(hInAcc,"Kaons with |y|<0.5");
  leg3->AddEntry(hOutAcc,"0.1*Kaons with |y|>0.5");
  leg3->SetTextSize(0.0456989);
  leg3->Draw();

  TH3F *fHistK0EDepositsVsPtInAcceptance = l->FindObject("fHistK0EDepositsVsPtInAcceptance");
  TH3F *fHistK0EGammaVsPtInAcceptance = l->FindObject("fHistK0EGammaVsPtInAcceptance");
  TH3F *fHistK0EDepositsVsPtOutOfAcceptance = l->FindObject("fHistK0EDepositsVsPtOutOfAcceptance");
  TH3F *fHistK0EGammaVsPtOutOfAcceptance = l->FindObject("fHistK0EGammaVsPtOutOfAcceptance");
  TObjArray hKaonDepositsInAcceptance(11);
  TObjArray hKaonGammaInAcceptance(11);
  TObjArray hKaonDepositsOutOfAcceptance(11);
  TObjArray hKaonGammaOutOfAcceptance(11);
  TObjArray hKaonDepositsInAcceptanceAvg(11);
  TObjArray hKaonGammaInAcceptanceAvg(11);
  TObjArray hKaonDepositsOutOfAcceptanceAvg(11);
  TObjArray hKaonGammaOutOfAcceptanceAvg(11);
  TObjArray hKaonDepositsInAcceptanceLimAvg(11);
  TObjArray hKaonGammaInAcceptanceLimAvg(11);
  TObjArray hKaonDepositsOutOfAcceptanceLimAvg(11);
  TObjArray hKaonGammaOutOfAcceptanceLimAvg(11);
  for(int bin=1;bin<=11;bin++){
    fHistK0EDepositsVsPtInAcceptance->GetZaxis()->SetRange(bin,bin);
    hKaonDepositsInAcceptance[bin-1] = (TH2D*) fHistK0EDepositsVsPtInAcceptance->Project3D("yx");
    ((TH2D*)hKaonDepositsInAcceptance[bin-1])->SetName(Form("hKaonDepositsInAcceptance%i",bin));
    hKaonDepositsInAcceptanceAvg[bin-1]  = (TProfile *) ((TH2D*)hKaonDepositsInAcceptance[bin-1])->ProfileX(Form("hKaonDepositsInAcceptanceAvg%i",bin));
    ((TH2D*)hKaonDepositsInAcceptance[bin-1])->GetYaxis()->SetRange(2,((TH2D*)hKaonDepositsInAcceptance[bin-1])->GetNbinsY());
    hKaonDepositsInAcceptanceLimAvg[bin-1]  = (TProfile *) ((TH2D*)hKaonDepositsInAcceptance[bin-1])->ProfileX(Form("hKaonDepositsInAcceptanceLimAvg%i",bin));

    fHistK0EGammaVsPtInAcceptance->GetZaxis()->SetRange(bin,bin);
    hKaonGammaInAcceptance[bin-1] = (TH2D*) fHistK0EGammaVsPtInAcceptance->Project3D("yx");
    ((TH2D*)hKaonGammaInAcceptance[bin-1])->SetName(Form("hKaonGammaInAcceptance%i",bin));
    hKaonGammaInAcceptanceAvg[bin-1]  = (TProfile *) ((TH2D*)hKaonGammaInAcceptance[bin-1])->ProfileX(Form("hKaonGammaInAcceptanceAvg%i",bin));
    ((TH2D*)hKaonGammaInAcceptance[bin-1])->GetYaxis()->SetRange(2,((TH2D*)hKaonGammaInAcceptance[bin-1])->GetNbinsY());
    hKaonGammaInAcceptanceLimAvg[bin-1]  = (TProfile *) ((TH2D*)hKaonGammaInAcceptance[bin-1])->ProfileX(Form("hKaonGammaInAcceptanceLimAvg%i",bin));


    fHistK0EDepositsVsPtOutOfAcceptance->GetZaxis()->SetRange(bin,bin);
    hKaonDepositsOutOfAcceptance[bin-1] = (TH2D*) fHistK0EDepositsVsPtOutOfAcceptance->Project3D("yx");
    ((TH2D*)hKaonDepositsOutOfAcceptance[bin-1])->SetName(Form("hKaonDepositsOutOfAcceptance%i",bin));
    hKaonDepositsOutOfAcceptanceAvg[bin-1]  = (TProfile *) ((TH2D*)hKaonDepositsOutOfAcceptance[bin-1])->ProfileX(Form("hKaonDepositsOutOfAcceptanceAvg%i",bin));
    ((TH2D*)hKaonDepositsOutOfAcceptance[bin-1])->GetYaxis()->SetRange(2,((TH2D*)hKaonDepositsOutOfAcceptance[bin-1])->GetNbinsY());
    hKaonDepositsOutOfAcceptanceLimAvg[bin-1]  = (TProfile *) ((TH2D*)hKaonDepositsOutOfAcceptance[bin-1])->ProfileX(Form("hKaonDepositsOutOfAcceptanceLimAvg%i",bin));

    fHistK0EGammaVsPtOutOfAcceptance->GetZaxis()->SetRange(bin,bin);
    hKaonGammaOutOfAcceptance[bin-1] = (TH2D*) fHistK0EGammaVsPtOutOfAcceptance->Project3D("yx");
    ((TH2D*)hKaonGammaOutOfAcceptance[bin-1])->SetName(Form("hKaonGammaOutOfAcceptance%i",bin));
    hKaonGammaOutOfAcceptanceAvg[bin-1]  = (TProfile *) ((TH2D*)hKaonGammaOutOfAcceptance[bin-1])->ProfileX(Form("hKaonGammaOutOfAcceptanceAvg%i",bin));
    ((TH2D*)hKaonGammaOutOfAcceptance[bin-1])->GetYaxis()->SetRange(2,((TH2D*)hKaonGammaOutOfAcceptance[bin-1])->GetNbinsY());
    hKaonGammaOutOfAcceptanceLimAvg[bin-1]  = (TProfile *) ((TH2D*)hKaonGammaOutOfAcceptance[bin-1])->ProfileX(Form("hKaonGammaOutOfAcceptanceLimAvg%i",bin));

  }
  int markers[] = {20,24,21,25,22, 26,23,32,33,27, 29};
  int colors[] = {TColor::kRed,TColor::kRed, TColor::kOrange-3, TColor::kOrange-3, TColor::kGreen+3,
		  TColor::kGreen+3, TColor::kBlue, TColor::kBlue, TColor::kViolet, TColor::kViolet,
		  TColor::kMagenta+3};
  for(int bin=0;bin<11;bin++){
    SetStyles((TH1*)hKaonDepositsInAcceptanceAvg[bin],markers[bin],colors[bin],"p_{T}","E");
    SetStyles((TH1*)hKaonDepositsInAcceptanceLimAvg[bin],markers[bin],colors[bin],"p_{T}","E");
    SetStyles((TH1*)hKaonGammaInAcceptanceAvg[bin],markers[bin],colors[bin],"p_{T}","E");
    SetStyles((TH1*)hKaonGammaInAcceptanceLimAvg[bin],markers[bin],colors[bin],"p_{T}","E");
    SetStyles((TH1*)hKaonDepositsOutOfAcceptanceAvg[bin],markers[bin],colors[bin],"p_{T}","E");
    SetStyles((TH1*)hKaonDepositsOutOfAcceptanceLimAvg[bin],markers[bin],colors[bin],"p_{T}","E");
    SetStyles((TH1*)hKaonGammaOutOfAcceptanceAvg[bin],markers[bin],colors[bin],"p_{T}","E");
    SetStyles((TH1*)hKaonGammaOutOfAcceptanceLimAvg[bin],markers[bin],colors[bin],"p_{T}","E");
  }
  float cuts[] = {0.00,0.05,0.10,0.15,0.20, 0.25,0.30,0.35,0.40,0.45, 0.50};
  TLegend *leg3 = new TLegend(0.454698,0.145161,0.620805,0.427419);
  leg3->SetFillStyle(0);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.03);
  leg3->SetTextSize(0.0456989);
  for(int bin=0;bin<6;bin++){
    leg3->AddEntry((TH1*)hKaonGammaInAcceptanceAvg[bin],Form("E_{cluster}>%1.2f",cuts[bin]));
  }
  TLegend *leg4 = new TLegend(0.696309,0.153226,0.862416,0.435484);
  leg4->SetFillStyle(0);
  leg4->SetFillColor(0);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.03);
  leg4->SetTextSize(0.0456989);
  for(int bin=6;bin<11;bin++){
    leg4->AddEntry((TH1*)hKaonGammaInAcceptanceAvg[bin],Form("E_{cluster}>%1.2f",cuts[bin]));
  }
  TH1F *frame = new TH1F("frame","frame",1,0,3);
  frame->SetMaximum(0.055);
  frame->SetMinimum(0.0);
  frame->GetXaxis()->SetTitle("p_{T}");
  frame->GetYaxis()->SetTitle("E");
  TH1F *frameClone = frame->Clone("frameClone");
  frameClone->SetMaximum(1.5);
 TCanvas *c5 = new TCanvas("c5","c5",600,400);
  c5->SetTopMargin(0.02);
  c5->SetRightMargin(0.02);
  c5->SetBorderSize(0);
  c5->SetFillColor(0);
  c5->SetFillColor(0);
  c5->SetBorderMode(0);
  c5->SetFrameFillColor(0);
  c5->SetFrameBorderMode(0);
  frame->Draw();
  for(int bin=0;bin<11;bin++){
      ((TH1*)hKaonGammaInAcceptanceAvg[bin])->Draw("same");
  }
  leg3->Draw();
  leg4->Draw();
  c5->SaveAs("/tmp/AvgEnergyDepositedInAcc.png");
 TCanvas *c6 = new TCanvas("c6","c6",600,400);
  c6->SetTopMargin(0.02);
  c6->SetRightMargin(0.02);
  c6->SetBorderSize(0);
  c6->SetFillColor(0);
  c6->SetFillColor(0);
  c6->SetBorderMode(0);
  c6->SetFrameFillColor(0);
  c6->SetFrameBorderMode(0);
  frameClone->Draw();
  for(int bin=0;bin<11;bin++){
      ((TH1*)hKaonGammaInAcceptanceLimAvg[bin])->Draw("same");
  }
  leg3->Draw();
  leg4->Draw();
  c6->SaveAs("/tmp/LimAvgEnergyDepositedInAcc.png");
 TCanvas *c7 = new TCanvas("c7","c7",600,400);
  c7->SetTopMargin(0.02);
  c7->SetRightMargin(0.02);
  c7->SetBorderSize(0);
  c7->SetFillColor(0);
  c7->SetFillColor(0);
  c7->SetBorderMode(0);
  c7->SetFrameFillColor(0);
  c7->SetFrameBorderMode(0);
  TH1F *frameclone4 = frame->Clone("frameclone4");
  frameclone4->SetMaximum(0.008);
  frameclone4->Draw();
  for(int bin=0;bin<11;bin++){
      ((TH1*)hKaonGammaOutOfAcceptanceAvg[bin])->Draw("same");
  }
  leg3->Draw();
  leg4->Draw();
  c7->SaveAs("/tmp/AvgEnergyDepositedOutAcc.png");
 TCanvas *c8 = new TCanvas("c8","c8",600,400);
  c8->SetTopMargin(0.02);
  c8->SetRightMargin(0.02);
  c8->SetBorderSize(0);
  c8->SetFillColor(0);
  c8->SetFillColor(0);
  c8->SetBorderMode(0);
  c8->SetFrameFillColor(0);
  c8->SetFrameBorderMode(0);
  frameClone->Draw();
  for(int bin=0;bin<11;bin++){
      ((TH1*)hKaonGammaOutOfAcceptanceLimAvg[bin])->Draw("same");
  }
  leg3->Draw();
  leg4->Draw();
  c8->SaveAs("/tmp/LimAvgEnergyDepositedOutAcc.png");

 TCanvas *c10 = new TCanvas("c10","c10",600,400);
  c10->SetTopMargin(0.02);
  c10->SetRightMargin(0.02);
  c10->SetBorderSize(0);
  c10->SetFillColor(0);
  c10->SetFillColor(0);
  c10->SetBorderMode(0);
  c10->SetFrameFillColor(0);
  c10->SetFrameBorderMode(0);
  TH1F *frameclone2 = (TH1F *)frame->Clone("frameclone2");
  frameclone2->SetMaximum(0.13);
  frameclone2->Draw();
  for(int bin=0;bin<11;bin++){
      ((TH1*)hKaonDepositsInAcceptanceAvg[bin])->Draw("same");
  }
  leg3->Draw();
  leg4->Draw();
  c10->SaveAs("/tmp/AvgEnergyDepositedCorrInAcc.png");
 TCanvas *c11 = new TCanvas("c11","c11",600,400);
  c11->SetTopMargin(0.02);
  c11->SetRightMargin(0.02);
  c11->SetBorderSize(0);
  c11->SetFillColor(0);
  c11->SetFillColor(0);
  c11->SetBorderMode(0);
  c11->SetFrameFillColor(0);
  c11->SetFrameBorderMode(0);
  TH1F *frameclone3 = (TH1F *)frame->Clone("frameclone2");
  frameclone3->SetMaximum(4);
  frameclone3->Draw();
  for(int bin=0;bin<11;bin++){
      ((TH1*)hKaonDepositsInAcceptanceLimAvg[bin])->Draw("same");
  }
  leg3->Draw();
  leg4->Draw();
  c11->SaveAs("/tmp/LimAvgEnergyDepositedCorrInAcc.png");
 TCanvas *c12 = new TCanvas("c12","c12",600,400);
  c12->SetTopMargin(0.02);
  c12->SetRightMargin(0.02);
  c12->SetBorderSize(0);
  c12->SetFillColor(0);
  c12->SetFillColor(0);
  c12->SetBorderMode(0);
  c12->SetFrameFillColor(0);
  c12->SetFrameBorderMode(0);
  TH1F *frameclone5 = frame->Clone("frameclone5");
  frameclone5->SetMaximum(0.013);
  frameclone5->Draw();
  for(int bin=0;bin<11;bin++){
      ((TH1*)hKaonDepositsOutOfAcceptanceAvg[bin])->Draw("same");
  }
  leg3->Draw();
  leg4->Draw();
  c12->SaveAs("/tmp/AvgEnergyDepositedCorrOutAcc.png");
 TCanvas *c13 = new TCanvas("c13","c13",600,400);
  c13->SetTopMargin(0.02);
  c13->SetRightMargin(0.02);
  c13->SetBorderSize(0);
  c13->SetFillColor(0);
  c13->SetFillColor(0);
  c13->SetBorderMode(0);
  c13->SetFrameFillColor(0);
  c13->SetFrameBorderMode(0);
  frameclone3->Draw();
  for(int bin=0;bin<11;bin++){
      ((TH1*)hKaonDepositsOutOfAcceptanceLimAvg[bin])->Draw("same");
  }
  leg3->Draw();
  leg4->Draw();
  c13->SaveAs("/tmp/LimAvgEnergyDepositedCorrOutAcc.png");

 TCanvas *c9 = new TCanvas("c9","c9",600,400);
  c9->SetTopMargin(0.02);
  c9->SetRightMargin(0.02);
  c9->SetBorderSize(0);
  c9->SetFillColor(0);
  c9->SetFillColor(0);
  c9->SetBorderMode(0);
  c9->SetFrameFillColor(0);
  c9->SetFrameBorderMode(0);
  ScaleByBinWidth((TH2D*) hKaonGammaInAcceptance[bin-1]);
  bin = 11;
  ((TH2D*) hKaonGammaInAcceptance[bin-1])->GetYaxis()->SetRange(2,((TH2D*) hKaonGammaInAcceptance[bin-1])->GetYaxis()->FindBin(3.0));
  ((TH2D*) hKaonGammaInAcceptance[bin-1])->GetXaxis()->SetRange(1,((TH2D*) hKaonGammaInAcceptance[bin-1])->GetXaxis()->FindBin(3.0));
  ((TH2D*) hKaonGammaInAcceptance[bin-1])->Draw("colz");


 TCanvas *c14 = new TCanvas("c14","c14",600,400);
  c14->SetTopMargin(0.02);
  c14->SetRightMargin(0.02);
  c14->SetBorderSize(0);
  c14->SetFillColor(0);
  c14->SetFillColor(0);
  c14->SetBorderMode(0);
  c14->SetFrameFillColor(0);
  c14->SetFrameBorderMode(0);
  ((TH2D*) hKaonGammaOutOfAcceptance[bin-1])->Draw("colz");

  bin = 7;
 TCanvas *c15 = new TCanvas("c15","c15",600,400);
  c15->SetTopMargin(0.02);
  c15->SetRightMargin(0.02);
  c15->SetBorderSize(0);
  c15->SetFillColor(0);
  c15->SetFillColor(0);
  c15->SetBorderMode(0);
  c15->SetFrameFillColor(0);
  c15->SetFrameBorderMode(0);
  ScaleByBinWidth((TH2D*) hKaonDepositsInAcceptance[bin-1]);

  ((TH2D*) hKaonDepositsInAcceptance[bin-1])->GetXaxis()->SetRange(1,((TH2D*) hKaonDepositsInAcceptance[bin-1])->GetXaxis()->FindBin(3.0));
  ((TH2D*) hKaonDepositsInAcceptance[bin-1])->GetYaxis()->SetRange(2,((TH2D*) hKaonDepositsInAcceptance[bin-1])->GetYaxis()->FindBin(3.0));
  ((TH2D*) hKaonDepositsInAcceptance[bin-1])->Draw("colz");

}
void ScaleByBinWidth(TH1 *histo){
  int nbins = histo->GetNbinsX();
  for(int i=1;i<=nbins;i++){
    float content = histo->GetBinContent(i);
    float scale = 1/histo->GetBinWidth(i);
    histo->SetBinContent(i,content*scale);
  }
}
void ScaleByBinWidth(TH2 *histo){
  int nbinsx = histo->GetNbinsX();
  int nbinsy = histo->GetNbinsY();
  for(int i=1;i<=nbinsx;i++){
    float scaleX = 1/histo->GetXaxis()->GetBinWidth(i);
    for(int j=1;j<=nbinsx;j++){
      float scaleY = 1/histo->GetXaxis()->GetBinWidth(j);
      float content = histo->GetBinContent(i,j);
      histo->SetBinContent(i,content*scaleX*scaleY);
    }
  }
}
void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle, Bool_t scale){
  histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetYaxis()->SetTitle(ytitle);
  if(scale) ScaleByBinWidth(histo);
}
