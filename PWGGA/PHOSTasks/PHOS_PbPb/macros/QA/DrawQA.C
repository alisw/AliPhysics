#include "TCanvas.h" // needed for some reason.

TFile* file;
const char* prefixToName = "imgs/";
const int kNCents = 1;

void Draw(const char* name, const char* options = "", double yFrom=0., double yTo=-1.)
{
  TH1* hist = ((TH1*)file->Get(name))->Clone();
  //hist->SetAxisRange(0, 30 );
  hist->GetXaxis()->SetLabelSize(0.011);
  //hist->GetXaxis()->SetTitle("Run");
  //hist->GetXaxis()->SetRange(1, );
  if( yFrom < yTo ) {
    if(yFrom > hist->GetMinimum())
      Printf("in hist %s, yFrom (%f), is larger then hist min (%f)", name, yFrom, hist->GetMinimum());
    if(yTo < hist->GetMaximum())
      Printf("in hist %s, yTo (%f), is smaller then hist max (%f)", name, yTo, hist->GetMaximum());

    hist->GetYaxis()->SetRangeUser(yFrom, yTo);
  }

  //if( ! TString(options).Contains("same") )
  TCanvas* canv = new TCanvas;
  //canv->SetGrid();
  //hist->GetYaxis()->SetNdivisions(16);

  canv->Divide(1,2);

  hist->GetXaxis()->SetLabelSize(0.051);

  if( TString(options).Contains("LINFIT") )
    hist->Fit("pol0", "Q");
  
  canv->cd(1);
  if( TString(name).Contains("grChi2RP") )
    gPad->SetLogy();
  hist->GetXaxis()->SetRange(1, hist->GetNbinsX()/2);
  hist->DrawCopy(options);

  canv->cd(2);
  if( TString(name).Contains("grChi2RP") )
    gPad->SetLogy();
  hist->GetXaxis()->SetRange(hist->GetNbinsX()/2+1,200);
  hist->DrawCopy(options);


  canv->SaveAs(Form("%s%s.pdf", prefixToName, hist->GetName()));
  canv->SaveAs(Form("%s%s.png", prefixToName, hist->GetName()));
  delete hist;
}

void DrawPID()
{
  for(int cent=0; cent<kNCents; ++cent) {
    TString key = Form("nPhotPID_cen%i", cent);
    
    TH1* grNPhotAllcore = (TH1*)file->Get(Form("grNPhotAllcore_cen%d", cent))->Clone();
    TH1* grNPhotCPVcore = (TH1*)file->Get(Form("grNPhotCPVcore_cen%d", cent))->Clone();
    TH1* grNPhotDisp2core = (TH1*)file->Get(Form("grNPhotDisp2core_cen%d", cent))->Clone();
    TH1* grNPhotBoth2core = (TH1*)file->Get(Form("grNPhotBoth2core_cen%d", cent))->Clone();

    grNPhotAllcore->GetXaxis()->SetLabelSize(0.045);
    grNPhotAllcore->GetYaxis()->SetRangeUser(0, grNPhotAllcore->GetMaximum()*1.2);

    grNPhotAllcore  ->SetMarkerColor(kBlack);
    grNPhotCPVcore ->SetMarkerColor(kCyan+1);
    grNPhotDisp2core->SetMarkerColor(kBlue);
    grNPhotBoth2core  ->SetMarkerColor(kOrange+1);

    TCanvas* canv = new TCanvas(key.Data(), key.Data(), 1024, 768);
    canv->Divide(1,2);

    canv->cd(1);
    grNPhotAllcore->SetTitle("#LTN_{clusters}^{PID}#GT");
    gPad->SetGridy();
    grNPhotAllcore->GetXaxis()->SetRange(0, 85);
    grNPhotAllcore->DrawCopy();
    grNPhotCPVcore->DrawCopy("same");
    grNPhotDisp2core->DrawCopy("same");
    grNPhotBoth2core->DrawCopy("same");

    canv->cd(2);
    gPad->SetGridy();
    grNPhotAllcore->SetTitle("");
    grNPhotAllcore->GetXaxis()->SetRange(85, 200);
    grNPhotAllcore->DrawCopy();
    grNPhotCPVcore->DrawCopy("same");
    grNPhotDisp2core->DrawCopy("same");
    grNPhotBoth2core->DrawCopy("same");

    canv->cd(1);
    leg = new TLegend(0.9,0.7,0.99,0.99);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->AddEntry(Form("grNPhotAllcore_cen%d",cent),"All","lP");
    leg->AddEntry(Form("grNPhotCPVcore_cen%d",cent),"CPVcore","lP");
    leg->AddEntry(Form("grNPhotDisp2core_cen%d",cent),"Disp2core","lP");
    leg->AddEntry(Form("grNPhotBoth2core_cen%d",cent),"Both2core","lP");
    
    leg->Draw();
  
    canv->SaveAs(Form("%s%s_squeezed.pdf", prefixToName, key.Data()));
    canv->SaveAs(Form("%s%s_squeezed.png", prefixToName, key.Data()));
  }
}

void DrawCPVRatio()
{
  for(int cent=0; cent<kNCents; ++cent) {
    TH1* grNPhotAll = (TH1*)file->Get(Form("grNPhotAll_cen0", cent))->Clone();
    TH1* grNPhotCPV = (TH1*)file->Get(Form("grNPhotCPV_cen0", cent))->Clone();
    TH1* grNPhotCPV2 = (TH1*)file->Get(Form("grNPhotCPV2_cen0", cent))->Clone();

    grNPhotCPV->Divide(grNPhotAll);
    grNPhotCPV->SetTitle(Form("%s / %s", grNPhotCPV->GetTitle(), grNPhotAll->GetTitle()));
    grNPhotCPV->GetYaxis()->SetRangeUser(0.7,0.85);

    grNPhotCPV2->Divide(grNPhotAll);

    TCanvas* canv = new TCanvas;
    canv->Divide(1,2);
    canv->cd(1);
    grNPhotCPV->GetXaxis()->SetRange(0, 85);
    grNPhotCPV->DrawCopy();

    canv->cd(2);
    grNPhotCPV->SetTitle("");
    grNPhotCPV->GetXaxis()->SetRange(85, 200);
    grNPhotCPV->DrawCopy();
  
    canv->SaveAs(Form("%sCPVtoAllRatio_cen%d.pdf", prefixToName, cent));
    canv->SaveAs(Form("%sCPVtoAllRatio_cen%d.png", prefixToName, cent));
  }
}

void DrawNPhotAllAndHigh()
{
  for(int cent=0; cent<kNCents; ++cent) {
    TH1* grNPhotAll = (TH1*)file->Get(Form("grNPhotAll_cen%d", cent))->Clone();
    TH1* grNPhotAllHigh = (TH1*)file->Get(Form("grNPhotAllHigh_cen%d", cent))->Clone();
    TH1* grNPhotAllcoreHigh = (TH1*)file->Get(Form("grNPhotAllcoreHigh_cen%d", cent))->Clone();

    double sizeAll = grNPhotAll->Integral();
    double sizeAllHigh = grNPhotAllHigh->Integral();
    int scale = TMath::Nint(sizeAll / sizeAllHigh);
    grNPhotAllHigh->Scale(scale);
    grNPhotAllcoreHigh->Scale(scale);
  
    grNPhotAllHigh->SetMarkerColor(kRed);
    grNPhotAllHigh->SetLineColor(kRed);
    grNPhotAllcoreHigh->SetMarkerColor(kGreen+1);
    grNPhotAllcoreHigh->SetLineColor(kGreen+1);

    TCanvas* canv = new TCanvas;
    canv->Divide(1,2);
  
    canv->cd(1);
    grNPhotAll->SetTitle("#LTN_{clusters}#GT");
    grNPhotAll->GetXaxis()->SetRange(0, 85);
    grNPhotAll->GetYaxis()->SetRange(15, 40);
    //grNPhotAll->GetYaxis()->SetRangeUser(15, 45);
    grNPhotAll->DrawCopy();
    grNPhotAllHigh->DrawCopy("same");
    grNPhotAllcoreHigh->DrawCopy("same");

    canv->cd(2);
    grNPhotAll->GetXaxis()->SetRange(85, 200);
    grNPhotAll->DrawCopy();
    grNPhotAllHigh->DrawCopy("same");
    grNPhotAllcoreHigh->DrawCopy("same");

    canv->cd(1);
    leg = new TLegend(0.8, 0.15, 0.99, 0.4);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->AddEntry(grNPhotAll, Form("All"),"lP");
    leg->AddEntry(grNPhotAllHigh, Form("AllHigh * %d", scale),"lP");
    leg->AddEntry(grNPhotAllcoreHigh, Form("AllcoreHigh * %d", scale),"lP");
    leg->Draw();

    canv->SaveAs( Form("%s%s.pdf", prefixToName, Form("nPhotAllAndHigh_cen%d", cent)));
    canv->SaveAs( Form("%s%s.png", prefixToName, Form("nPhotAllAndHigh_cen%d", cent)));
  }
}

void DrawPIDRatiosHighCore(const char* pidNames[], int nPids, const char* high)
{
  int kNColors = 8;
  const Int_t colors[8] = {kBlack, kRed-1, kRed+1, kBlue, kCyan, kGreen+3, kYellow+1, kMagenta};
  int currentColorID = 99999;
  Int_t color;
  
  for(int cent=0; cent<kNCents; ++cent) {
    TH1* hAll = (TH1*)file->Get( Form("grNPhot%s%s_cen%d", pidNames[0], high, cent) )->Clone();

    leg = new TLegend(0.91,0.6,0.99,0.99);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);

    TCanvas* canv = new TCanvas;
    canv->Divide(1,2);
    char* same = "";
    for(int ipid = 1; ipid < nPids; ++ipid) {
      TString name(Form("grNPhot%s%s_cen%d", pidNames[ipid], high, cent));
      TH1* hPID = (TH1*)file->Get( name.Data() )->Clone();
      hPID->Divide(hAll);

      if( ++currentColorID < kNColors )
	color = colors[currentColorID] ;
      else
	color = colors[ currentColorID = 0 ]; 

      hPID->SetMarkerColor(color);
      hPID->SetLineColor(color);
      hPID->GetYaxis()->SetRangeUser(0, 1.);

      leg->AddEntry( hPID , pidNames[ipid], "lP");

      canv->cd(1);
      hPID->SetTitle( Form("#LTN_{clusters}^{PID}#GT / #LTN_{clusters}^{%s}#GT, cent=%d, %s", pidNames[0], cent, high) );
      hPID->GetXaxis()->SetRange(0, 85);
      hPID->DrawCopy(same);

      canv->cd(2);
      hPID->SetTitle("");
      hPID->GetXaxis()->SetRange(85, 200);
      hPID->DrawCopy(same);      

      same = "same";
    }
    canv->cd(1);
    leg->Draw();
    TString fn(Form("CPVtoAllRatio%s_cen%d", high, cent));

    canv->SaveAs(Form("%s%s.pdf", prefixToName, fn.Data()));
    canv->SaveAs(Form("%s%s.png", prefixToName, fn.Data()));
  }
}

void DrawPIDRatios()
{
  const int nn = 8;
  const int nc = 4;
  const char* kPIDNames[nn] = {"All", "Allwou", "CPV", "CPV2", "Disp", "Disp2", "Dispwou", "Both"};
  const char* kPIDNamesCore[nc] = {"Allcore", "CPVcore", "Dispcore", "Bothcore"};
  const char* fillHigh[2] = {"", "High"};
  for(int ihigh = 0; ihigh < 2; ++ihigh) {
    DrawPIDRatiosHighCore(kPIDNames, nn, fillHigh[ihigh]);
    DrawPIDRatiosHighCore(kPIDNamesCore, nc, fillHigh[ihigh]);
  }
}

void DrawQA()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  file = TFile::Open("outputQA.root", "read");

  Draw("grVtxZ10Cent", "", 0.7, 1.);
  Draw("grNCellsM1", "E");
  Draw("grNCellsM2");
  Draw("grNCellsM3");
  Draw("grECluster", "", 0.5, 0.7);
  Draw("grNCluster", "", 0, 40);
  Draw("grNTracks0", "", 0 , 12000);
  Draw("grNPhotAll_cen0", "", 0, 40);
  Draw("grNPhotAllcore_cen0", "", 0, 40);
  Draw("grNPhotAllwou_cen0", "", 0, 40);
  Draw("grNPhotDisp_cen0", "", 0, 40);
  Draw("grNPhotDisp2_cen0", "", 0, 40);
  Draw("grNPhotDispwou_cen0", "", 0, 40);
  Draw("grNPhotCPV_cen0", "", 0, 40);
  Draw("grNPhotCPV2_cen0", "", 0, 40);
  Draw("grNPhotBoth_cen0", "", 0, 40);
  Draw("grEnAll_cen0", "", 0.4, 0.7);
  Draw("grEnAllcore_cen0", "", 0.4, 0.7);
  Draw("grEnAllwou_cen0", "", 0.4, 0.7);
  Draw("grEnDisp_cen0", "", 0.4, 0.7);
  Draw("grEnDisp2_cen0", "", 0.4, 0.7);
  Draw("grEnDispcore_cen0", "", 0.4, 0.7);
  Draw("grEnDispwou_cen0", "", 0.4, 0.7);
  Draw("grEnCPV_cen0", "", 0.4, 0.7);
  Draw("grEnCPVcore_cen0", "", 0.4, 0.7);
  Draw("grEnCPV2_cen0", "", 0.4, 0.7);
  Draw("grEnBoth_cen0", "", 0.4, 0.7);
  Draw("grEnBothcore_cen0", "", 0.4, 0.7);


  DrawPID();
  DrawCPVRatio();
  DrawNPhotAllAndHigh();
  DrawPIDRatios();

  Draw("grMPi0", "LINFIT", 0.13, 0.15 );
  Draw("grWPi0", "LINFIT");
  Draw("grNPi0", "LINFIT");

  Draw("grSERPV0Aflat", "", 0, 0.4);
  Draw("grSERPV0Cflat", "", 0, 0.4);
  Draw("grSERPTPCflat", "", 0, 0.4);

  Draw("grChi2RPV0A", "", 0.1, 500);
  Draw("grChi2RPV0C", "", 0.1, 500);
  Draw("grChi2RPTPC", "", 0.1, 500);
  Draw("grChi2RPV0Aflat", "", 0.1, 500);
  Draw("grChi2RPV0Cflat", "", 0.1, 500);
  Draw("grChi2RPTPCflat", "", 0.1, 500);

  Draw("grSinRPV0A1", "", 0, 0.7);
  Draw("grSinRPV0A2", "", 0, 10);
  Draw("grSinRPV0A3", "", -TMath::Pi(), TMath::Pi());
  Draw("grSinRPV0C1", "", 0, 0.7);
  Draw("grSinRPV0C2", "", 0, 10);
  Draw("grSinRPV0C3", "", -TMath::Pi(), TMath::Pi());
  Draw("grSinRPTPC1", "", 0, 0.7);
  Draw("grSinRPTPC2", "", 0, 10);
  Draw("grSinRPTPC3", "", -TMath::Pi(), TMath::Pi());
 
  file->Close();
}
