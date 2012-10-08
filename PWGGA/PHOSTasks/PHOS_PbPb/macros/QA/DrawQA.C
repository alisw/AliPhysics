#include "TCanvas.h" // needed for some reason.

TFile* file;
const char* prefixToName = "imgs/";
const char* appendToName = ".pdf";
const int kNCents = 1;

void Draw(const char* name, const char* options = "", double yFrom=0., double yTo=-1.)
{
  TH1* hist = ((TH1*)file->Get(name))->Clone();
  //hist->SetAxisRange(0, 30 );
  hist->GetXaxis()->SetLabelSize(0.011);
  //hist->GetXaxis()->SetTitle("Run");
  hist->GetXaxis()->SetRange(1,84);
  if( yFrom < yTo )
    hist->GetYaxis()->SetRangeUser(yFrom, yTo);


  //if( ! TString(options).Contains("same") )
  TCanvas* canv = new TCanvas;
  canv->SetGrid();
  hist->GetYaxis()->SetNdivisions(16);

  canv->Divide(1,2);

  hist->GetXaxis()->SetLabelSize(0.051);

  canv->cd(1);
  hist->GetXaxis()->SetRange(1,84);
  hist->DrawCopy(options);

  canv->cd(2);
  hist->GetXaxis()->SetRange(85,200);
  hist->DrawCopy(options);


  canv->SaveAs(Form("%s%s%s", prefixToName, hist->GetName(), appendToName ));
  delete hist;
}

void DrawPID()
{
  for(int cent=0; cent<kNCents; ++cent) {
    TH1* grNPhotAll = (TH1*)file->Get(Form("grNPhotAll_cen%d", cent))->Clone();
    TH1* grNPhotDisp = (TH1*)file->Get(Form("grNPhotDisp_cen%d", cent))->Clone();
    TH1* grNPhotDisp2 = (TH1*)file->Get(Form("grNPhotDisp2_cen%d", cent))->Clone();
    TH1* grNPhotCPV = (TH1*)file->Get(Form("grNPhotCPV_cen%d", cent))->Clone();
    TH1* grNPhotCPV2 = (TH1*)file->Get(Form("grNPhotCPV2_cen%d", cent))->Clone();

    grNPhotAll->GetXaxis()->SetLabelSize(0.045);
    grNPhotAll->GetYaxis()->SetRangeUser(0, 40);

    grNPhotAll  ->SetMarkerColor(kBlack);
    grNPhotDisp ->SetMarkerColor(kCyan+1);
    grNPhotDisp2->SetMarkerColor(kBlue);
    grNPhotCPV  ->SetMarkerColor(kOrange+1);
    grNPhotCPV2 ->SetMarkerColor(kRed);

    TCanvas* canv = new TCanvas;
    canv->Divide(1,2);

    canv->cd(1);
    grNPhotAll->SetTitle("#LTN_{clusters}^{PID}#GT");
    grNPhotAll->GetXaxis()->SetRange(0, 84);
    grNPhotAll->DrawCopy();
    grNPhotDisp->DrawCopy("same");
    grNPhotDisp2->DrawCopy("same");
    grNPhotCPV->DrawCopy("same");
    grNPhotCPV2->DrawCopy("same");

    canv->cd(2);
    grNPhotAll->SetTitle("");
    grNPhotAll->GetXaxis()->SetRange(85, 200);
    grNPhotAll->DrawCopy();
    grNPhotDisp->DrawCopy("same");
    grNPhotDisp2->DrawCopy("same");
    grNPhotCPV->DrawCopy("same");
    grNPhotCPV2->DrawCopy("same");

    canv->cd(1);
    leg = new TLegend(0.9,0.7,0.99,0.99);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->AddEntry(Form("grNPhotAll_cen%d",cent),"All","lP");
    leg->AddEntry(Form("grNPhotCPV_cen%d",cent),"CPV","lP");
    leg->AddEntry(Form("grNPhotCPV2_cen%d",cent),"CPV2","lP");
    leg->AddEntry(Form("grNPhotDisp_cen%d",cent),"Disp","lP");
    leg->AddEntry(Form("grNPhotDisp2_cen%d",cent),"Disp2","lP");
    leg->Draw();
  
    canv->SaveAs(Form("%s%s%s", prefixToName, Form("nPhotPID_cen%d", cent), appendToName ));
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
    grNPhotCPV->GetXaxis()->SetRange(0, 84);
    grNPhotCPV->DrawCopy();

    canv->cd(2);
    grNPhotCPV->SetTitle("");
    grNPhotCPV->GetXaxis()->SetRange(85, 200);
    grNPhotCPV->DrawCopy();
  
    canv->SaveAs(Form("%s%s%s", prefixToName, Form("CPVtoAllRatio_cen%d", cent), appendToName ));
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
    grNPhotAll->GetXaxis()->SetRange(0, 84);
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

    canv->SaveAs( Form("%s%s%s", prefixToName, Form("nPhotAllAndHigh_cen%d", cent), appendToName ));
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
      hPID->GetXaxis()->SetRange(0, 84);
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
    canv->SaveAs(Form("%s%s%s", prefixToName, fn.Data(), appendToName ));
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

  file = TFile::Open("outputQA.root", "read");

  // Draw("grVtxZ10Cent", "", 0.7, 1.);
  // Draw("grNCellsM1", "E");
  // Draw("grNCellsM2");
  // Draw("grNCellsM3");
  // Draw("grECluster", "", 0.5, 0.7);
  // Draw("grNCluster", "", 0, 40);
  // Draw("grNTracks0", "", 0 , 12000);
  // Draw("grNPhotAll_cen0", "", 0, 40);
  // Draw("grNPhotAllcore_cen0", "", 0, 40);
  // Draw("grNPhotAllwou_cen0", "", 0, 40);
  // Draw("grNPhotDisp_cen0", "", 0, 40);
  // Draw("grNPhotDisp2_cen0", "", 0, 40);
  // Draw("grNPhotDispwou_cen0", "", 0, 40);
  // Draw("grNPhotCPV_cen0", "", 0, 40);
  // Draw("grNPhotCPV2_cen0", "", 0, 40);
  // Draw("grNPhotBoth_cen0", "", 0, 40);
  // Draw("grEnAll_cen0", "", 0.4, 0.7);
  // Draw("grEnAllcore_cen0", "", 0.4, 0.7);
  // Draw("grEnAllwou_cen0", "", 0.4, 0.7);
  // Draw("grEnDisp_cen0", "", 0.4, 0.7);
  // Draw("grEnDisp2_cen0", "", 0.4, 0.7);
  // Draw("grEnDispcore_cen0", "", 0.4, 0.7);
  // Draw("grEnDispwou_cen0", "", 0.4, 0.7);
  // Draw("grEnCPV_cen0", "", 0.4, 0.7);
  // Draw("grEnCPVcore_cen0", "", 0.4, 0.7);
  // Draw("grEnCPV2_cen0", "", 0.4, 0.7);
  // Draw("grEnBoth_cen0", "", 0.4, 0.7);
  // Draw("grEnBothcore_cen0", "", 0.4, 0.7);


  // DrawPID();
  // DrawCPVRatio();
  DrawNPhotAllAndHigh();
  //DrawPIDRatios();

  Draw("grMPi0");
  Draw("grWPi0");
  Draw("grNPi0");

  file->Close();
}
