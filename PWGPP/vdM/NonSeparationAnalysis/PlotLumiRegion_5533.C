// -*- C++ -*-

const char* scans[] = {
  "root/5533/lumiRegion_Scan1X.root",
  "root/5533/lumiRegion_Scan1Y.root",
  "root/5533/lumiRegion_Scan2X.root",
  "root/5533/lumiRegion_Scan2Y.root",
  "root/5533/lumiRegion_ScanLSC.root",
};
const char* frameTitles[] = {
  "Scan1-X", "Scan1-Y",
  "Scan2-X", "Scan2-Y",
  "LSC"
};

void PlotLumiRegion()
{
  gROOT->LoadMacro("Util.C"); // DrawFrame

  Int_t   fill   = 5533;
  Float_t xMin   = -0.199;
  Float_t xMax   =  0.199;
  TString xTitle = "Separation [mm]";
  TString pn     = "pdf/5533/LumiRegion_5533.pdf";

  TCanvas *c1 = new TCanvas("c1", "", 600, 600);
  c1->SaveAs(pn+"[");
  TText *tt = new TText;
  tt->SetTextSize(0.04);
  tt->SetTextAlign(22);
  TString label = "ALICE p-Pb #sqrt{#it{s}_{NN}}=8.2 TeV";
  for (Int_t i=0; i<5; ++i) {
    Printf("i=%d", i);
    c1->Clear();
    c1->cd();
    tt->DrawTextNDC(0.5, 0.97, Form("%s - fill %d", frameTitles[i], fill));

    TPad *pad = new TPad("pad", "", 0,0, 1,0.95);
    pad->Draw();
    pad->Divide(3,3);
    TFile::Open(scans[i]);

    if (i == 4) {
      xMin = -0.5;
      xMax =  9.5;
      xTitle = "Step number";
      gStyle->SetOptFit(111);
      DrawFrame(pad->cd(1),  0.070, 0.089,  "<X>",  "", 0, xMin, xMax, label, xTitle); gX->Draw("P"); geX->Draw("PE");
      geX->Fit("pol1", "", "", 0, 4);
      gPad->Update();
      TPaveStats *ps = (TPaveStats*)geX->FindObject("stats");
      ps->SetX1NDC(0.35); ps->SetX2NDC(0.9);
      ps->SetY1NDC(0.2); ps->SetY2NDC(0.45);

      DrawFrame(pad->cd(2),  0.350, 0.372,  "<Y>",  "", 0, xMin, xMax, label, xTitle); gY->Draw("P"); geY->Draw("PE");
      geY->Fit("pol1", "", "", 5, 9);
      gPad->Update();
      TPaveStats *ps = (TPaveStats*)geY->FindObject("stats");
      ps->SetX1NDC(0.2); ps->SetX2NDC(0.7);
      ps->SetY1NDC(0.2); ps->SetY2NDC(0.45);

      DrawFrame(pad->cd(3), 0.5,      2.5,    "<Z>",  "", 0, xMin, xMax, label, xTitle); gZ->Draw("P"); geZ->Draw("PE");
    } else {
      DrawFrame(pad->cd(1),  0.077, 0.084,  "<X>",  "", 0, xMin, xMax, label, xTitle); gX->Draw("P"); geX->Draw("PE");
      DrawFrame(pad->cd(2),  0.352, 0.366,  "<Y>",  "", 0, xMin, xMax, label, xTitle); gY->Draw("P"); geY->Draw("PE");
      DrawFrame(pad->cd(3), -1,     4.5,    "<Z>",  "", 0, xMin, xMax, label, xTitle); gZ->Draw("P"); geZ->Draw("PE");
    }

    DrawFrame(pad->cd(4),  0.000, 0.007, "<XX>", "", 0, xMin, xMax, label, xTitle); gSX->Draw("P"); geSX->Draw("PE");
    DrawFrame(pad->cd(5),  0.000, 0.007, "<YY>", "", 0, xMin, xMax, label, xTitle); gSY->Draw("P"); geSY->Draw("PE");
    DrawFrame(pad->cd(6),  0.000, 7.500,  "<ZZ>", "", 0, xMin, xMax, label, xTitle); gSZ->Draw("P"); geSZ->Draw("PE");

    DrawFrame(pad->cd(7), -1,    +1,      "<XY>", "", 0, xMin, xMax, label, xTitle); gCXY->Draw("P"); geCXY->Draw("PE");
    gesY->SetLineColor(4);
    gesY->SetMarkerColor(4);
    DrawFrame(pad->cd(8),  0,     0.0005, "#color[2]{#mu_{X}'}, #color[4]{#mu_{Y}'}", "", 0, xMin, xMax, label, xTitle); gesX->Draw("PE"); gesY->Draw("PE");
    DrawFrame(pad->cd(9),  0.5,   2.5,    "k",    "", 0, xMin, xMax, label, xTitle); gek->Draw("PE");
    c1->SaveAs(pn);
  }
  c1->SaveAs(pn+"]");
}
