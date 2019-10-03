// -*- C++ -*-

const char* scans[] = {
  "root/4690/lumiRegion_Scan1X.root",
  "root/4690/lumiRegion_Scan1Y.root",
  "root/4690/lumiRegion_Scan2X.root",
  "root/4690/lumiRegion_Scan2Y.root",
  "root/4690/lumiRegion_ScanOffsetX.root",
  "root/4690/lumiRegion_ScanOffsetY.root",
  "root/4690/lumiRegion_ScanLSC.root",
};
const char* frameTitles[] = {
  "Scan1-X", "Scan1-Y",
  "Scan2-X", "Scan2-Y",
  "ScanOffset-X", "ScanOffset-Y",
  "LSC"
};

void PlotLumiRegion()
{
  gROOT->LoadMacro("Util.C"); // DrawFrame

  Int_t   fill   = 4690;
  Float_t xMin   = -0.14;
  Float_t xMax   =  0.14;
  TString pn     = "pdf/4690/LumiRegion_4690.pdf";

  TCanvas *c1 = new TCanvas("c1", "", 600, 600);
  c1->SaveAs(pn+"[");
  TText *tt = new TText;
  tt->SetTextSize(0.04);
  tt->SetTextAlign(22);
  TString label  = "ALICE Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV";
  TString xTitle = "Separation [mm]";
  for (Int_t i=0; i<7; ++i) {
    c1->Clear();
    c1->cd();
    tt->DrawTextNDC(0.5, 0.97, Form("%s - fill %d", frameTitles[i], fill));

    TPad *pad = new TPad("pad", "", 0,0, 1,0.95);
    pad->Draw();
    pad->Divide(3,3);
    TFile::Open(scans[i]);

    if (i == 6) {
      xMin = -0.5;
      xMax =  9.5;
      xTitle = "Step number";
      gStyle->SetOptFit(111);
      DrawFrame(pad->cd(1),  0.063, 0.075,  "<X>",  "", 0, xMin, xMax, label, xTitle); gX->Draw("P"); geX->Draw("PE");
      geX->Fit("pol1", "", "", 0, 4);
      gPad->Update();
      TPaveStats *ps = (TPaveStats*)geX->FindObject("stats");
      ps->SetX1NDC(0.35); ps->SetX2NDC(0.9);
      ps->SetY1NDC(0.2); ps->SetY2NDC(0.45);
      DrawFrame(pad->cd(2),  0.324, 0.337,  "<Y>",  "", 0, xMin, xMax, label, xTitle); gY->Draw("P"); geY->Draw("PE");
      geY->Fit("pol1", "", "", 5, 9);
      gPad->Update();
      TPaveStats *ps = (TPaveStats*)geY->FindObject("stats");
      ps->SetX1NDC(0.2); ps->SetX2NDC(0.7);
      ps->SetY1NDC(0.2); ps->SetY2NDC(0.45);
      DrawFrame(pad->cd(3), -2,     4,      "<Z>",  "", 0, xMin, xMax, label, xTitle); gZ->Draw("P"); geZ->Draw("PE");
    } else {
      DrawFrame(pad->cd(1),  0.067, 0.071,  "<X>",  "", 0, xMin, xMax, label, xTitle); gX->Draw("P"); geX->Draw("PE");
      DrawFrame(pad->cd(2),  0.329, 0.331,  "<Y>",  "", 0, xMin, xMax, label, xTitle); gY->Draw("P"); geY->Draw("PE");
      DrawFrame(pad->cd(3), -6,     6,      "<Z>",  "", 0, xMin, xMax, label, xTitle); gZ->Draw("P"); geZ->Draw("PE");
    }

    DrawFrame(pad->cd(4),  0.000, 0.0045, "<XX>", "", 0, xMin, xMax, label, xTitle); gSX->Draw("P"); geSX->Draw("PE");
    DrawFrame(pad->cd(5),  0.000, 0.0045, "<YY>", "", 0, xMin, xMax, label, xTitle); gSY->Draw("P"); geSY->Draw("PE");
    DrawFrame(pad->cd(6),  0.000, 7.500,  "<ZZ>", "", 0, xMin, xMax, label, xTitle); gSZ->Draw("P"); geSZ->Draw("PE");

    DrawFrame(pad->cd(7), -1,    +1,      "<XY>", "", 0, xMin, xMax, label, xTitle); gCXY->Draw("P"); geCXY->Draw("PE");
    gesY->SetLineColor(4);
    gesY->SetMarkerColor(4);
    DrawFrame(pad->cd(8),  0,     0.0005, "#color[2]{#mu_{X}'}, #color[4]{#mu_{Y}'}", "", 0, xMin, xMax, label, xTitle); gesX->Draw("PE"); gesY->Draw("PE");
    DrawFrame(pad->cd(9),  0.5,   2.5,    "k", "", 0, xMin, xMax, label, xTitle); gek->Draw("PE");
    c1->SaveAs(pn);
    gFile->Close();
  }
  c1->SaveAs(pn+"]");
}
