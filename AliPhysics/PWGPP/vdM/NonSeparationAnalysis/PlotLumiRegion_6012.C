// -*- C++ -*-

const char* scans[] = {
  "root/6012/lumiRegion_Scan1X.root",
  "root/6012/lumiRegion_Scan1Y.root",
  "root/6012/lumiRegion_Scan2X.root",
  "root/6012/lumiRegion_Scan2Y.root",
  "root/6012/lumiRegion_ScanLSC.root",
  "root/6012/lumiRegion_ScanOffsetX.root",
  "root/6012/lumiRegion_ScanOFfsetY.root",
};
const char* frameTitles[] = {
  "Scan1-X", "Scan1-Y",
  "Scan2-X", "Scan2-Y",
  "LSC",
  "Offset-X", "Offset-Y"
};

void PlotLumiRegion()
{
  gROOT->LoadMacro("Util.C"); // DrawFrame

  Int_t   fill   = 6012;
  Float_t xMin    = -0.62;
  Float_t xMax    =  0.62;
  TString pn     = "pdf/6012/LumiRegion_6012.pdf";

  TCanvas *c1 = new TCanvas("c1", "", 600, 600);
  c1->SaveAs(pn+"[");
  TText *tt = new TText;
  tt->SetTextSize(0.04);
  tt->SetTextAlign(22);
  TString label = "ALICE pp #sqrt{#it{s}}=5 TeV";
  for (Int_t i=0; i<7; ++i) {
    c1->Clear();
    c1->cd();
    tt->DrawTextNDC(0.5, 0.97, Form("%s - fill %d", frameTitles[i], fill));

    TPad *pad = new TPad("pad", "", 0,0, 1,0.95);
    pad->Draw();
    pad->Divide(3,3);
    TFile::Open(scans[i]);
    gesY->SetLineColor(4);
    gesY->SetMarkerColor(4);

    if (i == 4) {
      xMin = -0.5;
      xMax =  13.5;
      xTitle = "Step number";
      gStyle->SetOptFit(111);

      DrawFrame(pad->cd(1),  0.02, 0.14,  "<X>",  "", 0, xMin, xMax, label, xTitle); gX->Draw("P"); geX->Draw("PE");
      geX->Fit("pol1", "", "", 1.5, 5.5);
      gPad->Update();
      TPaveStats *ps = (TPaveStats*)geX->FindObject("stats");
      ps->SetX1NDC(0.35); ps->SetX2NDC(0.9);
      ps->SetY1NDC(0.2); ps->SetY2NDC(0.45);

      DrawFrame(pad->cd(2),  0.41, 0.51,  "<Y>",  "", 0, xMin, xMax, label, xTitle); gY->Draw("P"); geY->Draw("PE");
      geY->Fit("pol1", "", "", 6.5, 12.5);
      gPad->Update();
      TPaveStats *ps = (TPaveStats*)geY->FindObject("stats");
      ps->SetX1NDC(0.2); ps->SetX2NDC(0.7);
      ps->SetY1NDC(0.2); ps->SetY2NDC(0.45);

      DrawFrame(pad->cd(3),  0.3,   0.8,     "<Z>", "", 0, xMin, xMax, label, xTitle); gZ->Draw("P"); geZ->Draw("PE");
    } else {
      if (i > 4) {
         xMin    = -0.45;
         xMax    =  0.45;
      }
      DrawFrame(pad->cd(1),  0.068, 0.082,  "<X>",  "", 0, xMin, xMax, label); gX->Draw("P"); geX->Draw("PE");
      DrawFrame(pad->cd(2),  0.450, 0.460,  "<Y>",  "", 0, xMin, xMax, label); gY->Draw("P"); geY->Draw("PE");
      DrawFrame(pad->cd(3), -7,     7.0,    "<Z>",  "", 0, xMin, xMax, label); gZ->Draw("P"); geZ->Draw("PE");
    }

    DrawFrame(pad->cd(4),  0.000, 0.0135, "<XX>", "", 0, xMin, xMax, label); gSX->Draw("P"); geSX->Draw("PE");
    DrawFrame(pad->cd(5),  0.000, 0.0135, "<YY>", "", 0, xMin, xMax, label); gSY->Draw("P"); geSY->Draw("PE");
    DrawFrame(pad->cd(6),  0.000, 6.500,  "<ZZ>", "", 0, xMin, xMax, label); gSZ->Draw("P"); geSZ->Draw("PE");

    DrawFrame(pad->cd(7), -1,    +1,      "<XY>", "", 0, xMin, xMax, label); gCXY->Draw("P"); geCXY->Draw("PE");
    DrawFrame(pad->cd(8),  0,     0.0005, "#color[2]{#mu_{X}'}, #color[4]{#mu_{Y}'}", "", 0, xMin, xMax, label); gesX->Draw("PE"); gesY->Draw("PE");
    DrawFrame(pad->cd(9),  0.5,   2.5,    "k", "", 0, xMin, xMax, label); gek->Draw("PE");
    c1->SaveAs(pn);
  }
  c1->SaveAs(pn+"]");
}
