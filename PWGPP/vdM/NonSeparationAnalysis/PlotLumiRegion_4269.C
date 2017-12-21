// -*- C++ -*-

const char* scans[] = {
  "root/4269/lumiRegion_Scan1X.root",
  "root/4269/lumiRegion_Scan1Y.root",
  "root/4269/lumiRegion_Scan2X.root",
  "root/4269/lumiRegion_Scan2Y.root",
  "root/4269/lumiRegion_ScanOffsetX.root",
  "root/4269/lumiRegion_ScanOffsetY.root",
};
const char* frameTitles[] = {
  "Scan1-X", "Scan1-Y",
  "Scan2-X", "Scan2-Y",
  "ScanOffset-X", "ScanOffset-Y",
};

void PlotLumiRegion()
{
  gROOT->LoadMacro("Util.C"); // DrawFrame

  Int_t   fill   = 4269;
  Float_t xMin    = -0.6;
  Float_t xMax    =  0.6;
  TString pn     = "pdf/4269/LumiRegion_4269.pdf";

  TCanvas *c1 = new TCanvas("c1", "", 600, 600);
  c1->SaveAs(pn+"[");
  TText *tt = new TText;
  tt->SetTextSize(0.04);
  tt->SetTextAlign(22);
  TString label = "ALICE pp #sqrt{#it{s}}=13 TeV";
  for (Int_t i=0; i<6; ++i) {
    c1->Clear();
    c1->cd();
    tt->DrawTextNDC(0.5, 0.97, Form("%s - fill %d", frameTitles[i], fill));

    TPad *pad = new TPad("pad", "", 0,0, 1,0.95);
    pad->Draw();
    pad->Divide(3,3);
    TFile::Open(scans[i]);
    gesY->SetLineColor(4);
    gesY->SetMarkerColor(4);

    DrawFrame(pad->cd(1),  0.068, 0.082,  "<X>",  "", 0, xMin, xMax, label); gX->Draw("P"); geX->Draw("PE");
    DrawFrame(pad->cd(2),  0.532, 0.545,  "<Y>",  "", 0, xMin, xMax, label); gY->Draw("P"); geY->Draw("PE");
    DrawFrame(pad->cd(3), -5,     5.0,    "<Z>",  "", 0, xMin, xMax, label); gZ->Draw("P"); geZ->Draw("PE");

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
