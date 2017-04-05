// -*- C++ -*-

const char* scans[] = {
  "root/4690/lumiRegion_Scan1X.root",
  "root/4690/lumiRegion_Scan1Y.root",
  "root/4690/lumiRegion_ScanOffsetX.root",
  "root/4690/lumiRegion_ScanOffsetY.root"
};
const char* frameTitles[] = {
  "Scan1-X", "Scan1-Y",
  "ScanOffset-X", "ScanOffset-Y"
};

void PlotLumiRegion_4690()
{
  gROOT->LoadMacro("Util.C"); // DrawFrame

  Int_t   fill   = 4690;
  Float_t maxSep = 0.14;
  TString pn     = "pdf/4690/LumiRegion_4690.pdf";

  TCanvas *c1 = new TCanvas("c1", "", 600, 600);
  c1->SaveAs(pn+"[");
  TText *tt = new TText;
  tt->SetTextSize(0.04);
  tt->SetTextAlign(22);
  TString label = "ALICE Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV";
  for (Int_t i=0; i<4; ++i) {
    c1->Clear();
    c1->cd();
    tt->DrawTextNDC(0.5, 0.97, Form("%s - fill %d", frameTitles[i], fill));

    TPad *pad = new TPad("pad", "", 0,0, 1,0.95);
    pad->Draw();
    pad->Divide(3,3);
    TFile::Open(scans[i]);
    gesY->SetLineColor(4);
    gesY->SetMarkerColor(4);

    DrawFrame(pad->cd(1),  0.067, 0.071,  "<X>",  "", 0, maxSep, label); gX->Draw("P"); geX->Draw("PE");
    DrawFrame(pad->cd(2),  0.329, 0.331,  "<Y>",  "", 0, maxSep, label); gY->Draw("P"); geY->Draw("PE");
    DrawFrame(pad->cd(3), -6,     6,      "<Z>",  "", 0, maxSep, label); gZ->Draw("P"); geZ->Draw("PE");

    DrawFrame(pad->cd(4),  0.000, 0.0045, "<XX>", "", 0, maxSep, label); gSX->Draw("P"); geSX->Draw("PE");
    DrawFrame(pad->cd(5),  0.000, 0.0045, "<YY>", "", 0, maxSep, label); gSY->Draw("P"); geSY->Draw("PE");
    DrawFrame(pad->cd(6),  5.000, 6.500,  "<ZZ>", "", 0, maxSep, label); gSZ->Draw("P"); geSZ->Draw("PE");

    DrawFrame(pad->cd(7), -1,    +1,      "<XY>", "", 0, maxSep, label); gCXY->Draw("P"); geCXY->Draw("PE");
    DrawFrame(pad->cd(8),  0,     0.0005, "#color[2]{#mu_{X}'}, #color[4]{#mu_{Y}'}", "", 0, .12, label); gesX->Draw("PE"); gesY->Draw("PE");
    DrawFrame(pad->cd(9),  0.5,   2.5,    "k", "", 0, maxSep, label); gek->Draw("PE");
    c1->SaveAs(pn);
  }
  c1->SaveAs(pn+"]");
}
