void draw()
{
  //const char* files[] = { "trackCuts_normal.root", "trackCuts_increased.root", "trackCuts_decreased.root" };
  const char* files[] =
    { "Material-normal/trackCuts.root",
      "Material-increased-mcvtx/trackCuts.root",
      "Material-decreased-mcvtx/trackCuts.root" };
  const char* titles[] = { "default geometry", "+ 10% material", "- 10% material" };
  Int_t colors[] = { 1, 2, 4 };
  Int_t markers[] = {24, 25, 26, 27 };

  TCanvas* c = new TCanvas;
  TCanvas* c2 = new TCanvas;

  TLegend* legend = new TLegend(0.6, 0.5, 0.9, 0.9);
  TLegend* legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);

  for (Int_t i=0; i<3; i++) {
    if (!TFile::Open(files[i]))
      return;

    TH1* ptPrim = gFile->Get("fTrackCutsPrimaries/before_cuts/pt");
    TH1* ptPrimCut = gFile->Get("fTrackCutsPrimaries/after_cuts/pt_cut");

    TH1* ptSec = gFile->Get("fTrackCutsSecondaries/before_cuts/pt");
    TH1* ptSecCut = gFile->Get("fTrackCutsSecondaries/after_cuts/pt_cut");

    TH1* vertex = gFile->Get("fVertex");
    Int_t nEvents = vertex->GetEntries();

    ptPrim->Add(ptSec);
    ptPrimCut->Add(ptSecCut);

    ptPrim->Sumw2();
    ptPrimCut->Sumw2();
    ptSec->Sumw2();
    ptSecCut->Sumw2();

    Printf("%s", titles[i]);
    Printf("Total particles: %d", (Int_t) (ptPrim->GetEntries() + ptSec->GetEntries()));
    Printf("Total particles/event: %.2f", (ptPrim->GetEntries() + ptSec->GetEntries()) / nEvents);
    Printf("Primaries/event: %.2f", ptPrim->GetEntries() / nEvents);
    Printf("Secondaries/event: %.2f", ptSec->GetEntries() / nEvents);
    Printf("Primaries > 0.2 GeV/c/event: %.2f", ptPrim->Integral(ptPrim->GetXaxis()->FindBin(0.21), ptPrim->GetNbinsX()) / nEvents);
    Printf("Secondaries > 0.2 GeV/c/event: %.2f", ptSec->Integral(ptSec->GetXaxis()->FindBin(0.21), ptSec->GetNbinsX()) / nEvents);
    Printf("Primaries after cuts > 0.2 GeV/c/event: %.2f +- %.2f", ptPrimCut->Integral(ptPrimCut->GetXaxis()->FindBin(0.21), ptPrimCut->GetNbinsX()) / nEvents, TMath::Sqrt(ptPrimCut->Integral(ptPrimCut->GetXaxis()->FindBin(0.21), ptPrimCut->GetNbinsX())) / nEvents);
    Printf("Secondaries after cuts > 0.2 GeV/c/event: %.2f +- %.2f", ptSecCut->Integral(ptSecCut->GetXaxis()->FindBin(0.21), ptSecCut->GetNbinsX()) / nEvents, TMath::Sqrt(ptSecCut->Integral(ptSecCut->GetXaxis()->FindBin(0.21), ptSecCut->GetNbinsX())) / nEvents);
    Printf("%.2f %% secondaries before cuts", 100.0 * ptSec->GetEntries() / ptPrim->GetEntries());
    Printf("%.2f %% secondaries after cuts", 100.0 * ptSecCut->GetEntries() / ptPrimCut->GetEntries());
    Printf("");

    ptPrim->SetLineColor(colors[i]);
    ptPrimCut->SetLineColor(colors[i]);
    ptSec->SetLineColor(colors[i]);
    ptSecCut->SetLineColor(colors[i]);

    ptPrim->SetMarkerColor(colors[i]);
    ptPrimCut->SetMarkerColor(colors[i]);
    ptSec->SetMarkerColor(colors[i]);
    ptSecCut->SetMarkerColor(colors[i]);

    ptPrim->SetStats(kFALSE);
    ptPrim->SetTitle("");
    ptPrim->GetYaxis()->SetTitle("N");
    ptSec->SetStats(kFALSE);
    ptSec->SetTitle("");

    ptPrim->SetMarkerStyle(markers[0]);
    ptPrimCut->SetMarkerStyle(markers[1]);
    ptSec->SetMarkerStyle(markers[2]);
    ptSecCut->SetMarkerStyle(markers[3]);

    if (i == 0) {
      legend->AddEntry(ptPrim->Clone(), "Primaries");
      legend->AddEntry(ptPrimCut->Clone(), "Primaries after cuts");
      legend->AddEntry(ptSec->Clone(), "Secondaries");
      legend->AddEntry(ptSecCut->Clone(), "Secondaries after cuts");
    }

    ptPrim->GetXaxis()->SetRangeUser(0, 2);
    ptSec->GetXaxis()->SetRangeUser(0, 2);
    //ptPrim->GetYaxis()->SetRangeUser(1e-5, ptSec->GetMaximum() * 1.1);

    c->cd();
    ptPrim->DrawCopy((i > 0) ? "SAME" : "");
    ptPrimCut->DrawCopy("SAME");
    ptSec->DrawCopy("SAME");
    ptSecCut->DrawCopy("SAME");

    ptSec->Divide(ptSec, ptPrim, 1, 1, "B");
    ptSecCut->Divide(ptSecCut, ptPrimCut, 1, 1, "B");

    ptSec->SetMarkerStyle(1);
    ptSecCut->SetMarkerStyle(1);

    ptSec->GetYaxis()->SetRangeUser(0, 1);

    ptSec->GetYaxis()->SetTitle("N_{Secondaries} / N_{All}");

    c2->cd();
    ptSec->DrawCopy((i > 0) ? "SAME" : "");
    ptSecCut->DrawCopy("SAME");

    legend->AddEntry(ptSec, titles[i]);
    legend2->AddEntry(ptSec, titles[i]);
  }

  c->cd();
  legend->Draw();
  c->SaveAs("changedmaterial_absolute.gif");

  c2->cd();
  legend2->Draw();
  c2->SaveAs("changedmaterial_relative.gif");
}

void drawStats(const char* fileName = "trackCuts.root")
{
  TFile::Open(fileName);

  TH1* stat1 = gFile->Get("fTrackCutsPrimaries/cut_statistics");
  TH1* stat2 = gFile->Get("fTrackCutsSecondaries/cut_statistics");
  
  new TCanvas;
  stat1->Draw();
  stat2->SetLineColor(2);
  stat2->Draw("SAME");
}
