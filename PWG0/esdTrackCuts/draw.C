void draw()
{
  const char* files[] = { "trackCuts_normal.root", "trackCuts_increased.root", "trackCuts_decreased.root" };
  const char* titles[] = { "default geometry", "+ 10% material", "- 10% material" };
  Int_t colors[] = { 1, 2, 4 };

  TCanvas* c = new TCanvas;

  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

  for (Int_t i=0; i<3; i++) {
    TFile::Open(files[i]);
    
    TH1* ptPrim = gFile->Get("fTrackCutsPrimaries/before_cuts/pt");
    TH1* ptPrimCut = gFile->Get("fTrackCutsPrimaries/after_cuts/pt_cut");

    TH1* ptSec = gFile->Get("fTrackCutsSecondaries/before_cuts/pt");
    TH1* ptSecCut = gFile->Get("fTrackCutsSecondaries/after_cuts/pt_cut");

    ptPrim->Add(ptSec);
    ptPrimCut->Add(ptSecCut);

    ptPrim->Sumw2();
    ptPrimCut->Sumw2();
    ptSec->Sumw2();
    ptSecCut->Sumw2();
    
    Printf("%s", titles[i]);
    Printf("%.2f %% secondaries before cuts", 100.0 * ptSec->GetEntries() / ptPrim->GetEntries());
    Printf("%.2f %% secondaries after cuts", 100.0 * ptSecCut->GetEntries() / ptPrimCut->GetEntries());
    Printf("");

    ptSec->Divide(ptSec, ptPrim, 1, 1, "B");
    ptSecCut->Divide(ptSecCut, ptPrimCut, 1, 1, "B");

    ptSec->SetLineColor(colors[i]);
    ptSecCut->SetLineColor(colors[i]);
    ptSec->SetStats(kFALSE);

    ptSec->GetXaxis()->SetRangeUser(0, 2);
    ptSec->GetYaxis()->SetRangeUser(0, 1);

    ptSec->SetTitle("");
    ptSec->GetYaxis()->SetTitle("N_{Secondaries} / N_{All}");
    
    ptSec->DrawCopy((i > 0) ? "SAME" : "");
    ptSecCut->DrawCopy("SAME");

    legend->AddEntry(ptSec, titles[i]);
  }

  legend->Draw();

  c->SaveAs("secondaries_changedmaterial.gif");
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
