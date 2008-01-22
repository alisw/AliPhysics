void printStats(const char* fileName = "trackCuts.root")
{
  TFile::Open(fileName);
  
  fTriggerStats = (TH1*) gFile->Get("fTriggerStats");

  // !MB1 + MB1
  Int_t allEvents = fTriggerStats->GetBinContent(1) + fTriggerStats->GetBinContent(2);

  // MB2
  Int_t triggeredEvents = fTriggerStats->GetBinContent(3);

  Printf("Triggered %d out of %d events (%.2f %%)", triggeredEvents, allEvents, 100.0 * triggeredEvents / allEvents);

  fVertex = (TH1*) gFile->Get("fVertex");
  Int_t eventsVertex = fVertex->GetEntries();

  Printf("%d events have a vertex out of %d triggered (%.2f %%)", eventsVertex, triggeredEvents, 100.0 * eventsVertex / triggeredEvents);

  fStatsPrim = (TH1*) gFile->Get("fTrackCutsPrimaries/cut_statistics");

  Int_t tracksPrim = fStatsPrim->GetBinContent(1);
  Int_t tracksPrimAc = tracksPrim - fStatsPrim->GetBinContent(2);

  fStatsSec = (TH1*) gFile->Get("fTrackCutsSecondaries/cut_statistics");

  Int_t tracksSec = fStatsSec->GetBinContent(1);
  Int_t tracksSecAc = tracksSec - fStatsSec->GetBinContent(2);

  fPrimStats = (TH1*) gFile->Get("fPrimStats");

  if (fPrimStats->GetBinContent(5) + fPrimStats->GetBinContent(7) != fStatsPrim->GetBinContent(2))
    Printf("UNEXPECTED: %f != %f", fPrimStats->GetBinContent(5) + fPrimStats->GetBinContent(7), fStatsPrim->GetBinContent(2));

  Float_t notLostPrimaries = 100.0 * fPrimStats->GetBinContent(7) / (fPrimStats->GetBinContent(5) + fPrimStats->GetBinContent(7));
 
  Printf("Accepted %d out of %d primary tracks (%.2f %%)", tracksPrimAc, tracksPrim, 100.0 * tracksPrimAc / tracksPrim);

  Printf("Among the non accepted ones %.2f %% are from primaries that have been found with other tracks", notLostPrimaries);

  Printf("Accepted %d out of %d secondary tracks (%.2f %%)", tracksSecAc, tracksSec, 100.0 * tracksSecAc / tracksSec);

  Printf("Before cuts: %.2f %% of the tracks are secondaries", 100.0 * tracksSec / (tracksPrim + tracksSec));

  Printf("After cuts: %.2f %% of the tracks are secondaries", 100.0 * tracksSecAc / (tracksPrimAc + tracksSecAc));
}

void ComparePlots(const char* fileName1, const char* fileName2, const char* dir = "AliESDtrackCuts/before_cuts")
{
  file1 = TFile::Open(fileName1);
  file2 = TFile::Open(fileName2);

  dir1 = file1->GetDirectory(dir);
  dir2 = file2->GetDirectory(dir);

  keyList = dir1->GetListOfKeys();
  TIter iter(keyList);

  TKey* key = 0;
  while ((key = (TKey*) iter())) {
    TH1* hist1 = (TH1*) dir1->Get(key->GetName());
    TH1* hist2 = (TH1*) dir2->Get(key->GetName());

    if (hist1->Integral() > 0)
      hist1->Scale(1.0 / hist1->Integral());
    if (hist2->Integral())
      hist2->Scale(1.0 / hist2->Integral());

    TString name(key->GetName());

    c = new TCanvas(key->GetName(), key->GetName(), 600, 400);
    hist1->Draw();
    hist2->SetLineColor(2);
    hist2->Draw("SAME");

    Float_t min = 0.9 * TMath::Min(hist1->GetMinimum(), hist2->GetMinimum());

    if (name.Contains("cov") || name.Contains("dXY") || name.Contains("dZ")) {
      min += 1e-6;
      c->cd(1)->SetLogy();
      if (name.Contains("cov")) 
	hist1->GetXaxis()->SetRangeUser(0, 5);
    }

    hist1->GetYaxis()->SetRangeUser(min, 1.1 * TMath::Max(hist1->GetMaximum(), hist2->GetMaximum()));

    c->SaveAs(Form("%s.gif", key->GetName()));
    
    //break;
  }
}
      
