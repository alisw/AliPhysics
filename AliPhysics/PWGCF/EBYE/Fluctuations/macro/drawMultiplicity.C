void drawMultiplicity(const char* filename = "AnalysisResults.root") {
  //Draws the multiplicity for each centrality bin
  const Int_t nCentralityBins = 20;

  //_______________________________________________________________//
  //Open the input file
  TFile *f = TFile::Open(filename);
  if(!f->IsOpen()) {
    Printf("File not found!!!");
    break;
  }

  //_______________________________________________________________//
  //Get the TDirectoryFile
  TDirectoryFile *dir = dynamic_cast<TDirectoryFile *>(f->Get("outputFluctuationAnalysis.root"));
  if(!dir) {
    Printf("TDirectoryFile not found!!!");
    break;
  }

  //_______________________________________________________________//
  //Get the TList
  TList *list = dynamic_cast<TList *>(dir->Get("fluctuationsOutput"));
  if(!list) {
    Printf("TList not found!!!");
    break;
  }
  //list->ls();

  TH1F *gHistNMultiplicity[nCentralityBins];
  TH2F *gHistNPlusNMinus[nCentralityBins];
  TString histName;

  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
    histName = "fHistNPlusNMinusCentrality"; histName += iBin;
    gHistNPlusNMinus[iBin-1] = dynamic_cast<TH2F *>(list->FindObject(histName.Data()));
    gHistNPlusNMinus[iBin-1]->SetMarkerColor(iBin);

    histName = "fHistNMultCentrality"; histName += iBin;
    gHistNMultiplicity[iBin-1] = dynamic_cast<TH1F *>(list->FindObject(histName.Data()));
    gHistNMultiplicity[iBin-1]->SetMarkerColor(iBin);
    gHistNMultiplicity[iBin-1]->SetLineColor(iBin);
    gHistNMultiplicity[iBin-1]->Scale(5000./gHistNMultiplicity[iBin-1]->Integral());
  }

  //_______________________________________________________________//
  //Draw the multiplicity distributions
  TH2F *hEmpty = new TH2F("hEmpty","",200,0,2000,500,0,2000);
  hEmpty->SetStats(kFALSE);
  hEmpty->GetXaxis()->SetNdivisions(10);
  hEmpty->GetYaxis()->SetNdivisions(10);
  hEmpty->GetXaxis()->SetTitleOffset(1.4);
  hEmpty->GetYaxis()->SetTitleOffset(1.4);

  TCanvas *c1 = new TCanvas("c1","",0,0,600,500);
  c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->SetLeftMargin(0.15); c1->SetBottomMargin(0.15);
  hEmpty->GetXaxis()->SetTitle("N_{+}");
  hEmpty->GetYaxis()->SetTitle("N_{-}");
  hEmpty->GetXaxis()->SetRangeUser(0,1000);
  hEmpty->GetYaxis()->SetRangeUser(0,1000);
  hEmpty->DrawCopy();
  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++)
    gHistNPlusNMinus[iBin-1]->DrawCopy("same");

  TCanvas *c2 = new TCanvas("c2","",600,0,600,500);
  c2->SetFillColor(10); c2->SetHighLightColor(10);
  c2->SetLeftMargin(0.15); c2->SetBottomMargin(0.15);
  hEmpty->GetXaxis()->SetTitle("N_{tracks}");
  hEmpty->GetYaxis()->SetTitle("Entries");
  hEmpty->GetYaxis()->SetRangeUser(0,2000);
  hEmpty->GetXaxis()->SetRangeUser(0,2000);
  hEmpty->DrawCopy();
  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++)
    gHistNMultiplicity[iBin-1]->DrawCopy("esame");
}
