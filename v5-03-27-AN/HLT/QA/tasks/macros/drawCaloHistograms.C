drawCaloHistograms(const char* filename="HLT-OFFLINE-PHOS-comparison.root") {
  

  gStyle->SetPalette(1, 0);
  //gStyle->SetOptFit(1111);

  TFile * file = TFile::Open(filename);
  if(!file || file->IsZombie()){
     printf("file %s does not exist or there is an error opening it\n", filename);
     return;
  }
  
  TList * list = dynamic_cast<TList*>(file->Get("phos_histograms"));
  if(!list){
     printf("No list %s contained in your input file\n", list->GetName()); 
     return; 
  }
  
  TH1* hist = list->FindObject("fHistOfflResiduals");
  hist->SetLineColor(kRed);
  hist->DrawNormalized();
  hist = dynamic_cast<TH1*>(list->FindObject("fHistOnlResiduals"));
  hist->SetLineColor(kBlue);
  hist->DrawNormalized("same");


  TCanvas * c2 = new TCanvas("c2", "Matching_histos", 1000, 900);
  c2->Divide(3,2);

  int icd = 1;

  c2->cd(icd++);
  hist = dynamic_cast<TH1*>(list->FindObject("fHistOnlDzNeg"));
  hist->SetFillColor(25);
  hist->SetTitle("Dz of negative tracks");
  hist->DrawNormalized("");
  hist = dynamic_cast<TH1*>(list->FindObject("fHistOfflDzNeg"));
  hist->SetLineColor(kRed);
  hist->DrawNormalized("same");

  c2->cd(icd + 2);
  hist = dynamic_cast<TH1*>(list->FindObject("fHistOnlDzPos"));
  hist->SetFillColor(25);
  hist->SetTitle("dZ of positive tracks");
  hist->DrawNormalized("");
  hist = dynamic_cast<TH1*>(list->FindObject("fHistOfflDzPos"));
  hist->SetLineColor(kRed);
  hist->DrawNormalized("same");
  
  c2->cd(icd++);
  hist = dynamic_cast<TH1*>(list->FindObject("fHistOnlDxyNeg"));
  hist->SetTitle("dXY of negative tracks");
  hist->SetFillColor(25);
  hist->DrawNormalized("");
  hist = dynamic_cast<TH1*>(list->FindObject("fHistOfflDxyNeg"));
  hist->SetLineColor(kRed);
  hist->DrawNormalized("same");

  c2->cd(icd +2);
  hist = dynamic_cast<TH1*>(list->FindObject("fHistOnlDxyPos"));
  hist->SetTitle("dXY of positive tracks");
  hist->SetFillColor(25);
  hist->DrawNormalized("");
  hist = dynamic_cast<TH1*>(list->FindObject("fHistOfflDxyPos"));
  hist->SetLineColor(kRed);
  hist->DrawNormalized("same");

  c2->cd(icd++);
  TH2F * hist2d = dynamic_cast<TH2F*>(list->FindObject("_OFF_fHistdXYdZ"));
  hist2d->SetTitle("dXYdZ offline");
  hist2d->Draw("colz");
  c2->cd(icd +2);
  hist2d = dynamic_cast<TH2F*>(list->FindObject("_ON_fHistdXYdZ"));
  hist2d->SetTitle("dXYdZ hlt");
  hist2d->Draw("colz");

 

  TCanvas * c3 = new TCanvas("c3", "Energy_histos", 1000, 900);
  icd = 1;
  c3->Divide(3,2);

  c3->cd(icd++);
  hist = dynamic_cast<TH1*>(list->FindObject("_ON fHistClusterEnergy"));
  hist->SetTitle("Cluster energy");
  hist->SetAxisRange(0, 1, "X");
  hist->SetFillColor(25);
  hist->DrawNormalized();
  hist = dynamic_cast<TH1*>(list->FindObject("_OFF fHistClusterEnergy"));
  hist->DrawNormalized("same");

  c3->cd(icd+2);
  hist2d = dynamic_cast<TH2F*>(list->FindObject("fHistTotEVsTotE"));
  hist2d->SetTitle("Sum cluster Energies HLT vs offline");
  hist2d->SetAxisRange(0, 10, "X");
  hist2d->SetAxisRange(0, 10, "Y");
  hist2d->Draw("colz");

  c3->cd(icd++);
  hist2d = dynamic_cast<TH2F*>(list->FindObject("_ON fHistClusterEnergyVsNCells"));
  hist2d->SetTitle("Cluster energy vs nCells, HLT"); 
  hist2d->SetAxisRange(0, 5, "X");
  hist2d->SetAxisRange(0, 30, "Y");
  hist2d->Draw("colz");
  c3->cd(icd + 2);
  hist2d = dynamic_cast<TH2F*>(list->FindObject("_OFF fHistClusterEnergyVsNCells"));
  hist2d->SetTitle("Cluster energy vs nCells, Offline"); 
  hist2d->SetAxisRange(0, 5, "X");
  hist2d->SetAxisRange(0, 30, "Y");
  hist2d->Draw("colz");

  c3->cd(icd++);
  hist2d = dynamic_cast<TH2F*>(list->FindObject("_ON fHistClusterEnergyDepositedEtaPhi"));
  hist2d->SetTitle("Energy deposited eta vs phi HLT");
  hist2d->SetAxisRange(-0.2, 0.2, "Y");
  hist2d->SetAxisRange(4, 6, "X");
  hist2d->Draw("colz");
  c3->cd(icd+2);
  hist2d = dynamic_cast<TH2F*>(list->FindObject("_OFF fHistClusterEnergyDepositedEtaPhi"));
  hist2d->SetTitle("Energy deposited eta vs phi Offline");
  hist2d->SetAxisRange(-0.2, 0.2, "Y");
  hist2d->SetAxisRange(4, 6, "X");
  hist2d->Draw("colz");


  // TCanvas * c4 = new TCanvas("c4", "Energy_histos", 1100, 900);
  // icd = 1;
  // c4->Divide(3,2);
  // for( int i = 0; i < 6; i++) {
  //   hist = dynamic_cast<TH1*>(list->FindObject(Form("fHistDArr_%d", i)));
  //   c4->cd(icd++);
  //   hist->Draw();
  //   TFitResultPtr fr = hist->Fit("gaus");
  //   fr->Dump();
  //   //cout << fr->GetParams() << endl;
  // }
	


}
