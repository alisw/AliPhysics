void drawProtonQAResults(const char* filename = "Protons.QA.root") {
  gStyle->SetPalette(1,0);

  const Int_t NQAHISTOSPERLIST = 15;
  
  Double_t gEntriesQA2DList[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t gEntriesQAPrimaryProtonsAcceptedList[NQAHISTOSPERLIST], gEntriesQAPrimaryProtonsRejectedList[NQAHISTOSPERLIST];
  Double_t gEntriesQASecondaryProtonsAcceptedList[NQAHISTOSPERLIST], gEntriesQASecondaryProtonsRejectedList[NQAHISTOSPERLIST];
  Double_t gEntriesQAPrimaryAntiProtonsAcceptedList[NQAHISTOSPERLIST], gEntriesQAPrimaryAntiProtonsRejectedList[NQAHISTOSPERLIST];
  Double_t gEntriesQASecondaryAntiProtonsAcceptedList[NQAHISTOSPERLIST], gEntriesQASecondaryAntiProtonsRejectedList[NQAHISTOSPERLIST];

  for(Int_t i = 0; i < NQAHISTOSPERLIST; i++) {
    gEntriesQAPrimaryProtonsAcceptedList[i] = 0.0;
    gEntriesQAPrimaryProtonsRejectedList[i] = 0.0;
    gEntriesQASecondaryProtonsAcceptedList[i] = 0.0;
    gEntriesQASecondaryProtonsRejectedList[i] = 0.0;
    gEntriesQAPrimaryAntiProtonsAcceptedList[i] = 0.0;
    gEntriesQAPrimaryAntiProtonsRejectedList[i] = 0.0;
    gEntriesQASecondaryAntiProtonsAcceptedList[i] = 0.0;
    gEntriesQASecondaryAntiProtonsRejectedList[i] = 0.0;
  }

  TFile *f = TFile::Open(filename);
  TList *list = (TList *)f->Get("outputList1");

  TList *fQA2DList = (TList *)list->At(0);
  GetQAEntries(fQA2DList,gEntriesQA2DList);

  TList *fQAPrimaryProtonsAcceptedList = (TList *)list->At(1);
  GetQAEntries(fQAPrimaryProtonsAcceptedList,gEntriesQAPrimaryProtonsAcceptedList);

  TList *fQAPrimaryProtonsRejectedList = (TList *)list->At(2);
  GetQAEntries(fQAPrimaryProtonsRejectedList,gEntriesQAPrimaryProtonsRejectedList);

  TList *fQASecondaryProtonsAcceptedList = (TList *)list->At(3);
  GetQAEntries(fQASecondaryProtonsAcceptedList,gEntriesQASecondaryProtonsAcceptedList);

  TList *fQASecondaryProtonsRejectedList = (TList *)list->At(4);
  GetQAEntries(fQASecondaryProtonsRejectedList,gEntriesQASecondaryProtonsRejectedList);

  TList *fQAPrimaryAntiProtonsAcceptedList = (TList *)list->At(5);
  GetQAEntries(fQAPrimaryAntiProtonsAcceptedList,gEntriesQAPrimaryAntiProtonsAcceptedList);

  TList *fQAPrimaryAntiProtonsRejectedList = (TList *)list->At(6);
  GetQAEntries(fQAPrimaryAntiProtonsRejectedList,gEntriesQAPrimaryAntiProtonsRejectedList);

  TList *fQASecondaryAntiProtonsAcceptedList = (TList *)list->At(7);
  GetQAEntries(fQASecondaryAntiProtonsAcceptedList,gEntriesQASecondaryAntiProtonsAcceptedList);

  TList *fQASecondaryAntiProtonsRejectedList = (TList *)list->At(8);
  GetQAEntries(fQASecondaryAntiProtonsRejectedList,gEntriesQASecondaryAntiProtonsRejectedList);

  //_______________________________________________________//
  //Create the histograms
  const Int_t nx = 21;
  char *fCutName[nx] = {"Tracks","",
			"ITS Clusters",
			"#chi^{2}/N_{ITS-Clusters}",
			"TPC Clusters",
			"#chi^{2}/N_{TPC-Clusters}",
			"ExtCov11",
			"ExtCov22",
			"ExtCov33",
			"ExtCov44",
			"ExtCov55",
			"#sigma_{Vertex}",
			"#sigma_{Vertex-TPC}",
			"ITS refit",
			"TPC refit",
			"ESD pid",
			"TPC pid","",
			"N_{Secondaries}/N_{total}","",""};

  TH2F *hEmpty = new TH2F("hEmpty","",nx,0,nx,100,0,100); 
  hEmpty->SetStats(kFALSE); 
  hEmpty->SetMarkerStyle(kFullCircle);
  hEmpty->SetLineColor(2); hEmpty->SetMarkerColor(2);
  hEmpty->GetYaxis()->SetTitle("Influence of the cuts [%]");
  for(Int_t i = 1; i <= nx; i++) 
    hEmpty->GetXaxis()->SetBinLabel(i,fCutName[i-1]);
  hEmpty->GetXaxis()->SetLabelOffset(0.01); 
  hEmpty->GetXaxis()->SetLabelSize(0.045);

  //primary protons
  TH1F *hPrimaryProtons = new TH1F("hPrimaryProtons","",nx,0,nx); 
  hPrimaryProtons->SetStats(kFALSE); 
  hPrimaryProtons->SetMarkerStyle(kFullCircle);
  hPrimaryProtons->SetLineColor(4); hPrimaryProtons->SetLineWidth(2); 
  hPrimaryProtons->SetMarkerColor(4);
  hPrimaryProtons->GetYaxis()->SetTitle("Influence of the cuts [%]");
  for(Int_t i = 1; i <= nx; i++) {
    hPrimaryProtons->SetBinContent(i,-10.);
    hPrimaryProtons->GetXaxis()->SetBinLabel(i,fCutName[i-1]);
  }
  hPrimaryProtons->GetXaxis()->SetLabelOffset(0.01); 
  hPrimaryProtons->GetXaxis()->SetLabelSize(0.045);

  //secondary protons
  TH1F *hSecondaryProtons = new TH1F("hSecondaryProtons","",nx,0,nx); 
  hSecondaryProtons->SetStats(kFALSE); 
  hSecondaryProtons->SetMarkerStyle(22);
  hSecondaryProtons->SetMarkerSize(1.4);
  hSecondaryProtons->SetLineColor(2); hSecondaryProtons->SetLineWidth(2); 
  hSecondaryProtons->SetMarkerColor(2);
  hSecondaryProtons->GetYaxis()->SetTitle("Influence of the cuts [%]");
  for(Int_t i = 1; i <= nx; i++) {
    hSecondaryProtons->SetBinContent(i,-10.);
    hSecondaryProtons->GetXaxis()->SetBinLabel(i,fCutName[i-1]);
  }
  hSecondaryProtons->GetXaxis()->SetLabelOffset(0.01); 
  hSecondaryProtons->GetXaxis()->SetLabelSize(0.045);

  //primary antiprotons
  TH1F *hPrimaryAntiProtons = new TH1F("hPrimaryAntiProtons","",nx,0,nx); 
  hPrimaryAntiProtons->SetStats(kFALSE); 
  hPrimaryAntiProtons->SetMarkerStyle(kFullCircle);
  hPrimaryAntiProtons->SetLineColor(4); hPrimaryAntiProtons->SetLineWidth(2); 
  hPrimaryAntiProtons->SetMarkerColor(4);
  hPrimaryAntiProtons->GetYaxis()->SetTitle("Influence of the cuts [%]");
  for(Int_t i = 1; i <= nx; i++) {
    hPrimaryAntiProtons->SetBinContent(i,-10.);
    hPrimaryAntiProtons->GetXaxis()->SetBinLabel(i,fCutName[i-1]);
  }
  hPrimaryAntiProtons->GetXaxis()->SetLabelOffset(0.01); 
  hPrimaryAntiProtons->GetXaxis()->SetLabelSize(0.045);

  //secondary antiprotons
  TH1F *hSecondaryAntiProtons = new TH1F("hSecondaryAntiProtons","",nx,0,nx); 
  hSecondaryAntiProtons->SetStats(kFALSE); 
  hSecondaryAntiProtons->SetMarkerStyle(22);
  hSecondaryAntiProtons->SetMarkerSize(1.4);
  hSecondaryAntiProtons->SetLineColor(2); hSecondaryAntiProtons->SetLineWidth(2); 
  hSecondaryAntiProtons->SetMarkerColor(2);
  hSecondaryAntiProtons->GetYaxis()->SetTitle("Influence of the cuts [%]");
  for(Int_t i = 1; i <= nx; i++) {
    hSecondaryAntiProtons->SetBinContent(i,-10.);
    hSecondaryAntiProtons->GetXaxis()->SetBinLabel(i,fCutName[i-1]);
  }
  hSecondaryAntiProtons->GetXaxis()->SetLabelOffset(0.01); 
  hSecondaryAntiProtons->GetXaxis()->SetLabelSize(0.045);
  //_______________________________________________________//

  //1D for primary protons
  cout<<"_____________________________________________________"<<endl;
  cout<<"_______________PRIMARY PROTONS_______________________"<<endl;
  hPrimaryProtons->SetBinContent(1,GetPercentage(gEntriesQA2DList[0],gEntriesQA2DList[1]));

  for(Int_t i = 2; i < 17; i++) 
    hPrimaryProtons->SetBinContent(i+1,GetPercentage(gEntriesQAPrimaryProtonsAcceptedList[i-2],
						   gEntriesQAPrimaryProtonsRejectedList[i-2]));
  
  //1D for secondary protons
  cout<<"_____________________________________________________"<<endl;
  cout<<"_______________SECONDARY PROTONS_____________________"<<endl;
  hSecondaryProtons->SetBinContent(1,GetPercentage(gEntriesQA2DList[2],gEntriesQA2DList[3]));
  hSecondaryProtons->SetBinContent(19,GetPercentage(gEntriesQA2DList[0],gEntriesQA2DList[2]));

  for(Int_t i = 2; i < 17; i++) 
    hSecondaryProtons->SetBinContent(i+1,GetPercentage(gEntriesQASecondaryProtonsAcceptedList[i-2],
						     gEntriesQASecondaryProtonsRejectedList[i-2]));

  //1D for primary antiprotons
  cout<<"_________________________________________________________"<<endl;
  cout<<"_______________PRIMARY ANTIPROTONS_______________________"<<endl;
  hPrimaryAntiProtons->SetBinContent(1,GetPercentage(gEntriesQA2DList[4],gEntriesQA2DList[5]));

  for(Int_t i = 2; i < 17; i++) 
    hPrimaryAntiProtons->SetBinContent(i+1,GetPercentage(gEntriesQAPrimaryAntiProtonsAcceptedList[i-2],
						       gEntriesQAPrimaryAntiProtonsRejectedList[i-2]));

  //1D for secondary antiprotons
  cout<<"_________________________________________________________"<<endl;
  cout<<"_______________SECONDARY ANTIPROTONS_____________________"<<endl;
  hSecondaryAntiProtons->SetBinContent(1,GetPercentage(gEntriesQA2DList[6],gEntriesQA2DList[7]));
  hSecondaryAntiProtons->SetBinContent(19,GetPercentage(gEntriesQA2DList[4],gEntriesQA2DList[6]));

  for(Int_t i = 2; i < 17; i++) 
    hSecondaryAntiProtons->SetBinContent(i+1,GetPercentage(gEntriesQASecondaryAntiProtonsAcceptedList[i-2],
							 gEntriesQASecondaryAntiProtonsRejectedList[i-2]));

  TLatex *t1 = new TLatex();
  t1->SetTextSize(0.04);
  //_________________________________________________________//
  TCanvas *c1 = new TCanvas("c1","Cut Influence - Protons",0,0,450,450);
  c1->SetFillColor(10); c1->GetFrame()->SetFillColor(10);
  c1->SetHighLightColor(10); c1->SetBottomMargin(0.15);
  c1->SetGridx(); c1->SetGridy();

  for(Int_t i = 1; i <= nx; i++) {
    hPrimaryProtons->SetBinError(i,1.0);
    hSecondaryProtons->SetBinError(i,1.0);
  }
  hEmpty->Draw();
  hPrimaryProtons->Draw("EHISTSAME");
  hSecondaryProtons->Draw("EHISTSAME");
  DrawMarker(14.5, 80, 20, 1.2, 4);
  t1->DrawLatex(15,78,"Primary p");
  DrawMarker(14.5, 70, 22, 1.2, 2);
  t1->DrawLatex(15,68,"Secondary p");

  c1->SaveAs("CutInfluence-Protons,gif");

  //_________________________________________________________//
  TCanvas *c2 = new TCanvas("c2","Cut Influence - AntiProtons",450,0,450,450);
  c2->SetFillColor(10); c2->GetFrame()->SetFillColor(10);
  c2->SetHighLightColor(10); c2->SetBottomMargin(0.15);
  c2->SetGridx(); c2->SetGridy();

  for(Int_t i = 1; i <= nx; i++) {
    hPrimaryAntiProtons->SetBinError(i,1.0);
    hSecondaryAntiProtons->SetBinError(i,1.0);
  }
  hEmpty->Draw();
  hPrimaryAntiProtons->Draw("EHISTSAME");
  hSecondaryAntiProtons->Draw("EHISTSAME");
  DrawMarker(14.5, 80, 20, 1.2, 4);
  t1->DrawLatex(15,78,"Primary #bar{p}");
  DrawMarker(14.5, 70, 22, 1.2, 2);
  t1->DrawLatex(15,68,"Secondary #bar{p}");

  c2->SaveAs("CutInfluence-AntiProtons,gif");
}

//________________________________________________//
void GetQAEntries(TList *inputList, Double_t *entries) {
  //loops over the list entries
  //extracts the entries for each histogram
  cout<<"Extracting the entries for the histograms in the list: "<<
    inputList->GetName()<<"..."<<endl;

  for(Int_t i = 0; i < inputList->GetEntries(); i++) {
    TH1F *gHist = (TH1F *)inputList->At(i);
    entries[i] = gHist->GetEntries();
    cout<<"Histogram: "<<gHist->GetName()<<
      " - Entries: "<<entries[i]<<endl;
    gHist = 0;
  }
}

//________________________________________________//
Double_t GetPercentage(Double_t nPassEntries,
		       Double_t nRejectEntries) {
  //returns the percentage of tracks that were rejected by a cut
  Int_t nTotalEntries = nPassEntries + nRejectEntries;

  if(nTotalEntries == 0)
    return -1;

  return 100.*nRejectEntries/nTotalEntries;
}

//________________________________________________//
void DrawMarker(Double_t x, Double_t y, Int_t style, Double_t size, Int_t color) {
  TMarker *m = new TMarker(x,y,style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(color);
  m->Draw();
}
