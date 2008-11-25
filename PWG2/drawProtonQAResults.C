void drawProtonQAResults(const char* filename1 = "Protons.QA.root",
			 const char* filename2 = "Protons.MC.QA.root",
			 const char* filename3 = "Protons.Efficiency.root") {
  //Macro to visualize the results of the proton QA task
  //TCanvas objects: 15
  gStyle->SetPalette(1,0);
 
  TFile *fQA = TFile::Open(filename1);
  TList *listGlobalQA = (TList *)fQA->Get("globalQAList");
  drawCutStatistics(listGlobalQA);
 
  TFile *fMC = TFile::Open(filename2);
  TList *listPDG = (TList *)fMC->Get("pdgCodeList");  
  TList *listMCProcesses = (TList *)fMC->Get("mcProcessList");  
  drawMCQA(listPDG,listMCProcesses);
  
  TFile *fEfficiency = TFile::Open(filename3);
  TList *listEfficiency = (TList *)fEfficiency->Get("efficiencyList");
  drawEfficiency(listEfficiency);

  fQA->Close();
  fMC->Close();
  fEfficiency->Close();
}

//________________________________________//
void drawCutStatistics(TList *list) {
  //Function to display the statistics from the cuts
  //The histogram shows the influence of each cut on the primary
  //and secondary (anti)protons
  const Int_t NQAHISTOSPERLIST = 26;
  
  Double_t gEntriesQA2DList[12] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
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

  TList *fQA2DList = (TList *)list->At(0);
  GetQAEntries(fQA2DList,gEntriesQA2DList);
  TH3F *gHistYPtPDGProtonsPass = (TH3F *)fQA2DList->At(10);
  TH3F *gHistYPtPDGAntiProtonsPass = (TH3F *)fQA2DList->At(11);

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
  const Int_t nx = 26;
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
			"DCA_{xy}",
			"DCA_{xy}(TPC)",
			"DCA_{z}",
			"DCA_{z}(TPC)",
			"#chi^{2}(vertex)",
			"ITS refit",
			"TPC refit",
			"ESD pid",
			"TPC pid","",
			"N_{Secondaries}/N_{total}","",""};
  char *fCutITSName[6] = {"SPD_{1}","SPD_{2}",
			  "SDD_{1}","SDD_{2}",
			  "SSD_{1}","SSD_{2}"};

  //cut influence
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
  //cout<<"_____________________________________________________"<<endl;
  //cout<<"_______________PRIMARY PROTONS_______________________"<<endl;
  hPrimaryProtons->SetBinContent(1,GetPercentage(gEntriesQA2DList[0],gEntriesQA2DList[1]));

  for(Int_t i = 2; i < NQAHISTOSPERLIST-4; i++) 
    hPrimaryProtons->SetBinContent(i+1,GetPercentage(gEntriesQAPrimaryProtonsAcceptedList[i-2],
						     gEntriesQAPrimaryProtonsRejectedList[i-2]));
  
  //1D for secondary protons
  //cout<<"_____________________________________________________"<<endl;
  //cout<<"_______________SECONDARY PROTONS_____________________"<<endl;
  hSecondaryProtons->SetBinContent(1,GetPercentage(gEntriesQA2DList[2],gEntriesQA2DList[3]));
  
  for(Int_t i = 2; i < NQAHISTOSPERLIST-4; i++) 
    hSecondaryProtons->SetBinContent(i+1,GetPercentage(gEntriesQASecondaryProtonsAcceptedList[i-2],
						       gEntriesQASecondaryProtonsRejectedList[i-2]));
  hSecondaryProtons->SetBinContent(24,GetPercentage(gEntriesQA2DList[0],gEntriesQA2DList[2]));

  //1D for primary antiprotons
  //cout<<"_________________________________________________________"<<endl;
  //cout<<"_______________PRIMARY ANTIPROTONS_______________________"<<endl;
  hPrimaryAntiProtons->SetBinContent(1,GetPercentage(gEntriesQA2DList[4],gEntriesQA2DList[5]));
  
  for(Int_t i = 2; i < NQAHISTOSPERLIST-4; i++) 
    hPrimaryAntiProtons->SetBinContent(i+1,GetPercentage(gEntriesQAPrimaryAntiProtonsAcceptedList[i-2],
							 gEntriesQAPrimaryAntiProtonsRejectedList[i-2]));
  
  //1D for secondary antiprotons
  //cout<<"_________________________________________________________"<<endl;
  //cout<<"_______________SECONDARY ANTIPROTONS_____________________"<<endl;
  hSecondaryAntiProtons->SetBinContent(1,GetPercentage(gEntriesQA2DList[6],gEntriesQA2DList[7]));

  for(Int_t i = 2; i < NQAHISTOSPERLIST-4; i++) 
    hSecondaryAntiProtons->SetBinContent(i+1,GetPercentage(gEntriesQASecondaryAntiProtonsAcceptedList[i-2],
							   gEntriesQASecondaryAntiProtonsRejectedList[i-2]));
  hSecondaryAntiProtons->SetBinContent(24,GetPercentage(gEntriesQA2DList[4],gEntriesQA2DList[6]));

  TLatex *t1 = new TLatex();
  t1->SetTextSize(0.04);
  //_________________________________________________________//
  TCanvas *c1 = new TCanvas("c1","Cut Influence - Protons",0,0,700,400);
  c1->SetFillColor(10); c1->GetFrame()->SetFillColor(10);
  c1->SetHighLightColor(10); c1->SetBottomMargin(0.15);
  c1->SetGridx(); c1->SetGridy();

  for(Int_t i = 1; i <= nx; i++) {
    hPrimaryProtons->SetBinError(i,1.0);
    hSecondaryProtons->SetBinError(i,1.0);
  }
  hEmpty->DrawCopy();
  hPrimaryProtons->DrawCopy("EHISTSAME");
  hSecondaryProtons->DrawCopy("EHISTSAME");
  DrawMarker(20.5, 90, 20, 1.2, 4);
  t1->DrawLatex(21,88,"Primary p");
  DrawMarker(20.5, 80, 22, 1.2, 2);
  t1->DrawLatex(21,78,"Secondary p");

  c1->SaveAs("CutInfluence-Protons.gif");

  //_________________________________________________________//
  TCanvas *c2 = new TCanvas("c2","Cut Influence - AntiProtons",50,50,700,400);
  c2->SetFillColor(10); c2->GetFrame()->SetFillColor(10);
  c2->SetHighLightColor(10); c2->SetBottomMargin(0.15);
  c2->SetGridx(); c2->SetGridy();

  for(Int_t i = 1; i <= nx; i++) {
    hPrimaryAntiProtons->SetBinError(i,1.0);
    hSecondaryAntiProtons->SetBinError(i,1.0);
  }
  hEmpty->DrawCopy();
  hPrimaryAntiProtons->DrawCopy("EHISTSAME");
  hSecondaryAntiProtons->DrawCopy("EHISTSAME");
  DrawMarker(20.5, 90, 20, 1.2, 4);
  t1->DrawLatex(21,88,"Primary #bar{p}");
  DrawMarker(20.5, 80, 22, 1.2, 2);
  t1->DrawLatex(21,78,"Secondary #bar{p}");

  c2->SaveAs("CutInfluence-AntiProtons.gif");

  //_________________________________________________________//
  //ITS layers influence
  TH2F *hEmptyITS = new TH2F("hEmptyITS","",10,0,10,100,0,100); 
  hEmptyITS->SetStats(kFALSE); 
  hEmptyITS->SetMarkerStyle(kFullCircle);
  hEmptyITS->SetLineColor(2); hEmptyITS->SetMarkerColor(2);
  hEmptyITS->GetYaxis()->SetTitle("Influence of the ITS layers [%]");
  hEmptyITS->GetYaxis()->SetTitleOffset(1.3);
  for(Int_t i = 1; i <= 6; i++) 
    hEmptyITS->GetXaxis()->SetBinLabel(i,fCutITSName[i-1]);
  hEmptyITS->GetXaxis()->SetLabelOffset(0.01); 
  hEmptyITS->GetXaxis()->SetLabelSize(0.045);

  //_________________________________________________________//
  //primary protons
  TH1F *hPrimaryProtonsITS = new TH1F("hPrimaryProtonsITS","",7,0,7); 
  hPrimaryProtonsITS->SetStats(kFALSE); 
  hPrimaryProtonsITS->SetMarkerStyle(kFullCircle);
  hPrimaryProtonsITS->SetLineColor(4); hPrimaryProtonsITS->SetLineWidth(2); 
  hPrimaryProtonsITS->SetMarkerColor(4);
  hPrimaryProtonsITS->GetYaxis()->SetTitle("Influence of the cuts [%]");
  for(Int_t i = 1; i <= 6; i++) {
    hPrimaryProtonsITS->SetBinContent(i,-10.);
    hPrimaryProtonsITS->GetXaxis()->SetBinLabel(i,fCutITSName[i-1]);
  }
  hPrimaryProtonsITS->GetXaxis()->SetLabelOffset(0.01); 
  hPrimaryProtonsITS->GetXaxis()->SetLabelSize(0.045);

  //secondary protons
  TH1F *hSecondaryProtonsITS = new TH1F("hSecondaryProtonsITS","",7,0,7); 
  hSecondaryProtonsITS->SetStats(kFALSE); 
  hSecondaryProtonsITS->SetMarkerStyle(22);
  hSecondaryProtonsITS->SetMarkerSize(1.4);
  hSecondaryProtonsITS->SetLineColor(2); hSecondaryProtonsITS->SetLineWidth(2); 
  hSecondaryProtonsITS->SetMarkerColor(2);
  hSecondaryProtonsITS->GetYaxis()->SetTitle("Influence of the cuts [%]");
  for(Int_t i = 1; i <= 6; i++) {
    hSecondaryProtonsITS->SetBinContent(i,-10.);
    hSecondaryProtonsITS->GetXaxis()->SetBinLabel(i,fCutITSName[i-1]);
  }
  hSecondaryProtonsITS->GetXaxis()->SetLabelOffset(0.01); 
  hSecondaryProtonsITS->GetXaxis()->SetLabelSize(0.045);

  //primary antiprotons
  TH1F *hPrimaryAntiProtonsITS = new TH1F("hPrimaryAntiProtonsITS","",7,0,7); 
  hPrimaryAntiProtonsITS->SetStats(kFALSE); 
  hPrimaryAntiProtonsITS->SetMarkerStyle(kFullCircle);
  hPrimaryAntiProtonsITS->SetLineColor(4); hPrimaryAntiProtonsITS->SetLineWidth(2); 
  hPrimaryAntiProtonsITS->SetMarkerColor(4);
  hPrimaryAntiProtonsITS->GetYaxis()->SetTitle("Influence of the cuts [%]");
  for(Int_t i = 1; i <= 6; i++) {
    hPrimaryAntiProtonsITS->SetBinContent(i,-10.);
    hPrimaryAntiProtonsITS->GetXaxis()->SetBinLabel(i,fCutITSName[i-1]);
  }
  hPrimaryAntiProtonsITS->GetXaxis()->SetLabelOffset(0.01); 
  hPrimaryAntiProtonsITS->GetXaxis()->SetLabelSize(0.045);

  //secondary antiprotons
  TH1F *hSecondaryAntiProtonsITS = new TH1F("hSecondaryAntiProtonsITS","",7,0,7); 
  hSecondaryAntiProtonsITS->SetStats(kFALSE); 
  hSecondaryAntiProtonsITS->SetMarkerStyle(22);
  hSecondaryAntiProtonsITS->SetMarkerSize(1.4);
  hSecondaryAntiProtonsITS->SetLineColor(2); hSecondaryAntiProtonsITS->SetLineWidth(2); 
  hSecondaryAntiProtonsITS->SetMarkerColor(2);
  hSecondaryAntiProtonsITS->GetYaxis()->SetTitle("Influence of the cuts [%]");
  for(Int_t i = 1; i <= 6; i++) {
    hSecondaryAntiProtonsITS->SetBinContent(i,-10.);
    hSecondaryAntiProtonsITS->GetXaxis()->SetBinLabel(i,fCutITSName[i-1]);
  }
  hSecondaryAntiProtonsITS->GetXaxis()->SetLabelOffset(0.01); 
  hSecondaryAntiProtonsITS->GetXaxis()->SetLabelSize(0.045);

  //_______________________________________________________//
  TCanvas *c9 = new TCanvas("c9","ITS cluster map - (anti)protons",
			    100,100,700,400);
  c9->SetFillColor(10); c9->GetFrame()->SetFillColor(10);
  c9->SetHighLightColor(10); c9->Divide(2,1);
  for(Int_t i = 1; i <= 6; i++) {
    hPrimaryProtonsITS->SetBinError(i,1.0);
    hSecondaryProtonsITS->SetBinError(i,1.0);
    hPrimaryAntiProtonsITS->SetBinError(i,1.0);
    hSecondaryAntiProtonsITS->SetBinError(i,1.0);
  }
  c9->cd(1)->SetBottomMargin(0.15);
  c9->cd(1)->SetLeftMargin(0.15);
  c9->cd(1)->SetGridx(); c9->cd(1)->SetGridy();
  hEmptyITS->SetTitle("Protons");
  hEmptyITS->DrawCopy();

  for(Int_t i = 1; i < 7; i++) {
    hPrimaryProtonsITS->SetBinContent(i,GetPercentage(gEntriesQAPrimaryProtonsAcceptedList[19+i],
						      gEntriesQAPrimaryProtonsRejectedList[19+i]));
    hSecondaryProtonsITS->SetBinContent(i,GetPercentage(gEntriesQASecondaryProtonsAcceptedList[19+i],
							gEntriesQASecondaryProtonsRejectedList[19+i]));
  }
  hPrimaryProtonsITS->DrawCopy("EHISTSAME");
  hSecondaryProtonsITS->DrawCopy("EHISTSAME");
  DrawMarker(6.5, 90, 20, 1.2, 4);
  t1->DrawLatex(7,88,"Primary p");
  DrawMarker(6.5, 80, 22, 1.2, 2);
  t1->DrawLatex(7,78,"Secondary p");

  c9->cd(2)->SetBottomMargin(0.15);
  c9->cd(2)->SetLeftMargin(0.15);
  c9->cd(2)->SetGridx(); c9->cd(2)->SetGridy();
  hEmptyITS->SetTitle("Antiprotons");
  hEmptyITS->DrawCopy();
  for(Int_t i = 1; i < 7; i++) {
    hPrimaryAntiProtonsITS->SetBinContent(i,GetPercentage(gEntriesQAPrimaryAntiProtonsAcceptedList[19+i],
							  gEntriesQAPrimaryAntiProtonsRejectedList[19+i]));
    hSecondaryAntiProtonsITS->SetBinContent(i,GetPercentage(gEntriesQASecondaryAntiProtonsAcceptedList[19+i],
							    gEntriesQASecondaryAntiProtonsRejectedList[19+i]));
  }
  hPrimaryAntiProtonsITS->DrawCopy("EHISTSAME");
  hSecondaryAntiProtonsITS->DrawCopy("EHISTSAME");
  DrawMarker(6.5, 90, 20, 1.2, 4);
  t1->DrawLatex(7,88,"Primary #bar{p}");
  DrawMarker(6.5, 80, 22, 1.2, 2);
  t1->DrawLatex(7,78,"Secondary #bar{p}");

  c9->SaveAs("CutInfluence-ITS.gif");
  
  //Efficiency - Contamination plots
  DrawContamination(fQA2DList);
  DrawCutEfficiency(fQA2DList);
  DrawComposition(gHistYPtPDGProtonsPass,gHistYPtPDGAntiProtonsPass);
}

//________________________________________//
void DrawComposition(TH3F *gHistYPtPDGProtons,
		     TH3F *gHistYPtPDGAntiProtons) {
  //Function to display the composition of secondary (anti)protons
  //that survive the quality criteria
  Double_t nParticleCompositionProtonY[100], nParticleCompositionProtonPt[100];
  Double_t nParticleCompositionAntiProtonY[100], nParticleCompositionAntiProtonPt[100];
  Double_t gY[100], gPt[100];
  for(Int_t iBins = 0; iBins < 100; iBins++) {
    nParticleCompositionProtonY[iBins] = 0;
    nParticleCompositionProtonPt[iBins] = 0;
    nParticleCompositionAntiProtonY[iBins] = 0;
    nParticleCompositionAntiProtonPt[iBins] = 0;
    gY[iBins] = 0;
    gPt[iBins] = 0;
  }
  
  TGraph *gParticleProtonY[14];
  TGraph *gParticleProtonPt[14];
  TGraph *gParticleAntiProtonY[14];
  TGraph *gParticleAntiProtonPt[14];
  for(Int_t iParticle = 0; iParticle < 14; iParticle++) {
    GetComposition(iParticle,
		   gHistYPtPDGProtons,
		   nParticleCompositionProtonY,
		   gY, nParticleCompositionProtonPt, gPt);
    gParticleProtonY[iParticle] = new TGraph(gHistYPtPDGProtons->GetNbinsX(),
				       gY,nParticleCompositionProtonY);
    gParticleProtonY[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleProtonY[iParticle]->SetMarkerSize(1.2);

    gParticleProtonPt[iParticle] = new TGraph(gHistYPtPDGProtons->GetNbinsY(),
					gPt,nParticleCompositionProtonPt);
    gParticleProtonPt[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleProtonPt[iParticle]->SetMarkerSize(1.2);

    GetComposition(iParticle,
		   gHistYPtPDGAntiProtons,
		   nParticleCompositionAntiProtonY,
		   gY, nParticleCompositionAntiProtonPt, gPt);
    gParticleAntiProtonY[iParticle] = new TGraph(gHistYPtPDGAntiProtons->GetNbinsX(),
				       gY,nParticleCompositionAntiProtonY);
    gParticleAntiProtonY[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleAntiProtonY[iParticle]->SetMarkerSize(1.2);

    gParticleAntiProtonPt[iParticle] = new TGraph(gHistYPtPDGAntiProtons->GetNbinsY(),
					gPt,nParticleCompositionAntiProtonPt);
    gParticleAntiProtonPt[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleAntiProtonPt[iParticle]->SetMarkerSize(1.2);
  }

  //_________________________________________________________//
  char *fParticleName[14] = {"Primary","K_{L}","#pi","K_{S}","K",
			     "n","p","#Sigma^{-}","#Lambda","#Sigma^{+}",
			     "#Xi^{-}","#Xi^{0}","#Omega^{-}"};
  TLatex *t1 = new TLatex();
  t1->SetTextSize(0.04);

  TH2F *hEmptyY = new TH2F("hEmptyY","",100,-1.2,1.2,100,0,120); 
  hEmptyY->SetStats(kFALSE); 
  hEmptyY->GetYaxis()->SetTitle("Particle composition [%]");
  hEmptyY->GetXaxis()->SetTitle("y");

  TCanvas *c12 = new TCanvas("c12",
			     "Composition of accepted secondaries vs y",
			    550,550,700,400);
  c12->SetFillColor(10); c12->GetFrame()->SetFillColor(10); c12->Divide(2,1);
  c12->SetHighLightColor(10); c12->cd(1)->SetBottomMargin(0.15);
  c12->cd(1)->SetGridx(); c12->cd(1)->SetGridy();
  hEmptyY->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    //if((iParticle == 0)||(iParticle == 2)||(iParticle == 5)||(iParticle == 6)||(iParticle == 8))
    if((iParticle == 0)||(iParticle == 2)||(iParticle == 6)||(iParticle == 8))
      gParticleProtonY[iParticle]->Draw("P");
    if(iParticle < 5) {
      DrawMarker(-1.1, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(-1.0,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(0.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*(iParticle-5),fParticleName[iParticle]);
    }
  }

  c12->SetHighLightColor(10); c12->cd(2)->SetBottomMargin(0.15);
  c12->cd(2)->SetGridx(); c12->cd(2)->SetGridy();
  hEmptyY->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    if((iParticle == 0)||(iParticle == 6)||(iParticle == 8))
      gParticleAntiProtonY[iParticle]->Draw("P");
    if(iParticle < 5) {
      DrawMarker(-1.1, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(-1.0,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(0.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*(iParticle-5),fParticleName[iParticle]);
    }
  }
  c12->SaveAs("SurvivedSecondaries-Composition-Rapidity.gif");

  TH2F *hEmptyPt = new TH2F("hEmptyPt","",100,0.0,4.0,100,0,120); 
  hEmptyPt->SetStats(kFALSE); 
  hEmptyPt->GetYaxis()->SetTitle("Particle composition [%]");
  hEmptyPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");

  TCanvas *c13 = new TCanvas("c13",
			    "Composition of accepted secondaries vs pT",
			    600,600,700,400);
  c13->SetFillColor(10); c13->GetFrame()->SetFillColor(10); c13->Divide(2,1);
  c13->SetHighLightColor(10); c13->cd(1)->SetBottomMargin(0.15);
  c13->cd(1)->SetGridx(); c13->cd(1)->SetGridy();
  hEmptyPt->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    if(iParticle < 5) {
      DrawMarker(0.2, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(2.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(2.3,113-5*(iParticle-5),fParticleName[iParticle]);
    }
    if((iParticle == 0)||(iParticle == 2)||(iParticle == 6)||(iParticle == 8))
      gParticleProtonPt[iParticle]->Draw("P");
  }

  c13->SetHighLightColor(10); c13->cd(2)->SetBottomMargin(0.15);
  c13->cd(2)->SetGridx(); c13->cd(2)->SetGridy();
  hEmptyPt->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    if(iParticle < 5) {
      DrawMarker(0.2, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(2.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(2.3,113-5*(iParticle-5),fParticleName[iParticle]);
    }
    if((iParticle == 0)||(iParticle == 6)||(iParticle == 8))
      gParticleAntiProtonPt[iParticle]->Draw("P");
  }
  c13->SaveAs("SurvivedSecondaries-Composition-Pt.gif");
}

//________________________________________________//
void DrawContamination(TList *inputList) {
  //loops over the list entries and
  //draws the rapidity and pT dependence
  //of the percentage of primary and secondary
  //protons and antiprotons after the track cuts
  cout<<"Extracting the entries for the histograms in the list: "<<
    inputList->GetName()<<"..."<<endl;

  TLatex *t1 = new TLatex();
  t1->SetTextSize(0.04);

  TH2F *hPrimaryProtons = (TH2F *)inputList->At(0);
  TH2F *hSecondaryProtons = (TH2F *)inputList->At(2);
  TH2F *hPrimaryAntiProtons = (TH2F *)inputList->At(4);
  TH2F *hSecondaryAntiProtons = (TH2F *)inputList->At(6);

  //rapidity dependence
  //Protons
  TH1D *gYPrimaryProtons = (TH1D *)hPrimaryProtons->ProjectionX("gYPrimaryProtons",0,hPrimaryProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYSecondaryProtons = (TH1D *)hSecondaryProtons->ProjectionX("gYSecondaryProtons",0,hSecondaryProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYTotalProtons = (TH1D *)hPrimaryProtons->ProjectionX("gYTotalProtons",0,hPrimaryProtons->GetXaxis()->GetNbins(),"e");
  gYTotalProtons->Add(gYSecondaryProtons);

  TH1D *gYPrimaryProtonsPercentage = new TH1D("gYPrimaryProtonsPercentage",
					      "",
					      hPrimaryProtons->GetXaxis()->GetNbins(),
					      hPrimaryProtons->GetXaxis()->GetXmin(),
					      hPrimaryProtons->GetXaxis()->GetXmax());					      
  gYPrimaryProtonsPercentage->Divide(gYPrimaryProtons,
				     gYTotalProtons,1.,1.0);
  SetError(gYPrimaryProtonsPercentage,gYTotalProtons);
  gYPrimaryProtonsPercentage->Scale(100.);
  gYPrimaryProtonsPercentage->SetMarkerStyle(kFullCircle);

  TH1D *gYSecondaryProtonsPercentage = new TH1D("gYSecondaryProtonsPercentage",
						"",
						hSecondaryProtons->GetXaxis()->GetNbins(),
						hSecondaryProtons->GetXaxis()->GetXmin(),
						hSecondaryProtons->GetXaxis()->GetXmax());					      
  gYSecondaryProtonsPercentage->Divide(gYSecondaryProtons,
				       gYTotalProtons,1.,1.0);
  SetError(gYSecondaryProtonsPercentage,gYTotalProtons);
  gYSecondaryProtonsPercentage->Scale(100.);
  gYSecondaryProtonsPercentage->SetMarkerStyle(kOpenCircle);

  //Antiprotons
  TH1D *gYPrimaryAntiProtons = (TH1D *)hPrimaryAntiProtons->ProjectionX("gYPrimaryAntiProtons",0,hPrimaryAntiProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYSecondaryAntiProtons = (TH1D *)hSecondaryAntiProtons->ProjectionX("gYSecondaryAntiProtons",0,hSecondaryAntiProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYTotalAntiProtons = (TH1D *)hPrimaryAntiProtons->ProjectionX("gYTotalAntiProtons",0,hPrimaryAntiProtons->GetXaxis()->GetNbins(),"e");
  gYTotalAntiProtons->Add(gYSecondaryAntiProtons);
  
  TH1D *gYPrimaryAntiProtonsPercentage = new TH1D("gYPrimaryAntiProtonsPercentage",
						  "",
						  hPrimaryAntiProtons->GetXaxis()->GetNbins(),
						  hPrimaryAntiProtons->GetXaxis()->GetXmin(),
						  hPrimaryAntiProtons->GetXaxis()->GetXmax());					      
  gYPrimaryAntiProtonsPercentage->Divide(gYPrimaryAntiProtons,
					 gYTotalAntiProtons,1.,1.0);
  SetError(gYPrimaryAntiProtonsPercentage,gYTotalAntiProtons);
  gYPrimaryAntiProtonsPercentage->Scale(100.);
  gYPrimaryAntiProtonsPercentage->SetMarkerStyle(kFullCircle);
  
  TH1D *gYSecondaryAntiProtonsPercentage = new TH1D("gYSecondaryAntiProtonsPercentage",
						    "",
						    hSecondaryAntiProtons->GetXaxis()->GetNbins(),
						    hSecondaryAntiProtons->GetXaxis()->GetXmin(),
						    hSecondaryAntiProtons->GetXaxis()->GetXmax());					      
  gYSecondaryAntiProtonsPercentage->Divide(gYSecondaryAntiProtons,
					   gYTotalAntiProtons,1.,1.0);
  SetError(gYSecondaryAntiProtonsPercentage,gYTotalAntiProtons);
  gYSecondaryAntiProtonsPercentage->Scale(100.);
  gYSecondaryAntiProtonsPercentage->SetMarkerStyle(kOpenCircle);
  
  
  TH2F *hEmptyY = new TH2F("hEmptyCompositionY","",
			   100,-1.2,1.2,100,-10.0,130); 
  hEmptyY->SetStats(kFALSE); 
  hEmptyY->GetYaxis()->SetTitle("Particle composition [%]");
  hEmptyY->GetYaxis()->SetTitleOffset(1.3);
  hEmptyY->GetXaxis()->SetTitle("y");

  TCanvas *c7 = new TCanvas("c7","(Anti)Proton contamination vs y",
			    150,150,700,400);
  c7->SetFillColor(10); c7->GetFrame()->SetFillColor(10); 
  c7->SetHighLightColor(10); c7->Divide(2,1);

  c7->cd(1)->SetBottomMargin(0.15); 
  c7->cd(1)->SetLeftMargin(0.15); 
  c7->cd(1)->SetGridx(); c7->cd(1)->SetGridy();
  hEmptyY->SetTitle("Protons");
  hEmptyY->DrawCopy();
  gYPrimaryProtonsPercentage->DrawCopy("ESAME");
  gYSecondaryProtonsPercentage->DrawCopy("ESAME");

  DrawMarker(0, 55, kFullCircle, 1.2, 1);
  t1->DrawLatex(0.1,53,"Primaries");
  DrawMarker(0, 45, kOpenCircle, 1.2, 1);
  t1->DrawLatex(0.1,43,"Secondaries");

  c7->cd(2)->SetBottomMargin(0.15); 
  c7->cd(2)->SetLeftMargin(0.15); 
  c7->cd(2)->SetGridx(); c7->cd(2)->SetGridy();
  hEmptyY->SetTitle("Antiprotons");
  hEmptyY->DrawCopy();
  gYPrimaryAntiProtonsPercentage->DrawCopy("ESAME");
  gYSecondaryAntiProtonsPercentage->DrawCopy("ESAME");

  DrawMarker(0, 55, kFullCircle, 1.2, 1);
  t1->DrawLatex(0.1,53,"Primaries");
  DrawMarker(0, 45, kOpenCircle, 1.2, 1);
  t1->DrawLatex(0.1,43,"Secondaries");

  c7->SaveAs("Contamination-Protons-Rapidity.gif");

  //pT dependence
  //Protons
  TH1D *gPtPrimaryProtons = (TH1D *)hPrimaryProtons->ProjectionY("gPtPrimaryProtons",0,hPrimaryProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtSecondaryProtons = (TH1D *)hSecondaryProtons->ProjectionY("gPtSecondaryProtons",0,hSecondaryProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtTotalProtons = (TH1D *)hPrimaryProtons->ProjectionY("gPtTotalProtons",0,hPrimaryProtons->GetYaxis()->GetNbins(),"e");
  gPtTotalProtons->Add(gPtSecondaryProtons);

  TH1D *gPtPrimaryProtonsPercentage = new TH1D("gPtPrimaryProtonsPercentage",
					       "",
					       hPrimaryProtons->GetYaxis()->GetNbins(),
					       hPrimaryProtons->GetYaxis()->GetXmin(),
					       hPrimaryProtons->GetYaxis()->GetXmax());					      
  gPtPrimaryProtonsPercentage->Divide(gPtPrimaryProtons,
				      gPtTotalProtons,1.,1.0);
  SetError(gPtPrimaryProtonsPercentage,gPtTotalProtons);
  gPtPrimaryProtonsPercentage->Scale(100.);
  gPtPrimaryProtonsPercentage->SetMarkerStyle(kFullCircle);
  
  TH1D *gPtSecondaryProtonsPercentage = new TH1D("gPtSecondaryProtonsPercentage",
						 "",
						 hSecondaryProtons->GetYaxis()->GetNbins(),
						 hSecondaryProtons->GetYaxis()->GetXmin(),
						 hSecondaryProtons->GetYaxis()->GetXmax());					      
  gPtSecondaryProtonsPercentage->Divide(gPtSecondaryProtons,
					gPtTotalProtons,1.,1.0);
  SetError(gPtSecondaryProtonsPercentage,gPtTotalProtons);
  gPtSecondaryProtonsPercentage->Scale(100.);
  gPtSecondaryProtonsPercentage->SetMarkerStyle(kOpenCircle);


  //Antiprotons
  TH1D *gPtPrimaryAntiProtons = (TH1D *)hPrimaryAntiProtons->ProjectionY("gPtPrimaryAntiProtons",0,hPrimaryAntiProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtSecondaryAntiProtons = (TH1D *)hSecondaryAntiProtons->ProjectionY("gPtSecondaryAntiProtons",0,hSecondaryAntiProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtTotalAntiProtons = (TH1D *)hPrimaryAntiProtons->ProjectionY("gPtTotalAntiProtons",0,hPrimaryAntiProtons->GetYaxis()->GetNbins(),"e");
  gPtTotalAntiProtons->Add(gPtSecondaryAntiProtons);
  
  TH1D *gPtPrimaryAntiProtonsPercentage = new TH1D("gPtPrimaryAntiProtonsPercentage",
						   "",
						   hPrimaryAntiProtons->GetYaxis()->GetNbins(),
						   hPrimaryAntiProtons->GetYaxis()->GetXmin(),
						   hPrimaryAntiProtons->GetYaxis()->GetXmax());					      
  gPtPrimaryAntiProtonsPercentage->Divide(gPtPrimaryAntiProtons,
					  gPtTotalAntiProtons,1.,1.0);
  SetError(gPtPrimaryAntiProtonsPercentage,gPtTotalAntiProtons);
  gPtPrimaryAntiProtonsPercentage->Scale(100.);
  gPtPrimaryAntiProtonsPercentage->SetMarkerStyle(kFullCircle);
  
  TH1D *gPtSecondaryAntiProtonsPercentage = new TH1D("gPtSecondaryAntiProtonsPercentage",
						     "",
						     hSecondaryAntiProtons->GetYaxis()->GetNbins(),
						     hSecondaryAntiProtons->GetYaxis()->GetXmin(),
						     hSecondaryAntiProtons->GetYaxis()->GetXmax());					      
  gPtSecondaryAntiProtonsPercentage->Divide(gPtSecondaryAntiProtons,
					    gPtTotalAntiProtons,1.,1.0);
  SetError(gPtSecondaryAntiProtonsPercentage,gPtTotalAntiProtons);
  gPtSecondaryAntiProtonsPercentage->Scale(100.);
  gPtSecondaryAntiProtonsPercentage->SetMarkerStyle(kOpenCircle);
  
  TH2F *hEmptyPt = new TH2F("hEmptyCompositionPt","",
			   100,0.0,4.0,100,-10.0,130); 
  hEmptyPt->SetStats(kFALSE); 
  hEmptyPt->GetYaxis()->SetTitle("Particle composition [%]");
  hEmptyPt->GetYaxis()->SetTitleOffset(1.3);
  hEmptyPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");

  TCanvas *c8 = new TCanvas("c8","(Anti)Proton comtamination vs pT",
			    200,200,700,400);
  c8->SetFillColor(10); c8->GetFrame()->SetFillColor(10); 
  c8->SetHighLightColor(10); c8->Divide(2,1);

  c8->cd(1)->SetBottomMargin(0.15); 
  c8->cd(1)->SetLeftMargin(0.15); 
  c8->cd(1)->SetGridx(); c8->cd(1)->SetGridy();
  hEmptyPt->SetTitle("Protons");
  hEmptyPt->DrawCopy();
  gPtPrimaryProtonsPercentage->DrawCopy("ESAME");
  gPtSecondaryProtonsPercentage->DrawCopy("ESAME");

  DrawMarker(2.0, 55, kFullCircle, 1.2, 1);
  t1->DrawLatex(2.1,53,"Primaries");
  DrawMarker(2.0, 45, kOpenCircle, 1.2, 1);
  t1->DrawLatex(2.1,43,"Secondaries");

  c8->cd(2)->SetBottomMargin(0.15); 
  c8->cd(2)->SetLeftMargin(0.15); 
  c8->cd(2)->SetGridx(); c8->cd(2)->SetGridy();
  hEmptyPt->SetTitle("Antirotons");
  hEmptyPt->DrawCopy();
  gPtPrimaryAntiProtonsPercentage->DrawCopy("ESAME");
  gPtSecondaryAntiProtonsPercentage->DrawCopy("ESAME");

  DrawMarker(2.0, 55, kFullCircle, 1.2, 1);
  t1->DrawLatex(2.1,53,"Primaries");
  DrawMarker(2.0, 45, kOpenCircle, 1.2, 1);
  t1->DrawLatex(2.1,43,"Secondaries");

  c8->SaveAs("Contamination-Protons-Pt.gif");

  TFile *fout = TFile::Open("Contamination.root","recreate");
  gYPrimaryProtonsPercentage->Write();
  gYSecondaryProtonsPercentage->Write();
  gPtPrimaryProtonsPercentage->Write();
  gPtSecondaryProtonsPercentage->Write();
  gYPrimaryAntiProtonsPercentage->Write();
  gYSecondaryAntiProtonsPercentage->Write();
  gPtPrimaryAntiProtonsPercentage->Write();
  gPtSecondaryAntiProtonsPercentage->Write();
  fout->Close();
}

//________________________________________________//
void DrawCutEfficiency(TList *inputList) {
  //loops over the list entries and
  //draws the rapidity and pT dependence
  //of the percentage of primary and secondary
  //protons and antiprotons after the track cuts
  cout<<"Extracting the entries for the histograms in the list: "<<
    inputList->GetName()<<"..."<<endl;

  TLatex *t1 = new TLatex();
  t1->SetTextSize(0.04);

  TH2F *hPrimaryESDProtons = (TH2F *)inputList->At(0);
  TH2F *hPrimaryESDAntiProtons = (TH2F *)inputList->At(4);
  TH2F *hPrimaryMCProtons = (TH2F *)inputList->At(8);
  TH2F *hPrimaryMCAntiProtons = (TH2F *)inputList->At(9);

  //rapidity dependence
  //Protons
  TH1D *gYPrimaryESDProtons = (TH1D *)hPrimaryESDProtons->ProjectionX("gYPrimaryESDProtons",0,hPrimaryESDProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYPrimaryMCProtons = (TH1D *)hPrimaryMCProtons->ProjectionX("gYPrimaryMCProtons",0,hPrimaryMCProtons->GetXaxis()->GetNbins(),"e");
  gYPrimaryESDProtons->Divide(gYPrimaryMCProtons);
  SetError(gYPrimaryESDProtons,gYPrimaryMCProtons);
  gYPrimaryESDProtons->Scale(100.);
  gYPrimaryESDProtons->SetMarkerStyle(kFullCircle);

  //Antiprotons
  TH1D *gYPrimaryESDAntiProtons = (TH1D *)hPrimaryESDAntiProtons->ProjectionX("gYPrimaryESDAntiProtons",0,hPrimaryESDAntiProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYPrimaryMCAntiProtons = (TH1D *)hPrimaryMCAntiProtons->ProjectionX("gYPrimaryMCAntiProtons",0,hPrimaryMCAntiProtons->GetXaxis()->GetNbins(),"e");
  gYPrimaryESDAntiProtons->Divide(gYPrimaryMCAntiProtons);
  SetError(gYPrimaryESDAntiProtons,gYPrimaryMCAntiProtons);
  gYPrimaryESDAntiProtons->Scale(100.);
  gYPrimaryESDAntiProtons->SetMarkerStyle(kFullCircle);
  
  TH2F *hEmptyY = new TH2F("hEmptyEfficiencyY","",
			   100,-1.2,1.2,100,-10.0,130); 
  hEmptyY->SetStats(kFALSE); 
  hEmptyY->GetYaxis()->SetTitle("#epsilon [%]");
  hEmptyY->GetYaxis()->SetTitleOffset(1.3);
  hEmptyY->GetXaxis()->SetTitle("y");

  TCanvas *c10 = new TCanvas("c10","(Anti)Proton cut efficiency vs y",
			    250,250,700,400);
  c10->SetFillColor(10); c10->GetFrame()->SetFillColor(10); 
  c10->SetHighLightColor(10); c10->Divide(2,1);

  c10->cd(1)->SetBottomMargin(0.15); 
  c10->cd(1)->SetLeftMargin(0.15); 
  c10->cd(1)->SetGridx(); c10->cd(1)->SetGridy();
  hEmptyY->SetTitle("Protons");
  hEmptyY->DrawCopy();
  gYPrimaryESDProtons->DrawCopy("ESAME");

  c10->cd(2)->SetBottomMargin(0.15); 
  c10->cd(2)->SetLeftMargin(0.15); 
  c10->cd(2)->SetGridx(); c10->cd(2)->SetGridy();
  hEmptyY->SetTitle("Antiprotons");
  hEmptyY->DrawCopy();
  gYPrimaryESDAntiProtons->DrawCopy("ESAME");

  c10->SaveAs("CutEfficiency-Protons-Rapidity.gif");

  //pT dependence
  //Protons
  TH1D *gPtPrimaryESDProtons = (TH1D *)hPrimaryESDProtons->ProjectionY("gPtPrimaryESDProtons",0,hPrimaryESDProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtPrimaryMCProtons = (TH1D *)hPrimaryMCProtons->ProjectionY("gPtPrimaryMCProtons",0,hPrimaryMCProtons->GetYaxis()->GetNbins(),"e");
  gPtPrimaryESDProtons->Divide(gPtPrimaryMCProtons);
  SetError(gPtPrimaryESDProtons,gPtPrimaryMCProtons);
  gPtPrimaryESDProtons->Scale(100.);
  gPtPrimaryESDProtons->SetMarkerStyle(kFullCircle);

  //Antiprotons
  TH1D *gPtPrimaryESDAntiProtons = (TH1D *)hPrimaryESDAntiProtons->ProjectionY("gPtPrimaryESDAntiProtons",0,hPrimaryESDAntiProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtPrimaryMCAntiProtons = (TH1D *)hPrimaryMCAntiProtons->ProjectionY("gPtPrimaryMCAntiProtons",0,hPrimaryMCAntiProtons->GetYaxis()->GetNbins(),"e");
  gPtPrimaryESDAntiProtons->Divide(gPtPrimaryMCAntiProtons);
  SetError(gPtPrimaryESDAntiProtons,gPtPrimaryMCAntiProtons);
  gPtPrimaryESDAntiProtons->Scale(100.);
  gPtPrimaryESDAntiProtons->SetMarkerStyle(kFullCircle);

  TH2F *hEmptyPt = new TH2F("hEmptyEfficiencyPt","",
			   100,0.0,4.0,100,-10.0,130); 
  hEmptyPt->SetStats(kFALSE); 
  hEmptyPt->GetYaxis()->SetTitle("#epsilon [%]");
  hEmptyPt->GetYaxis()->SetTitleOffset(1.3);
  hEmptyPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");

  TCanvas *c11 = new TCanvas("c11","(Anti)Proton cut efficiency vs pT",
			    300,300,700,400);
  c11->SetFillColor(10); c11->GetFrame()->SetFillColor(10); 
  c11->SetHighLightColor(10); c11->Divide(2,1);

  c11->cd(1)->SetBottomMargin(0.15); 
  c11->cd(1)->SetLeftMargin(0.15); 
  c11->cd(1)->SetGridx(); c11->cd(1)->SetGridy();
  hEmptyPt->SetTitle("Protons");
  hEmptyPt->DrawCopy();
  gPtPrimaryESDProtons->DrawCopy("ESAME");

  c11->cd(2)->SetBottomMargin(0.15); 
  c11->cd(2)->SetLeftMargin(0.15); 
  c11->cd(2)->SetGridx(); c11->cd(2)->SetGridy();
  hEmptyPt->SetTitle("Antirotons");
  hEmptyPt->DrawCopy();
  gPtPrimaryESDAntiProtons->DrawCopy("ESAME");

  c11->SaveAs("CutEfficiency-Protons-Pt.gif");

  TFile *fout = TFile::Open("Efficiency.root","recreate");
  gYPrimaryESDProtons->Write();
  gYPrimaryESDAntiProtons->Write();
  gPtPrimaryESDProtons->Write();
  gPtPrimaryESDAntiProtons->Write();
  fout->Close();
}

//________________________________________________//
void GetQAEntries(TList *inputList, Double_t *entries) {
  //loops over the list entries
  //extracts the entries for each histogram
  //cout<<"Extracting the entries for the histograms in the list: "<<
  //inputList->GetName()<<"..."<<endl;

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

//________________________________________//
void drawMCQA(TList *listPDG, TList *listMCProcesses) {
  //Function to display the composition of secondary (anti)protons
  //The histogram shows the percentage of secondary (anti)protons
  //originating from each particle species.
  //The box summarizes the MC process that gave these secondary (anti)protons
  TDatabasePDG *db = TDatabasePDG::Instance();
  TParticlePDG *p = 0x0;
  
  TH3F *gHistYPtPDGProtons = (TH3F *)listPDG->At(0);
  TH3F *gHistYPtPDGAntiProtons = (TH3F *)listPDG->At(1);
  readProcesses(listMCProcesses);
  Double_t nParticleCompositionProtonY[100], nParticleCompositionProtonPt[100];
  Double_t nParticleCompositionAntiProtonY[100], nParticleCompositionAntiProtonPt[100];
  Double_t gY[100], gPt[100];
  for(Int_t iBins = 0; iBins < 100; iBins++) {
    nParticleCompositionProtonY[iBins] = 0;
    nParticleCompositionProtonPt[iBins] = 0;
    nParticleCompositionAntiProtonY[iBins] = 0;
    nParticleCompositionAntiProtonPt[iBins] = 0;
    gY[iBins] = 0;
    gPt[iBins] = 0;
  }
  
  TGraph *gParticleProtonY[14];
  TGraph *gParticleProtonPt[14];
  TGraph *gParticleAntiProtonY[14];
  TGraph *gParticleAntiProtonPt[14];
  for(Int_t iParticle = 0; iParticle < 14; iParticle++) {
    GetComposition(iParticle,
		   gHistYPtPDGProtons,
		   nParticleCompositionProtonY,
		   gY, nParticleCompositionProtonPt, gPt);
    gParticleProtonY[iParticle] = new TGraph(gHistYPtPDGProtons->GetNbinsX(),
				       gY,nParticleCompositionProtonY);
    gParticleProtonY[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleProtonY[iParticle]->SetMarkerSize(1.2);

    gParticleProtonPt[iParticle] = new TGraph(gHistYPtPDGProtons->GetNbinsY(),
					gPt,nParticleCompositionProtonPt);
    gParticleProtonPt[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleProtonPt[iParticle]->SetMarkerSize(1.2);

    GetComposition(iParticle,
		   gHistYPtPDGAntiProtons,
		   nParticleCompositionAntiProtonY,
		   gY, nParticleCompositionAntiProtonPt, gPt);
    gParticleAntiProtonY[iParticle] = new TGraph(gHistYPtPDGAntiProtons->GetNbinsX(),
				       gY,nParticleCompositionAntiProtonY);
    gParticleAntiProtonY[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleAntiProtonY[iParticle]->SetMarkerSize(1.2);

    gParticleAntiProtonPt[iParticle] = new TGraph(gHistYPtPDGAntiProtons->GetNbinsY(),
					gPt,nParticleCompositionAntiProtonPt);
    gParticleAntiProtonPt[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleAntiProtonPt[iParticle]->SetMarkerSize(1.2);
  }

  //_________________________________________________________//
  char *fParticleName[14] = {"Primary","K_{L}","#pi","K_{S}","K",
			     "n","p","#Sigma^{-}","#Lambda","#Sigma^{+}",
			     "#Xi^{-}","#Xi^{0}","#Omega^{-}"};
  TLatex *t1 = new TLatex();
  t1->SetTextSize(0.04);

  TH2F *hEmptyY = new TH2F("hEmptyY","",100,-1.2,1.2,100,0,120); 
  hEmptyY->SetStats(kFALSE); 
  hEmptyY->GetYaxis()->SetTitle("Particle composition [%]");
  hEmptyY->GetXaxis()->SetTitle("y");

  TCanvas *c3 = new TCanvas("c3","MC secondary composition vs y - Protons",
			    350,350,700,400);
  c3->SetFillColor(10); c3->GetFrame()->SetFillColor(10);
  c3->SetHighLightColor(10); c3->SetBottomMargin(0.15);
  c3->SetGridx(); c3->SetGridy();
  hEmptyY->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    //if((iParticle == 0)||(iParticle == 2)||(iParticle == 5)||(iParticle == 6)||(iParticle == 8))
    if((iParticle == 0)||(iParticle == 2)||(iParticle == 6)||(iParticle == 8))
      gParticleProtonY[iParticle]->Draw("P");
    if(iParticle < 5) {
      DrawMarker(-1.1, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(-1.0,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(0.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*(iParticle-5),fParticleName[iParticle]);
    }
  }

  TCanvas *c5 = new TCanvas("c5","MC secondary composition vs y - antiProtons",
			    400,400,700,400);
  c5->SetFillColor(10); c5->GetFrame()->SetFillColor(10);
  c5->SetHighLightColor(10); c5->SetBottomMargin(0.15);
  c5->SetGridx(); c5->SetGridy();
  hEmptyY->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    if((iParticle == 0)||(iParticle == 6)||(iParticle == 8))
      gParticleAntiProtonY[iParticle]->Draw("P");
    if(iParticle < 5) {
      DrawMarker(-1.1, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(-1.0,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(0.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*(iParticle-5),fParticleName[iParticle]);
    }
  }

  TH2F *hEmptyPt = new TH2F("hEmptyPt","",100,0.0,4.0,100,0,120); 
  hEmptyPt->SetStats(kFALSE); 
  hEmptyPt->GetYaxis()->SetTitle("Particle composition [%]");
  hEmptyPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");

  TCanvas *c4 = new TCanvas("c4","MC secondary composition vs pT - Protons",
			    450,450,700,400);
  c4->SetFillColor(10); c4->GetFrame()->SetFillColor(10);
  c4->SetHighLightColor(10); c4->SetBottomMargin(0.15);
  c4->SetGridx(); c4->SetGridy();
  hEmptyPt->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    if(iParticle < 5) {
      DrawMarker(0.2, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(2.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(2.3,113-5*(iParticle-5),fParticleName[iParticle]);
    }
    if((iParticle == 0)||(iParticle == 2)||(iParticle == 6)||(iParticle == 8))
      gParticleProtonPt[iParticle]->Draw("P");
  }

  TCanvas *c6 = new TCanvas("c6",
			    "MC secondary composition vs pT - AntiProtons",
			    500,500,700,400);
  c6->SetFillColor(10); c6->GetFrame()->SetFillColor(10);
  c6->SetHighLightColor(10); c6->SetBottomMargin(0.15);
  c6->SetGridx(); c6->SetGridy();
  hEmptyPt->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    if(iParticle < 5) {
      DrawMarker(0.2, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(2.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(2.3,113-5*(iParticle-5),fParticleName[iParticle]);
    }
    if((iParticle == 0)||(iParticle == 6)||(iParticle == 8))
      gParticleAntiProtonPt[iParticle]->Draw("P");
  }
}

//________________________________________//
void GetComposition(Int_t iSpecies,
		    TH3F *gHist, 
		    Double_t *nParticleCompositionY,
		    Double_t *gY,
		    Double_t *nParticleCompositionPt,
		    Double_t *gPt) {
  //Returns the pT and y dependence of the MC composition
  Double_t ymin = gHist->GetXaxis()->GetXmin();
  Double_t ymax = gHist->GetXaxis()->GetXmax();
  Double_t nybins = gHist->GetNbinsX();
  Double_t ptmin = gHist->GetYaxis()->GetXmin();
  Double_t ptmax = gHist->GetYaxis()->GetXmax();
  Double_t nptbins = gHist->GetNbinsY();
  Double_t nTotalY[100], nTotalPt[100];
  for(Int_t iBins = 0; iBins < 100; iBins++) {
    nParticleCompositionY[iBins] = 0;
    nParticleCompositionPt[iBins] = 0;
    nTotalY[iBins] = 0.0;
    nTotalPt[iBins] = 0.0;
  }

  //rapidity dependence
  //cout<<"Ymin: "<<ymin<<" - Ymax: "<<ymax<<" - Ybins: "<<nybins<<endl;
  for(Int_t iXbins = 1; iXbins <= gHist->GetNbinsX(); iXbins++) {
    for(Int_t iZbins = 1; iZbins <= gHist->GetNbinsZ(); iZbins++) {
      for(Int_t iYbins = 1; iYbins <= gHist->GetNbinsY(); iYbins++) {
	nTotalY[iXbins-1] += gHist->GetBinContent(iXbins,iYbins,iZbins);
      }
    }
  }

  for(Int_t iXbins = 1; iXbins <= gHist->GetNbinsX(); iXbins++) {
    for(Int_t iYbins = 1; iYbins <= gHist->GetNbinsY(); iYbins++) {
      if(nTotalY[iXbins-1] > 0)
	nParticleCompositionY[iXbins-1] += 100.*gHist->GetBinContent(iXbins,iYbins,iSpecies+1)/nTotalY[iXbins-1];
      //if(nParticleCompositionY[iXbins-1] == 0) 
      //nParticleCompositionY[iXbins-1] = -10.0;
    }//pt loop
    gY[iXbins-1] = ymin + (iXbins-1)*(ymax - ymin)/nybins + 0.5*(ymax - ymin)/nybins;
    //cout<<"y: "<<gY[iXbins-1]<<
    //" - test: "<<ymin + (iXbins-1)*(ymax - ymin)/nybins + 0.5*(ymax - ymin)/nybins<<
    //" - Number of protons: "<<nY[iXbins-1]<<
    //" - Total: "<<nTotalY[iXbins-1]<<
    //" - Percentage: "<<nParticleCompositionY[iXbins-1]<<endl;
  }//y loop

  //pt dependence
  //cout<<"Ptmin: "<<ptmin<<" - Ptmax: "<<ptmax<<" - Ptbins: "<<nptbins<<endl;
  for(Int_t iYbins = 1; iYbins <= gHist->GetNbinsY(); iYbins++) {
    for(Int_t iZbins = 1; iZbins <= gHist->GetNbinsZ(); iZbins++) {
      for(Int_t iXbins = 1; iXbins <= gHist->GetNbinsX(); iXbins++) {
	nTotalPt[iYbins-1] += gHist->GetBinContent(iXbins,iYbins,iZbins);
      }
    }
  }

  for(Int_t iYbins = 1; iYbins <= gHist->GetNbinsY(); iYbins++) {
    for(Int_t iXbins = 1; iXbins <= gHist->GetNbinsX(); iXbins++) {
      if(nTotalPt[iYbins-1] > 0)
	nParticleCompositionPt[iYbins-1] += 100.*gHist->GetBinContent(iXbins,iYbins,iSpecies+1)/nTotalPt[iYbins-1];
      //if(nParticleCompositionPt[iYbins-1] == 0) 
      //nParticleCompositionPt[iYbins-1] = -10.0;
    }//pt loop
    gPt[iYbins-1] = ptmin + (iYbins-1)*(ptmax - ptmin)/nptbins + 0.5*(ptmax - ptmin)/nptbins;
    //cout<<"Pt: "<<gPt[iYbins-1]<<
    //" - test: "<<ptmin + (iYbins-1)*(ptmax - ptmin)/nptbins + 0.5*(ptmax - ptmin)/nptbins<<
    //" - Number of protons: "<<nY[iXbins-1]<<
    //" - Total: "<<nTotalPt[iYbins-1]<<
    //" - Percentage: "<<nParticleCompositionPt[iYbins-1]<<endl;
  }//pt loop
}

//________________________________________//
void readProcesses(TList *list) {
  char *fParticleProtonName[12] = {"K_{L}","#pi","K_{S}","K",
				   "n","p","#Sigma^{-}","#Lambda","#Sigma^{+}",
				   "#Xi^{-}","#Xi^{0}","#Omega^{-}"};
  char *fParticleAntiProtonName[8] = {"K_{L}","#pi","K_{S}","K",
				      "n","p","#Lambda","#Sigma^{+}"};
  Int_t iProtonCounter = 0, iAntiProtonCounter = 0;
  TH1F *gMCProcesses;
  for(Int_t iEntry = 0; iEntry < list->GetEntries(); iEntry++) {
    gMCProcesses = (TH1F *)list->At(iEntry);
    TString histName = gMCProcesses->GetName();
    if(histName.Contains("gHistProtons")) {
      cout<<"Protons coming from "<<fParticleProtonName[iProtonCounter]<<endl;

      iProtonCounter += 1;
    }
    if(histName.Contains("gHistAntiProtons")) {
      cout<<"Antiprotons coming from "<<fParticleAntiProtonName[iAntiProtonCounter]<<endl;

      iAntiProtonCounter += 1;
    }
    for(Int_t iBin = 1; iBin < gMCProcesses->GetNbinsX(); iBin++) {
      Double_t binContent = gMCProcesses->GetBinContent(iBin);
      if(binContent > 0) {
	Int_t processId = gMCProcesses->GetBinCenter(iBin);
	cout<<"\t Process ID: "<<processId<<" - "<<
	  gMCProcessName[processId]<<" - "<<
	  100.*binContent/gMCProcesses->GetEntries()<<"%"<<endl;
      }
    }
  }
}

//________________________________________________//
void SetError(TH1D *hEff, TH1D *hGen) {
  for(Int_t iBin = 1; iBin <= hEff->GetNbinsX(); iBin++) {
    Double_t error = 0.0;
    if(hEff->GetBinContent(iBin) <= 1  .) 
      error = TMath::Sqrt(hEff->GetBinContent(iBin)*(1  . - hEff->GetBinContent(iBin))/hGen->GetBinContent(iBin));
    hEff->SetBinError(iBin,error);
  }
}

//________________________________________//
void drawEfficiency(TList *list) {
  //Function to display the reconstruction and PID efficiencies
  //for protons and antiprotons vs y and pT

  TH2F *hEmpty = new TH2F("hEmptyReconstructionEfficiency","",
			   100,-1.2,3.5,100,-10.0,130); 
  hEmpty->SetStats(kFALSE); 
  hEmpty->GetYaxis()->SetTitle("#epsilon [%]");
  hEmpty->GetYaxis()->SetTitleOffset(1.3);

  //Reconstruction efficiency
  TH2D *gHistMCYPtProtons = (TH2D *)list->At(0);
  TH2D *gHistMCYPtAntiProtons = (TH2D *)list->At(1);
  TH2D *gHistESDYPtProtons = (TH2D *)list->At(2);
  TH2D *gHistESDYPtAntiProtons = (TH2D *)list->At(3);

  //rapidity dependence
  TCanvas *c14 = new TCanvas("c14",
			     "(Anti)Proton reconstruction efficiency vs y",
			     600,600,700,400);
  c14->SetFillColor(10); c14->GetFrame()->SetFillColor(10); 
  c14->SetHighLightColor(10); c14->Divide(2,1);

  //Protons
  TH1D *gYESDProtons = (TH1D *)gHistESDYPtProtons->ProjectionX("gYESDProtons",0,gHistESDYPtProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYMCProtons = (TH1D *)gHistMCYPtProtons->ProjectionX("gYMCProtons",0,gHistMCYPtProtons->GetXaxis()->GetNbins(),"e");
  gYESDProtons->Divide(gYMCProtons);
  SetError(gYESDProtons,gYMCProtons);
  gYESDProtons->Scale(100.);
  gYESDProtons->SetMarkerStyle(kFullCircle);

  //AntiProtons
  TH1D *gYESDAntiProtons = (TH1D *)gHistESDYPtAntiProtons->ProjectionX("gYESDAntiProtons",0,gHistESDYPtAntiProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYMCAntiProtons = (TH1D *)gHistMCYPtAntiProtons->ProjectionX("gYMCAntiProtons",0,gHistMCYPtProtons->GetXaxis()->GetNbins(),"e");
  gYESDAntiProtons->Divide(gYMCAntiProtons);
  SetError(gYESDAntiProtons,gYMCAntiProtons);
  gYESDAntiProtons->Scale(100.);
  gYESDAntiProtons->SetMarkerStyle(kFullCircle);

  c14->cd(1)->SetBottomMargin(0.15); 
  c14->cd(1)->SetLeftMargin(0.15); 
  c14->cd(1)->SetGridx(); c14->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(-1.0,1.0);
  hEmpty->GetXaxis()->SetTitle("y");
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gYESDProtons->DrawCopy("ESAME");

  c14->cd(2)->SetBottomMargin(0.15); 
  c14->cd(2)->SetLeftMargin(0.15); 
  c14->cd(2)->SetGridx(); c14->cd(2)->SetGridy();
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gYESDAntiProtons->DrawCopy("ESAME");
  c14->SaveAs("ReconstructionEfficiency-Protons-Rapidity.gif");

  //pT dependence
  TCanvas *c15 = new TCanvas("c15",
			     "(Anti)Proton reconstruction efficiency vs pT",
			     650,650,700,400);
  c15->SetFillColor(10); c15->GetFrame()->SetFillColor(10); 
  c15->SetHighLightColor(10); c15->Divide(2,1);

  //Protons
  TH1D *gPtESDProtons = (TH1D *)gHistESDYPtProtons->ProjectionY("gPtESDProtons",0,gHistESDYPtProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtMCProtons = (TH1D *)gHistMCYPtProtons->ProjectionY("gPtMCProtons",0,gHistMCYPtProtons->GetYaxis()->GetNbins(),"e");
  gPtESDProtons->Divide(gPtMCProtons);
  SetError(gPtESDProtons,gPtMCProtons);
  gPtESDProtons->Scale(100.);
  gPtESDProtons->SetMarkerStyle(kFullCircle);

  //AntiProtons
  TH1D *gPtESDAntiProtons = (TH1D *)gHistESDYPtAntiProtons->ProjectionY("gPtESDAntiProtons",0,gHistESDYPtAntiProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtMCAntiProtons = (TH1D *)gHistMCYPtAntiProtons->ProjectionY("gPtMCAntiProtons",0,gHistMCYPtProtons->GetYaxis()->GetNbins(),"e");
  gPtESDAntiProtons->Divide(gPtMCAntiProtons);
  SetError(gPtESDAntiProtons,gPtMCAntiProtons);
  gPtESDAntiProtons->Scale(100.);
  gPtESDAntiProtons->SetMarkerStyle(kFullCircle);


  c15->cd(1)->SetBottomMargin(0.15); 
  c15->cd(1)->SetLeftMargin(0.15); 
  c15->cd(1)->SetGridx(); c15->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(0.,1.2);
  hEmpty->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gPtESDProtons->DrawCopy("ESAME");

  c15->cd(2)->SetBottomMargin(0.15); 
  c15->cd(2)->SetLeftMargin(0.15); 
  c15->cd(2)->SetGridx(); c15->cd(2)->SetGridy();
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gPtESDAntiProtons->DrawCopy("ESAME");
  c15->SaveAs("ReconstructionEfficiency-Protons-Pt.gif");

  //PID efficiency
  TH2D *gHistESDInitYPtProtons = (TH2D *)list->At(4);
  TH2D *gHistESDIdYPtProtons = (TH2D *)list->At(5);
  TH2D *gHistESDRecIdYPtProtons = (TH2D *)list->At(6);
  TH2D *gHistESDContamYPtProtons = (TH2D *)list->At(7);

  TCanvas *c16 = new TCanvas("c16",
			     "(Anti)Proton PID efficiency vs y and pT",
			     700,700,700,400);
  c16->SetFillColor(10); c16->GetFrame()->SetFillColor(10); 
  c16->SetHighLightColor(10); c16->Divide(2,1);

  //rapidity dependence
  //protons pid efficiency
  TH1D *gYESDIdProtons = (TH1D *)gHistESDIdYPtProtons->ProjectionX("gYESDIdProtons",0,gHistESDIdYPtProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYESDInitProtons = (TH1D *)gHistESDInitYPtProtons->ProjectionX("gYESDInitProtons",0,gHistESDInitYPtProtons->GetXaxis()->GetNbins(),"e");
  gYESDIdProtons->Divide(gYESDInitProtons);
  SetError(gYESDIdProtons,gYESDInitProtons);
  gYESDIdProtons->Scale(100.);
  gYESDIdProtons->SetMarkerStyle(kFullCircle);

  //protons pid contamination
  TH1D *gYESDContamProtons = (TH1D *)gHistESDContamYPtProtons->ProjectionX("gYESDContamProtons",0,gHistESDContamYPtProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYESDRecIdProtons = (TH1D *)gHistESDRecIdYPtProtons->ProjectionX("gYESDRecIdProtons",0,gHistESDRecIdYPtProtons->GetXaxis()->GetNbins(),"e");
  gYESDContamProtons->Divide(gYESDRecIdProtons);
  SetError(gYESDContamProtons,gYESDRecIdProtons);
  gYESDContamProtons->Scale(100.);
  gYESDContamProtons->SetMarkerStyle(kOpenCircle);

  c16->cd(1)->SetBottomMargin(0.15); 
  c16->cd(1)->SetLeftMargin(0.15); 
  c16->cd(1)->SetGridx(); c16->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(-1.0.,1.0);
  hEmpty->GetXaxis()->SetTitle("y");
  hEmpty->DrawCopy();
  gYESDIdProtons->DrawCopy("ESAME");
  gYESDContamProtons->DrawCopy("ESAME");

  //pT dependence
  //protons pid efficiency
  TH1D *gPtESDIdProtons = (TH1D *)gHistESDIdYPtProtons->ProjectionY("gPtESDIdProtons",0,gHistESDIdYPtProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtESDInitProtons = (TH1D *)gHistESDInitYPtProtons->ProjectionY("gPtESDInitProtons",0,gHistESDInitYPtProtons->GetYaxis()->GetNbins(),"e");
  gPtESDIdProtons->Divide(gPtESDInitProtons);
  SetError(gPtESDIdProtons,gPtESDInitProtons);
  gPtESDIdProtons->Scale(100.);
  gPtESDIdProtons->SetMarkerStyle(kFullCircle);

  //protons pid contamination
  TH1D *gPtESDContamProtons = (TH1D *)gHistESDContamYPtProtons->ProjectionY("gPtESDContamProtons",0,gHistESDContamYPtProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtESDRecIdProtons = (TH1D *)gHistESDRecIdYPtProtons->ProjectionY("gPtESDRecIdProtons",0,gHistESDRecIdYPtProtons->GetYaxis()->GetNbins(),"e");
  gPtESDContamProtons->Divide(gPtESDRecIdProtons);
  SetError(gPtESDContamProtons,gPtESDRecIdProtons);
  gPtESDContamProtons->Scale(100.);
  gPtESDContamProtons->SetMarkerStyle(kOpenCircle);

  c16->cd(2)->SetBottomMargin(0.15); 
  c16->cd(2)->SetLeftMargin(0.15); 
  c16->cd(2)->SetGridx(); c16->cd(2)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(0.0,1.2);
  hEmpty->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  hEmpty->DrawCopy();
  gPtESDIdProtons->DrawCopy("ESAME");
  gPtESDContamProtons->DrawCopy("ESAME");

  c16->SaveAs("PIDEfficiency-Protons.gif");
}

//________________________________________________//
void DrawMarker(Double_t x, Double_t y, Int_t style, 
		Double_t size, Int_t color) {
  TMarker *m = new TMarker(x,y,style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(color);
  m->Draw();
}

//________________________________________________//
const char * const gMCProcessName[45] = {
  "Primary particle emission",
  "Multiple scattering",
  "Energy loss",
  "Bending in magnetic field",
  "Decay",
  "Lepton pair production",
  "Compton scattering",
  "Photoelectric effect",
  "Bremstrahlung",
  "Delta ray",
  "Positron annihilation",
  "Positron annihilation at rest",
  "Positron annihilation in flight",
  "Hadronic interaction",
  "Nuclear evaporation",
  "Nuclear fission",
  "Nuclear absorbtion",
  "Antiproton annihilation",
  "Antineutron annihilation",
  "Neutron capture",
  "Hadronic elastic",
  "Hadronic incoherent elastic",
  "Hadronic coherent elastic",
  "Hadronic inelastic",
  "Photon inelastic",
  "Muon nuclear interaction",
  "Electron nuclear interaction",
  "Positron nuclear interaction",
  "Time of flight limit",
  "Nuclear photofission",
  "Rayleigh effect",
  "No active process",
  "Energy threshold",
  "Light absorption",
  "Light detection",
  "Light scattering",
  "Maximum allowed step",
  "Cerenkov production",
  "Cerenkov feed back photon",
  "Cerenkov photon reflection",
  "Cerenkov photon refraction",
  "Synchrotron radiation",
  "Scintillation",
  "Transportation",
  "Unknown process"
};
