void drawProtonQAResults(const char *analysisType = "TPC") {
  //Macro to visualize the results of the proton QA task
  //TCanvas objects: 20
  gStyle->SetPalette(1,0);
  gStyle->SetCanvasColor(41);
  gStyle->SetFrameFillColor(10);
  
  PrintHelpMenu();
  
  TString filename1 = "Protons.QA.";
  filename1 += analysisType; filename1 += ".root";
  TString filename2 = "Protons.MC.QA.";
  filename2 += analysisType; filename2 += ".root";
  TString filename3 = "Protons.Efficiency.";
  filename3 += analysisType; filename3 += ".root";
  
  TFile *fQA = TFile::Open(filename1.Data());
  TList *listGlobalQA = (TList *)fQA->Get("globalQAList");
  drawCutStatistics(listGlobalQA,analysisType);
  fQA->Close();
 
  TFile *fMC = TFile::Open(filename2.Data());
  TList *listPDG = (TList *)fMC->Get("pdgCodeList");  
  TList *listMCProcesses = (TList *)fMC->Get("mcProcessList");  
  drawMCQA(listPDG,listMCProcesses);
  fMC->Close();
  
  TFile *fEfficiency = TFile::Open(filename3.Data());
  TList *listEfficiency = (TList *)fEfficiency->Get("efficiencyList");
  drawEfficiency(listEfficiency,analysisType);
  fEfficiency->Close();

  TFile *fVertex = TFile::Open("Vertex.QA.root");
  TList *listVertex = (TList *)fVertex->Get("vertexList");
  drawVertexQA(listVertex);
  fVertex->Close();
}

//________________________________________//
void PrintHelpMenu() {
  //Function that prints the different possibilities the user has 
  //with this macro (what to display, what is the input etc)
  Printf("==============================HELP MENU===========================");
  Printf("Function drawCutStatistics: Takes as an argument the list obtained from the Protons.QA.root file and the corresponding analysis type. Displays the following plots:");
  Printf("\t 1. Influence of the track cuts to the primary and secondary protons. For each bin we show the percentage of particles rejected by the corresponding cut. The first column shows the percentage of tracks excluded from the analysis.");
  Printf("\t 2. Influence of the track cuts to the primary and secondary antiprotons. For each bin we show the percentage of particles rejected by the corresponding cut. The first column shows the percentage of tracks excluded from the analysis.");
  Printf("\t 3. The ITS cluster map for primary and secondary (anti)protons.");
  Printf("\t 4. The (anti)proton purity and contamination as a function of eta or y.");
  Printf("\t 5. The (anti)proton purity and contamination as a function of Pt.");
  Printf("\t 6. The (anti)proton cut efficiency as a function of eta or y.");
  Printf("\t 7. The (anti)proton cut efficiency as a function of Pt.");
  Printf("\t 8. The composition of accepted secondary (anti)protons as a function of eta or y.");
  Printf("\t 9. The composition of accepted secondary (anti)protons as a function of Pt.");
  Printf("\t ***This function calls the DrawContamination function that creates an output root file with the name Contamination.root holding all the purity and contamination plots (eta or y and Pt dependence) and the DrawCutEfficiency function that creates an output root file with the name CutEfficiency.root with the y or eta and Pt dependence of the cut efficiency for protons and antiprotons.***");

  Printf("\nFunction drawMCQA: Takes as an argument the two lists obtained from the Protons.MC.QA.root file holding pure MC information. Displays the following plots:");
  Printf("\t 1. MC composition of secondary protons as a function of either eta or y.");
  Printf("\t 2. MC composition of secondary antiprotons as a function of either eta or y.");
  Printf("\t 3. MC composition of secondary protons as a function of Pt.");
  Printf("\t 4. MC composition of secondary antiprotons as a function of Pt");

  Printf("\nFunction drawEfficiency: Takes as an argument the list obtained from the Protons.Efficiency.root file along with the analysis type. Displays the following plots:");
  Printf("\t 1. Reconstruction efficiencies for primary and secondary (anti)protons as a function of either eta or y. The secondaries are categorized by the MC process as originating from a weak decay or from a hadronic interaction.");
  Printf("\t 2. Reconstruction efficiencies for primary and secondary (anti)protons as a function of Pt. The secondaries are categorized by the MC process as originating from a weak decay or from a hadronic interaction.");
  Printf("\t 3. Particle identifiction efficiencies for (anti)protons as a function of either eta or y and of Pt.");
  Printf("\t ***This function produces an output root file with the name Reconstruction-PID-Efficiency.root with all the efficiency and contamination plots.***");

  Printf("\nFunction drawVertexQA: Takes as an argument the list obtained from the Vertex.QA.root file holding pure vertex information. Displays the following plots:");
  Printf("\t 1. The vertex efficiencies for the TPC, SPD and track ones.");
  Printf("\t 2. The QA plots for the TPC vertex.");
  Printf("\t 3. The QA plots for the SPD vertex.");
  Printf("\t 4. The QA plots for the tracks vertex.");

  Printf("\nFunction drawKineQA: Takes as an argument the file Protons.QA.Histograms.root with all the histograms about kinematic or not variables");
  Printf("\t 1. The phi vs eta correlation plots for primary and secondary (anti)protons.");
  Printf("\t 2. The Nclusters TPC vs eta correlation plots for primary and secondary (anti)protons.");
  Printf("\t 3. The Nclusters TPC vs phi correlation plots for primary and secondary (anti)protons.");

 Printf("\nFunction drawEfficiencies: Takes as arguments the filename produced by the drawEfficiency (see before - usually the name is Reconstruction-PID-Efficiency.root) and three booleans indicating whether the points for primaries, secondaries from weak decays and secondaries from hadronic interactions will be shown.It displays:");
  Printf("\t 1. Reconstruction efficiencies for primary and secondary (anti)protons as a function of either eta or y. The secondaries are categorized by the MC process as originating from a weak decay or from a hadronic interaction.");
  Printf("\t 2. Reconstruction efficiencies for primary and secondary (anti)protons as a function of Pt. The secondaries are categorized by the MC process as originating from a weak decay or from a hadronic interaction.");
  Printf("\t 3. Particle identifiction efficiencies for (anti)protons as a function of either eta or y and of Pt.");

 Printf("\nFunction compareEfficiencies: Takes as arguments the filenames produced by the drawEfficiency (see before - usually the name is Reconstruction-PID-Efficiency.root) and three booleans indicating whether the points for primaries, secondaries from weak decays and secondaries from hadronic interactions will be shown. The two files correspond to two different tracking methods (e.g. TPC, Hybrid, Global). It displays:");
  Printf("\t 1. Reconstruction efficiencies for primary and secondary (anti)protons as a function of either eta or y. The secondaries are categorized by the MC process as originating from a weak decay or from a hadronic interaction.");
  Printf("\t 2. Reconstruction efficiencies for primary and secondary (anti)protons as a function of Pt. The secondaries are categorized by the MC process as originating from a weak decay or from a hadronic interaction.");

  Printf("\nFunction drawCutParametersDistributions: Takes as an argument the third file created by the QA code (usually the name is Protons.QA.Histograms.<AnalysisMode>.root) and displays the distributions of accepted primaries (in blue) and secondaries (in orange) for the different cut parameters.");
  Printf("==================================================================\n\n\n");
}

//________________________________________//
void drawCutStatistics(TList *list,
		       const char* analysisType) {
  //Function to display the statistics from the cuts
  //The histogram shows the influence of each cut on the primary
  //and secondary (anti)protons
  const Int_t NQAHISTOSPERLIST = 27;

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
  const Int_t nx = 27;
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
			"TPC pid",
			"N_{points} (dE/dx)","",
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
    hPrimaryProtonsITS->SetBinContent(i,GetPercentage(gEntriesQAPrimaryProtonsAcceptedList[20+i],
						      gEntriesQAPrimaryProtonsRejectedList[20+i]));
    hSecondaryProtonsITS->SetBinContent(i,GetPercentage(gEntriesQASecondaryProtonsAcceptedList[20+i],
							gEntriesQASecondaryProtonsRejectedList[20+i]));
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
    hPrimaryAntiProtonsITS->SetBinContent(i,GetPercentage(gEntriesQAPrimaryAntiProtonsAcceptedList[20+i],
							  gEntriesQAPrimaryAntiProtonsRejectedList[20+i]));
    hSecondaryAntiProtonsITS->SetBinContent(i,GetPercentage(gEntriesQASecondaryAntiProtonsAcceptedList[20+i],
							    gEntriesQASecondaryAntiProtonsRejectedList[20+i]));
  }
  hPrimaryAntiProtonsITS->DrawCopy("EHISTSAME");
  hSecondaryAntiProtonsITS->DrawCopy("EHISTSAME");
  DrawMarker(6.5, 90, 20, 1.2, 4);
  t1->DrawLatex(7,88,"Primary #bar{p}");
  DrawMarker(6.5, 80, 22, 1.2, 2);
  t1->DrawLatex(7,78,"Secondary #bar{p}");

  c9->SaveAs("CutInfluence-ITS.gif");
  
  //Efficiency - Contamination plots
  DrawContamination(fQA2DList,analysisType);
  DrawCutEfficiency(fQA2DList,analysisType);
  DrawComposition(gHistYPtPDGProtonsPass,gHistYPtPDGAntiProtonsPass);
}

//________________________________________//
void DrawComposition(TH3F *gHistYPtPDGProtons,
		     TH3F *gHistYPtPDGAntiProtons) {
  //Function to display the composition of secondary (anti)protons
  //that survive the quality criteria
  Double_t nParticleCompositionProtonY[200], nParticleCompositionProtonPt[200];
  Double_t nParticleCompositionProtonYError[200], nParticleCompositionProtonPtError[200];
  Double_t nParticleCompositionAntiProtonY[200], nParticleCompositionAntiProtonPt[200];
  Double_t nParticleCompositionAntiProtonYError[200], nParticleCompositionAntiProtonPtError[200];
  Double_t gY[200], gPt[200];
  Double_t gYError[200], gPtError[200];
  for(Int_t iBins = 0; iBins < 200; iBins++) {
    nParticleCompositionProtonY[iBins] = 0;
    nParticleCompositionProtonPt[iBins] = 0;
    nParticleCompositionProtonYError[iBins] = 0;
    nParticleCompositionProtonPtError[iBins] = 0;
    nParticleCompositionAntiProtonY[iBins] = 0;
    nParticleCompositionAntiProtonPt[iBins] = 0;
    nParticleCompositionAntiProtonYError[iBins] = 0;
    nParticleCompositionAntiProtonPtError[iBins] = 0;
    gY[iBins] = 0;
    gPt[iBins] = 0;
    gYError[iBins] = 0;
    gPtError[iBins] = 0;
  }
  
  TGraphErrors *gParticleProtonY[14];
  TGraphErrors *gParticleProtonPt[14];
  TGraphErrors *gParticleAntiProtonY[14];
  TGraphErrors *gParticleAntiProtonPt[14];
  for(Int_t iParticle = 0; iParticle < 14; iParticle++) {
    GetComposition(iParticle,
		   gHistYPtPDGProtons,
		   nParticleCompositionProtonY,nParticleCompositionProtonYError,gY, gYError,
		   nParticleCompositionProtonPt, nParticleCompositionProtonPtError, gPt, gPtError);
    gParticleProtonY[iParticle] = new TGraphErrors(gHistYPtPDGProtons->GetNbinsX(),
						  gY,nParticleCompositionProtonY,gYError,nParticleCompositionProtonYError);
    gParticleProtonY[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleProtonY[iParticle]->SetMarkerSize(1.2);

    gParticleProtonPt[iParticle] = new TGraphErrors(gHistYPtPDGProtons->GetNbinsY(),
						   gPt,nParticleCompositionProtonPt,gPtError,nParticleCompositionProtonPtError);
    gParticleProtonPt[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleProtonPt[iParticle]->SetMarkerSize(1.2);

    GetComposition(iParticle,
		   gHistYPtPDGAntiProtons,
		   nParticleCompositionAntiProtonY,nParticleCompositionAntiProtonYError,gY, gYError, 
		   nParticleCompositionAntiProtonPt, nParticleCompositionAntiProtonPtError, gPt, gPtError);
    gParticleAntiProtonY[iParticle] = new TGraphErrors(gHistYPtPDGAntiProtons->GetNbinsX(),
						      gY,nParticleCompositionAntiProtonY,gYError,nParticleCompositionAntiProtonYError);
    gParticleAntiProtonY[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleAntiProtonY[iParticle]->SetMarkerSize(1.2);

    gParticleAntiProtonPt[iParticle] = new TGraphErrors(gHistYPtPDGAntiProtons->GetNbinsY(),
						       gPt,nParticleCompositionAntiProtonPt,gPtError,nParticleCompositionAntiProtonPtError);
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
  hEmptyY->GetXaxis()->SetTitle(gHistYPtPDGProtons->GetXaxis()->GetTitle());

  TCanvas *c12 = new TCanvas("c12",
			     "Composition of accepted secondaries vs y",
			     350,350,700,400);
  c12->Divide(2,1);
  c12->SetHighLightColor(10); c12->cd(1)->SetBottomMargin(0.15);
  c12->cd(1)->SetGridx(); c12->cd(1)->SetGridy();
  hEmptyY->SetTitle("Protons");
  hEmptyY->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    //if((iParticle == 0)||(iParticle == 2)||(iParticle == 5)||(iParticle == 6)||(iParticle == 8))
    if((iParticle == 2)||(iParticle == 6)||(iParticle == 8)) {
      gParticleProtonY[iParticle]->Draw("P");
      //if(iParticle < 5) {
      //DrawMarker(-1.1, 115-5*iParticle, 20+iParticle, 1.2, 1);
      //t1->DrawLatex(-1.0,113-5*iParticle,fParticleName[iParticle]);
    }
    /*else {
      DrawMarker(0.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*(iParticle-5),fParticleName[iParticle]);
      }*/
  }
  DrawMarker(0.0, 115, 22, 1.2, 1);
  t1->DrawLatex(0.1,113,fParticleName[2]);
  DrawMarker(0.0, 105, 26, 1.2, 1);
  t1->DrawLatex(0.1,103,fParticleName[6]);
  DrawMarker(0.0, 95, 28, 1.2, 1);
  t1->DrawLatex(0.1,93,fParticleName[8]);

  c12->SetHighLightColor(10); c12->cd(2)->SetBottomMargin(0.15);
  c12->cd(2)->SetGridx(); c12->cd(2)->SetGridy();
  hEmptyY->SetTitle("Antiprotons");
  hEmptyY->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    if((iParticle == 6)||(iParticle == 8))
      gParticleAntiProtonY[iParticle]->Draw("P");
    /*if(iParticle < 5) {
      DrawMarker(-1.1, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(-1.0,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(0.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*(iParticle-5),fParticleName[iParticle]);
      }*/
  }
  DrawMarker(0.0, 115, 26, 1.2, 1);
  t1->DrawLatex(0.1,113,fParticleName[6]);
  DrawMarker(0.0, 105, 28, 1.2, 1);
  t1->DrawLatex(0.1,103,fParticleName[8]);
  c12->SaveAs("SurvivedSecondaries-Composition-Rapidity.gif");

  TH2F *hEmptyPt = new TH2F("hEmptyPt","",100,0.0,4.0,100,0,120); 
  hEmptyPt->SetStats(kFALSE); 
  hEmptyPt->GetYaxis()->SetTitle("Particle composition [%]");
  hEmptyPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");

  TCanvas *c13 = new TCanvas("c13",
			     "Composition of accepted secondaries vs pT",
			     400,400,700,400);
  c13->Divide(2,1);
  c13->SetHighLightColor(10); c13->cd(1)->SetBottomMargin(0.15);
  c13->cd(1)->SetGridx(); c13->cd(1)->SetGridy();
  hEmptyPt->GetXaxis()->SetRangeUser(gParticleProtonPt[0]->GetXaxis()->GetXmin()-0.2,gParticleProtonPt[0]->GetXaxis()->GetXmax()+0.2);
  hEmptyPt->SetTitle("Protons");
  hEmptyPt->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    /*if(iParticle < 5) {
      DrawMarker(gParticleProtonPt[0]->GetXaxis()->GetXmin()+0.1, 
		 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(gParticleProtonPt[0]->GetXaxis()->GetXmin()+0.2,
		    113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(gParticleProtonPt[0]->GetXaxis()->GetXmax()*0.5, 
		 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(gParticleProtonPt[0]->GetXaxis()->GetXmax()*0.5+0.1,
		    113-5*(iParticle-5),fParticleName[iParticle]);
		    }*/
    if((iParticle == 2)||(iParticle == 6)||(iParticle == 8))
      gParticleProtonPt[iParticle]->Draw("P");
  }
  DrawMarker(0.5, 115, 22, 1.2, 1);
  t1->DrawLatex(0.6,113,fParticleName[2]);
  DrawMarker(0.5, 105, 26, 1.2, 1);
  t1->DrawLatex(0.6,103,fParticleName[6]);
  DrawMarker(0.5, 95, 28, 1.2, 1);
  t1->DrawLatex(0.6,93,fParticleName[8]);

  c13->SetHighLightColor(10); c13->cd(2)->SetBottomMargin(0.15);
  c13->cd(2)->SetGridx(); c13->cd(2)->SetGridy();
  hEmptyPt->SetTitle("Antiprotons");
  hEmptyPt->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    /*if(iParticle < 5) {
      DrawMarker(gParticleProtonPt[0]->GetXaxis()->GetXmin()+0.1, 
		 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(gParticleProtonPt[0]->GetXaxis()->GetXmin()+0.2,
		    113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(gParticleProtonPt[0]->GetXaxis()->GetXmax()*0.5, 
		 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(gParticleProtonPt[0]->GetXaxis()->GetXmax()*0.5+0.1,
		    113-5*(iParticle-5),fParticleName[iParticle]);
		    }*/
    if((iParticle == 6)||(iParticle == 8))
      gParticleAntiProtonPt[iParticle]->Draw("P");
  }
  DrawMarker(0.5, 115, 26, 1.2, 1);
  t1->DrawLatex(0.6,113,fParticleName[6]);
  DrawMarker(0.5, 105, 28, 1.2, 1);
  t1->DrawLatex(0.6,103,fParticleName[8]);
  c13->SaveAs("SurvivedSecondaries-Composition-Pt.gif");
}

//________________________________________________//
void DrawContamination(TList *inputList,
		       const char* analysisType) {
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
  hEmptyY->GetXaxis()->SetTitle(hPrimaryProtons->GetXaxis()->GetTitle());

  TCanvas *c7 = new TCanvas("c7","(Anti)Proton contamination vs y",
			    150,150,700,400);
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
  t1->DrawLatex(0.1,41,"Secondaries");

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
  c8->SetHighLightColor(10); c8->Divide(2,1);

  c8->cd(1)->SetBottomMargin(0.15); 
  c8->cd(1)->SetLeftMargin(0.15); 
  c8->cd(1)->SetGridx(); c8->cd(1)->SetGridy();
  hEmptyPt->SetTitle("Protons");
  hEmptyPt->GetXaxis()->SetRangeUser(gPtPrimaryProtonsPercentage->GetXaxis()->GetXmin()-0.2,gPtPrimaryProtonsPercentage->GetXaxis()->GetXmax()+0.2);
  hEmptyPt->DrawCopy();
  gPtPrimaryProtonsPercentage->DrawCopy("ESAME");
  gPtSecondaryProtonsPercentage->DrawCopy("ESAME");

  DrawMarker(0.5, 55, kFullCircle, 1.2, 1);
  t1->DrawLatex(0.6,53,"Primaries");
  DrawMarker(0.5, 45, kOpenCircle, 1.2, 1);
  t1->DrawLatex(0.6,41,"Secondaries");

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

  TString outputFileName = "Contamination."; 
  outputFileName += analysisType; outputFileName += ".root";
  TFile *fout = TFile::Open(outputFileName.Data(),"recreate");
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
void DrawCutEfficiency(TList *inputList,
		       const char* analysisType) {
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
  hEmptyY->GetXaxis()->SetTitle(hPrimaryESDProtons->GetXaxis()->GetTitle());

  TCanvas *c10 = new TCanvas("c10","(Anti)Proton cut efficiency vs y",
			    250,250,700,400);
  c10->SetHighLightColor(10); c10->Divide(2,1);

  c10->cd(1)->SetBottomMargin(0.15); 
  c10->cd(1)->SetLeftMargin(0.15); 
  c10->cd(1)->SetGridx(); c10->cd(1)->SetGridy();
  hEmptyY->SetTitle("Protons");
  hEmptyY->GetXaxis()->SetRangeUser(gYPrimaryESDAntiProtons->GetXaxis()->GetXmin()-0.2,gYPrimaryESDAntiProtons->GetXaxis()->GetXmax()+0.2);
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
  c11->SetHighLightColor(10); c11->Divide(2,1);

  c11->cd(1)->SetBottomMargin(0.15); 
  c11->cd(1)->SetLeftMargin(0.15); 
  c11->cd(1)->SetGridx(); c11->cd(1)->SetGridy();
  hEmptyPt->SetTitle("Protons");
  hEmptyPt->GetXaxis()->SetRangeUser(gPtPrimaryESDAntiProtons->GetXaxis()->GetXmin()-0.2,gPtPrimaryESDAntiProtons->GetXaxis()->GetXmax()+0.2);
  hEmptyPt->DrawCopy();
  gPtPrimaryESDProtons->DrawCopy("ESAME");

  c11->cd(2)->SetBottomMargin(0.15); 
  c11->cd(2)->SetLeftMargin(0.15); 
  c11->cd(2)->SetGridx(); c11->cd(2)->SetGridy();
  hEmptyPt->SetTitle("Antirotons");
  hEmptyPt->DrawCopy();
  gPtPrimaryESDAntiProtons->DrawCopy("ESAME");

  c11->SaveAs("CutEfficiency-Protons-Pt.gif");

  TString outputFileName = "CutEfficiency.";
  outputFileName += analysisType; outputFileName += ".root";
  TFile *fout = TFile::Open(outputFileName.Data(),"recreate");
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
    cout<<"Position: "<<i+1<<" - Histogram: "<<gHist->GetName()<<
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
  Double_t nParticleCompositionProtonY[200], nParticleCompositionProtonPt[200];
  Double_t nParticleCompositionProtonYError[200], nParticleCompositionProtonPtError[200];
  Double_t nParticleCompositionAntiProtonY[200], nParticleCompositionAntiProtonPt[200];
  Double_t nParticleCompositionAntiProtonYError[200], nParticleCompositionAntiProtonPtError[200];
  Double_t gY[200], gPt[200];
  Double_t gYError[200], gPtError[200];
  for(Int_t iBins = 0; iBins < 200; iBins++) {
    nParticleCompositionProtonY[iBins] = 0;
    nParticleCompositionProtonPt[iBins] = 0;
    nParticleCompositionProtonYError[iBins] = 0;
    nParticleCompositionProtonPtError[iBins] = 0;
    nParticleCompositionAntiProtonY[iBins] = 0;
    nParticleCompositionAntiProtonPt[iBins] = 0;
    nParticleCompositionAntiProtonYError[iBins] = 0;
    nParticleCompositionAntiProtonPtError[iBins] = 0;
    gY[iBins] = 0;
    gPt[iBins] = 0;
    gYError[iBins] = 0;
    gPtError[iBins] = 0;
  }
  
  TGraphErrors *gParticleProtonY[14];
  TGraphErrors *gParticleProtonPt[14];
  TGraphErrors *gParticleAntiProtonY[14];
  TGraphErrors *gParticleAntiProtonPt[14];
  for(Int_t iParticle = 0; iParticle < 14; iParticle++) {
    GetComposition(iParticle,
		   gHistYPtPDGProtons,
		   nParticleCompositionProtonY,
		   nParticleCompositionProtonYError, gY, gYError, 
		   nParticleCompositionProtonPt, 
		   nParticleCompositionProtonPtError, gPt, gPtError);
    gParticleProtonY[iParticle] = new TGraphErrors(gHistYPtPDGProtons->GetNbinsX(),
						   gY,nParticleCompositionProtonY,
						   gYError,nParticleCompositionProtonYError);
    gParticleProtonY[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleProtonY[iParticle]->SetMarkerSize(1.2);

    gParticleProtonPt[iParticle] = new TGraphErrors(gHistYPtPDGProtons->GetNbinsY(),
						    gPt,nParticleCompositionProtonPt,
						    gPtError,nParticleCompositionProtonPtError);
    gParticleProtonPt[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleProtonPt[iParticle]->SetMarkerSize(1.2);

    GetComposition(iParticle,
		   gHistYPtPDGAntiProtons,
		   nParticleCompositionAntiProtonY,
		   nParticleCompositionAntiProtonYError, gY, gYError, 
		   nParticleCompositionAntiProtonPt, 
		   nParticleCompositionAntiProtonPtError, gPt, gPtError);
    gParticleAntiProtonY[iParticle] = new TGraphErrors(gHistYPtPDGAntiProtons->GetNbinsX(),
						       gY,nParticleCompositionAntiProtonY,
						       gYError,nParticleCompositionAntiProtonYError);
    gParticleAntiProtonY[iParticle]->SetMarkerStyle(iParticle+20);
    gParticleAntiProtonY[iParticle]->SetMarkerSize(1.2);

    gParticleAntiProtonPt[iParticle] = new TGraphErrors(gHistYPtPDGAntiProtons->GetNbinsY(),
							gPt,nParticleCompositionAntiProtonPt,
							gPtError,nParticleCompositionAntiProtonPtError);
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
  hEmptyY->GetXaxis()->SetTitle(gHistYPtPDGProtons->GetXaxis()->GetTitle());

  TCanvas *c3 = new TCanvas("c3","MC secondary composition vs y - Protons",
			    450,450,700,400);
  c3->SetHighLightColor(10); c3->SetBottomMargin(0.15);
  c3->SetGridx(); c3->SetGridy();
  hEmptyY->GetXaxis()->SetRangeUser(gParticleProtonY[0]->GetXaxis()->GetXmin()-0.2,gParticleProtonY[0]->GetXaxis()->GetXmax()+0.2);
  hEmptyY->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    //if((iParticle == 0)||(iParticle == 2)||(iParticle == 5)||(iParticle == 6)||(iParticle == 8))
    if((iParticle == 0)||(iParticle == 2)||(iParticle == 6)||(iParticle == 8))
      gParticleProtonY[iParticle]->Draw("P");
    /*if(iParticle < 5) {
      DrawMarker(-1.1, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(-1.0,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(0.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*(iParticle-5),fParticleName[iParticle]);
      }*/
  }
  DrawMarker(0.0, 115, 20, 1.2, 1);
  t1->DrawLatex(0.1,113,fParticleName[0]);
  DrawMarker(0.0, 108, 22, 1.2, 1);
  t1->DrawLatex(0.1,106,fParticleName[2]);
  DrawMarker(0.0, 101, 26, 1.2, 1);
  t1->DrawLatex(0.1,99,fParticleName[6]);
  DrawMarker(0.0, 94, 28, 1.2, 1);
  t1->DrawLatex(0.1,92,fParticleName[8]);

  TCanvas *c5 = new TCanvas("c5","MC secondary composition vs y - antiProtons",
			    500,500,700,400);
  c5->SetHighLightColor(10); c5->SetBottomMargin(0.15);
  c5->SetGridx(); c5->SetGridy();
  hEmptyY->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    if((iParticle == 0)||(iParticle == 6)||(iParticle == 8))
      gParticleAntiProtonY[iParticle]->Draw("P");
    /*if(iParticle < 5) {
      DrawMarker(-1.1, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(-1.0,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(0.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*(iParticle-5),fParticleName[iParticle]);
      }*/
  }
  DrawMarker(0.0, 115, 20, 1.2, 1);
  t1->DrawLatex(0.1,113,fParticleName[0]);
  DrawMarker(0.0, 108, 26, 1.2, 1);
  t1->DrawLatex(0.1,106,fParticleName[6]);
  DrawMarker(0.0, 101, 28, 1.2, 1);
  t1->DrawLatex(0.1,99,fParticleName[8]);

  TH2F *hEmptyPt = new TH2F("hEmptyPt","",100,0.0,4.0,100,0,120); 
  hEmptyPt->SetStats(kFALSE); 
  hEmptyPt->GetYaxis()->SetTitle("Particle composition [%]");
  hEmptyPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");

  TCanvas *c4 = new TCanvas("c4","MC secondary composition vs pT - Protons",
			    550,550,700,400);
  c4->SetHighLightColor(10); c4->SetBottomMargin(0.15);
  c4->SetGridx(); c4->SetGridy();
  hEmptyPt->GetXaxis()->SetRangeUser(gParticleProtonPt[0]->GetXaxis()->GetXmin()-0.2,gParticleProtonPt[0]->GetXaxis()->GetXmax()+0.2);
  hEmptyPt->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    /*if(iParticle < 5) {
      DrawMarker(0.2, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(2.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(2.3,113-5*(iParticle-5),fParticleName[iParticle]);
      }*/
    if((iParticle == 0)||(iParticle == 2)||(iParticle == 6)||(iParticle == 8))
      gParticleProtonPt[iParticle]->Draw("P");
  }
  DrawMarker(0.5, 115, 20, 1.2, 1);
  t1->DrawLatex(0.6,113,fParticleName[0]);
  DrawMarker(0.5, 108, 22, 1.2, 1);
  t1->DrawLatex(0.6,106,fParticleName[2]);
  DrawMarker(0.5, 101, 26, 1.2, 1);
  t1->DrawLatex(0.6,99,fParticleName[6]);
  DrawMarker(0.5, 94, 28, 1.2, 1);
  t1->DrawLatex(0.6,92,fParticleName[8]);

  TCanvas *c6 = new TCanvas("c6",
			    "MC secondary composition vs pT - AntiProtons",
			    600,600,700,400);
  c6->SetHighLightColor(10); c6->SetBottomMargin(0.15);
  c6->SetGridx(); c6->SetGridy();
  hEmptyPt->DrawCopy();
  for(Int_t iParticle = 0; iParticle < 10; iParticle++) {
    /*if(iParticle < 5) {
      DrawMarker(0.2, 115-5*iParticle, 20+iParticle, 1.2, 1);
      t1->DrawLatex(0.3,113-5*iParticle,fParticleName[iParticle]);
    }
    else {
      DrawMarker(2.2, 115-5*(iParticle-5), 20+iParticle, 1.2, 1);
      t1->DrawLatex(2.3,113-5*(iParticle-5),fParticleName[iParticle]);
      }*/
    if((iParticle == 0)||(iParticle == 6)||(iParticle == 8))
      gParticleAntiProtonPt[iParticle]->Draw("P");
  }
  DrawMarker(0.5, 115, 20, 1.2, 1);
  t1->DrawLatex(0.6,113,fParticleName[0]);
  DrawMarker(0.5, 108, 26, 1.2, 1);
  t1->DrawLatex(0.6,106,fParticleName[6]);
  DrawMarker(0.5, 101, 28, 1.2, 1);
  t1->DrawLatex(0.6,99,fParticleName[8]);
}

//________________________________________//
void GetComposition(Int_t iSpecies,
		    TH3F *gHist, 
		    Double_t *nParticleCompositionY,
		    Double_t *nParticleCompositionYError,
		    Double_t *gY, Double_t *gYError,
		    Double_t *nParticleCompositionPt,
		    Double_t *nParticleCompositionPtError,
		    Double_t *gPt, Double_t *gPtError) {
  //Returns the pT and y dependence of the MC composition
  Double_t ymin = gHist->GetXaxis()->GetXmin();
  Double_t ymax = gHist->GetXaxis()->GetXmax();
  Double_t nybins = gHist->GetNbinsX();
  Double_t ptmin = gHist->GetYaxis()->GetXmin();
  Double_t ptmax = gHist->GetYaxis()->GetXmax();
  Double_t nptbins = gHist->GetNbinsY();
  Double_t nTotalY[200], nTotalPt[200];
  for(Int_t iBins = 0; iBins < 200; iBins++) {
    nParticleCompositionY[iBins] = 0;
    nParticleCompositionPt[iBins] = 0;
    nParticleCompositionYError[iBins] = 0;
    nParticleCompositionPtError[iBins] = 0;
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
      //nCompositionY[iXbins-1] += gHist->GetBinContent(iXbins,iYbins,iSpecies+1);
      //if(nParticleCompositionY[iXbins-1] == 0) 
      //nParticleCompositionY[iXbins-1] = -10.0;
    }//pt loop
    if((nParticleCompositionY[iXbins-1] <= 100.)&&(nTotalY[iXbins-1] != 0))
      nParticleCompositionYError[iXbins-1] = TMath::Sqrt(nParticleCompositionY[iXbins-1]*(100. - nParticleCompositionY[iXbins-1])/nTotalY[iXbins-1]);
    gY[iXbins-1] = ymin + (iXbins-1)*(ymax - ymin)/nybins + 0.5*(ymax - ymin)/nybins;
    gYError[iXbins-1] = 0.5*(ymax - ymin)/nybins;
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
    if((nParticleCompositionPt[iYbins-1] <= 100.)&&(nTotalPt[iYbins-1] != 0))
      nParticleCompositionPtError[iYbins-1] = TMath::Sqrt(nParticleCompositionPt[iYbins-1]*(100. - nParticleCompositionPt[iYbins-1])/nTotalPt[iYbins-1]);
    gPt[iYbins-1] = ptmin + (iYbins-1)*(ptmax - ptmin)/nptbins + 0.5*(ptmax - ptmin)/nptbins;
    gPtError[iYbins-1] = 0.5*(ptmax - ptmin)/nptbins;
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
void SetError(TH1 *hEff, TH1 *hGen) {
  for(Int_t iBin = 1; iBin <= hEff->GetNbinsX(); iBin++) {
    Double_t error = 0.0;
    if((hEff->GetBinContent(iBin) <= 1  .)&&(hGen->GetBinContent(iBin) != 0))
      error = TMath::Sqrt(hEff->GetBinContent(iBin)*(1  . - hEff->GetBinContent(iBin))/hGen->GetBinContent(iBin));
    hEff->SetBinError(iBin,error);
  }
}

//________________________________________//
void drawVertexQA(TList *list) {
  //Function to display the vertex QA plots
  TH1I *gHistMCPrimaryMultiplicity = (TH1I *)list->At(0);
  //TPC
  TH1I *gHistMCPrimaryMultiplicityTPC = (TH1I *)list->At(1);
  TH2F *gHistTPCESDVx = (TH2F *)list->At(2);
  TH2F *gHistTPCESDVy = (TH2F *)list->At(3);
  TH2F *gHistTPCESDVz = (TH2F *)list->At(4);
  TH1F *gHistTPCDiffVx = (TH1F *)list->At(5);
  TH1F *gHistTPCDiffVy = (TH1F *)list->At(6);
  TH1F *gHistTPCDiffVz = (TH1F *)list->At(7);
  TH1F *gHistTPCResolutionVx = (TH1F *)list->At(8);
  TH1F *gHistTPCResolutionVy = (TH1F *)list->At(9);
  TH1F *gHistTPCResolutionVz = (TH1F *)list->At(10);
  //SPD
  TH1I *gHistMCPrimaryMultiplicitySPD = (TH1I *)list->At(11);
  TH2F *gHistSPDESDVx = (TH2F *)list->At(12);
  TH2F *gHistSPDESDVy = (TH2F *)list->At(13);
  TH2F *gHistSPDESDVz = (TH2F *)list->At(14);
  TH1F *gHistSPDDiffVx = (TH1F *)list->At(15);
  TH1F *gHistSPDDiffVy = (TH1F *)list->At(16);
  TH1F *gHistSPDDiffVz = (TH1F *)list->At(17);
  TH1F *gHistSPDResolutionVx = (TH1F *)list->At(18);
  TH1F *gHistSPDResolutionVy = (TH1F *)list->At(19);
  TH1F *gHistSPDResolutionVz = (TH1F *)list->At(20);
  //Tracks
  TH1I *gHistMCPrimaryMultiplicityTracks = (TH1I *)list->At(21);
  TH2F *gHistTracksESDVx = (TH2F *)list->At(22);
  TH2F *gHistTracksESDVy = (TH2F *)list->At(23);
  TH2F *gHistTracksESDVz = (TH2F *)list->At(24);
  TH1F *gHistTracksDiffVx = (TH1F *)list->At(25);
  TH1F *gHistTracksDiffVy = (TH1F *)list->At(26);
  TH1F *gHistTracksDiffVz = (TH1F *)list->At(27);
  TH1F *gHistTracksResolutionVx = (TH1F *)list->At(28);
  TH1F *gHistTracksResolutionVy = (TH1F *)list->At(29);
  TH1F *gHistTracksResolutionVz = (TH1F *)list->At(30);

  TCanvas *c17 = new TCanvas("c17",
			     "Vertex efficiency",
			     300,0,900,400);
  c17->SetHighLightColor(10); c17->Divide(3,1);
  c17->cd(1)->SetLeftMargin(0.15); c17->cd(1)->SetBottomMargin(0.15);  
  c17->cd(1)->SetRightMargin(0.2);
  c17->cd(1)->SetGridx(); c17->cd(1)->SetGridy();
  gHistMCPrimaryMultiplicityTPC->SetTitle("TPC vertex");
  gHistMCPrimaryMultiplicityTPC->Divide(gHistMCPrimaryMultiplicity);
  SetError(gHistMCPrimaryMultiplicityTPC,gHistMCPrimaryMultiplicity);
  gHistMCPrimaryMultiplicityTPC->SetMaximum(110.);
  gHistMCPrimaryMultiplicityTPC->Scale(100.);
  gHistMCPrimaryMultiplicityTPC->SetMarkerStyle(20);
  gHistMCPrimaryMultiplicityTPC->SetMarkerColor(1);
  gHistMCPrimaryMultiplicityTPC->GetYaxis()->SetTitle("#epsilon [%]");
  gHistMCPrimaryMultiplicityTPC->SetStats(kFALSE);
  gHistMCPrimaryMultiplicityTPC->Draw("E");
  c17->cd(2)->SetLeftMargin(0.15); c17->cd(2)->SetBottomMargin(0.15);  
  c17->cd(2)->SetRightMargin(0.2);
  c17->cd(2)->SetGridx(); c17->cd(2)->SetGridy();
  gHistMCPrimaryMultiplicitySPD->SetTitle("SPD vertex");
  gHistMCPrimaryMultiplicitySPD->Divide(gHistMCPrimaryMultiplicity);
  SetError(gHistMCPrimaryMultiplicitySPD,gHistMCPrimaryMultiplicity);
  gHistMCPrimaryMultiplicitySPD->SetMaximum(110.);
  gHistMCPrimaryMultiplicitySPD->Scale(100.);
  gHistMCPrimaryMultiplicitySPD->SetMarkerStyle(20);
  gHistMCPrimaryMultiplicitySPD->SetMarkerColor(1);
  gHistMCPrimaryMultiplicitySPD->GetYaxis()->SetTitle("#epsilon [%]");
  gHistMCPrimaryMultiplicitySPD->SetStats(kFALSE);
  gHistMCPrimaryMultiplicitySPD->Draw("E");
  c17->cd(3)->SetLeftMargin(0.15); c17->cd(3)->SetBottomMargin(0.15);  
  c17->cd(3)->SetRightMargin(0.2);
  c17->cd(3)->SetGridx(); c17->cd(3)->SetGridy();
  gHistMCPrimaryMultiplicityTracks->SetTitle("Vertex from tracks");
  gHistMCPrimaryMultiplicityTracks->Divide(gHistMCPrimaryMultiplicity);
  SetError(gHistMCPrimaryMultiplicityTracks,gHistMCPrimaryMultiplicity);
  gHistMCPrimaryMultiplicityTracks->SetMaximum(110.);
  gHistMCPrimaryMultiplicityTracks->Scale(100.);
  gHistMCPrimaryMultiplicityTracks->SetMarkerStyle(20);
  gHistMCPrimaryMultiplicityTracks->SetMarkerColor(1);
  gHistMCPrimaryMultiplicityTracks->GetYaxis()->SetTitle("#epsilon [%]");
  gHistMCPrimaryMultiplicityTracks->SetStats(kFALSE);
  gHistMCPrimaryMultiplicityTracks->Draw("E");
  c17->SaveAs("VertexEfficiency.gif");

  //TPC vertex
  TCanvas *c18 = new TCanvas("c18",
			     "TPC vertex",
			     350,50,700,700);
  c18->SetHighLightColor(10); c18->Divide(3,3);
  c18->cd(1)->SetLeftMargin(0.15); c18->cd(1)->SetBottomMargin(0.15);  
  c18->cd(1)->SetRightMargin(0.2); c18->cd(1)->SetLogy();
  gHistTPCESDVx->Draw("col");
  c18->cd(2)->SetLeftMargin(0.15); c18->cd(2)->SetBottomMargin(0.15);  
  c18->cd(2)->SetRightMargin(0.2); c18->cd(2)->SetLogy();
  gHistTPCESDVy->Draw("col");
  c18->cd(3)->SetLeftMargin(0.15); c18->cd(3)->SetBottomMargin(0.15);  
  c18->cd(3)->SetRightMargin(0.2); c18->cd(3)->SetLogy();
  gHistTPCESDVz->Draw("col");
  c18->cd(4)->SetLeftMargin(0.15); c18->cd(4)->SetBottomMargin(0.15);  
  c18->cd(4)->SetRightMargin(0.2); c18->cd(4)->SetLogy();
  gHistTPCDiffVx->Draw();
  c18->cd(5)->SetLeftMargin(0.15); c18->cd(5)->SetBottomMargin(0.15);  
  c18->cd(5)->SetRightMargin(0.2); c18->cd(5)->SetLogy();
  gHistTPCDiffVy->Draw();
  c18->cd(6)->SetLeftMargin(0.15); c18->cd(6)->SetBottomMargin(0.15);  
  c18->cd(6)->SetRightMargin(0.2); c18->cd(6)->SetLogy();
  gHistTPCDiffVz->Draw();
  c18->cd(7)->SetLeftMargin(0.15); c18->cd(7)->SetBottomMargin(0.15);  
  c18->cd(7)->SetRightMargin(0.2); c18->cd(7)->SetLogy();
  gHistTPCResolutionVx->Draw();
  c18->cd(8)->SetLeftMargin(0.15); c18->cd(8)->SetBottomMargin(0.15);  
  c18->cd(8)->SetRightMargin(0.2); c18->cd(8)->SetLogy();
  gHistTPCResolutionVy->Draw();
  c18->cd(9)->SetLeftMargin(0.15); c18->cd(9)->SetBottomMargin(0.15);  
  c18->cd(9)->SetRightMargin(0.2); c18->cd(9)->SetLogy();
  gHistTPCResolutionVz->Draw();
  c18->SaveAs("VertexTPC.gif");

  //SPD vertex
  TCanvas *c19 = new TCanvas("c19",
			     "SPD vertex",
			     400,100,700,700);
  c19->SetHighLightColor(10); c19->Divide(3,3);
  c19->cd(1)->SetLeftMargin(0.15); c19->cd(1)->SetBottomMargin(0.15);  
  c19->cd(1)->SetRightMargin(0.2); c19->cd(1)->SetLogy();
  gHistSPDESDVx->Draw("col");
  c19->cd(2)->SetLeftMargin(0.15); c19->cd(2)->SetBottomMargin(0.15);  
  c19->cd(2)->SetRightMargin(0.2); c19->cd(2)->SetLogy();
  gHistSPDESDVy->Draw("col");
  c19->cd(3)->SetLeftMargin(0.15); c19->cd(3)->SetBottomMargin(0.15);  
  c19->cd(3)->SetRightMargin(0.2); c19->cd(3)->SetLogy();
  gHistSPDESDVz->Draw("col");
  c19->cd(4)->SetLeftMargin(0.15); c19->cd(4)->SetBottomMargin(0.15);  
  c19->cd(4)->SetRightMargin(0.2); c19->cd(4)->SetLogy();
  gHistSPDDiffVx->Draw();
  c19->cd(5)->SetLeftMargin(0.15); c19->cd(5)->SetBottomMargin(0.15);  
  c19->cd(5)->SetRightMargin(0.2); c19->cd(5)->SetLogy();
  gHistSPDDiffVy->Draw();
  c19->cd(6)->SetLeftMargin(0.15); c19->cd(6)->SetBottomMargin(0.15);  
  c19->cd(6)->SetRightMargin(0.2); c19->cd(6)->SetLogy();
  gHistSPDDiffVz->Draw();
  c19->cd(7)->SetLeftMargin(0.15); c19->cd(7)->SetBottomMargin(0.15);  
  c19->cd(7)->SetRightMargin(0.2); c19->cd(7)->SetLogy();
  gHistSPDResolutionVx->Draw();
  c19->cd(8)->SetLeftMargin(0.15); c19->cd(8)->SetBottomMargin(0.15);  
  c19->cd(8)->SetRightMargin(0.2); c19->cd(8)->SetLogy();
  gHistSPDResolutionVy->Draw();
  c19->cd(9)->SetLeftMargin(0.15); c19->cd(9)->SetBottomMargin(0.15);  
  c19->cd(9)->SetRightMargin(0.2); c19->cd(9)->SetLogy();
  gHistSPDResolutionVz->Draw();
  c19->SaveAs("VertexSPD.gif");

  //Tracks vertex
  TCanvas *c20 = new TCanvas("c20",
			     "Tracks vertex",
			     450,150,700,700);
  c20->SetHighLightColor(10); c20->Divide(3,3);
  c20->cd(1)->SetLeftMargin(0.15); c20->cd(1)->SetBottomMargin(0.15);  
  c20->cd(1)->SetRightMargin(0.2); c20->cd(1)->SetLogy();
  gHistTracksESDVx->Draw("col");
  c20->cd(2)->SetLeftMargin(0.15); c20->cd(2)->SetBottomMargin(0.15);  
  c20->cd(2)->SetRightMargin(0.2); c20->cd(2)->SetLogy();
  gHistTracksESDVy->Draw("col");
  c20->cd(3)->SetLeftMargin(0.15); c20->cd(3)->SetBottomMargin(0.15);  
  c20->cd(3)->SetRightMargin(0.2); c20->cd(3)->SetLogy();
  gHistTracksESDVz->Draw("col");
  c20->cd(4)->SetLeftMargin(0.15); c20->cd(4)->SetBottomMargin(0.15);  
  c20->cd(4)->SetRightMargin(0.2); c20->cd(4)->SetLogy();
  gHistTracksDiffVx->Draw();
  c20->cd(5)->SetLeftMargin(0.15); c20->cd(5)->SetBottomMargin(0.15);  
  c20->cd(5)->SetRightMargin(0.2); c20->cd(5)->SetLogy();
  gHistTracksDiffVy->Draw();
  c20->cd(6)->SetLeftMargin(0.15); c20->cd(6)->SetBottomMargin(0.15);  
  c20->cd(6)->SetRightMargin(0.2); c20->cd(6)->SetLogy();
  gHistTracksDiffVz->Draw();
  c20->cd(7)->SetLeftMargin(0.15); c20->cd(7)->SetBottomMargin(0.15);  
  c20->cd(7)->SetRightMargin(0.2); c20->cd(7)->SetLogy();
  gHistTracksResolutionVx->Draw();
  c20->cd(8)->SetLeftMargin(0.15); c20->cd(8)->SetBottomMargin(0.15);  
  c20->cd(8)->SetRightMargin(0.2); c20->cd(8)->SetLogy();
  gHistTracksResolutionVy->Draw();
  c20->cd(9)->SetLeftMargin(0.15); c20->cd(9)->SetBottomMargin(0.15);  
  c20->cd(9)->SetRightMargin(0.2); c20->cd(9)->SetLogy();
  gHistTracksResolutionVz->Draw();
  c20->SaveAs("VertexTracks.gif");
}

//________________________________________//
void drawEfficiency(TList *list,
		    const char* analysisType) {
  //Function to display the reconstruction and PID efficiencies
  //for protons and antiprotons vs y and pT

  TH2F *hEmpty = new TH2F("hEmptyReconstructionEfficiency","",
			   100,-1.2,3.5,100,-10.0,130); 
  hEmpty->SetStats(kFALSE); 
  hEmpty->GetYaxis()->SetTitle("#epsilon [%]");
  hEmpty->GetYaxis()->SetTitleOffset(1.3);

  //Reconstruction efficiency
  TH2D *gHistPrimariesMCYPtProtons = (TH2D *)list->At(0);
  TH2D *gHistPrimariesMCYPtAntiProtons = (TH2D *)list->At(1);
  TH2D *gHistMCYPtProtonsFromWeak = (TH2D *)list->At(2);
  TH2D *gHistMCYPtAntiProtonsFromWeak = (TH2D *)list->At(3);
  TH2D *gHistMCYPtProtonsFromHadronic = (TH2D *)list->At(4);
  TH2D *gHistMCYPtAntiProtonsFromHadronic = (TH2D *)list->At(5);
  TH2D *gHistPrimariesESDYPtProtons = (TH2D *)list->At(6);
  TH2D *gHistPrimariesESDYPtAntiProtons = (TH2D *)list->At(7);
  TH2D *gHistESDYPtProtonsFromWeak = (TH2D *)list->At(8);
  TH2D *gHistESDYPtAntiProtonsFromWeak = (TH2D *)list->At(9);
  TH2D *gHistESDYPtProtonsFromHadronic = (TH2D *)list->At(10);
  TH2D *gHistESDYPtAntiProtonsFromHadronic = (TH2D *)list->At(11);

  //rapidity dependence
  TCanvas *c14 = new TCanvas("c14",
			     "(Anti)Proton reconstruction efficiency vs y",
			     650,650,700,400);
  c14->SetHighLightColor(10); c14->Divide(2,1);

  //Primary Protons
  TH1D *gYPrimariesESDProtons = (TH1D *)gHistPrimariesESDYPtProtons->ProjectionX("gYPrimariesESDProtons",0,gHistPrimariesESDYPtProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYPrimariesMCProtons = (TH1D *)gHistPrimariesMCYPtProtons->ProjectionX("gYPrimariesMCProtons",0,gHistPrimariesMCYPtProtons->GetXaxis()->GetNbins(),"e");
  gYPrimariesESDProtons->Divide(gYPrimariesMCProtons);
  SetError(gYPrimariesESDProtons,gYPrimariesMCProtons);
  gYPrimariesESDProtons->Scale(100.);
  gYPrimariesESDProtons->SetMarkerStyle(kFullCircle);

  //Primary AntiProtons
  TH1D *gYPrimariesESDAntiProtons = (TH1D *)gHistPrimariesESDYPtAntiProtons->ProjectionX("gYPrimariesESDAntiProtons",0,gHistPrimariesESDYPtAntiProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYPrimariesMCAntiProtons = (TH1D *)gHistPrimariesMCYPtAntiProtons->ProjectionX("gYPrimariesMCAntiProtons",0,gHistPrimariesMCYPtProtons->GetXaxis()->GetNbins(),"e");
  gYPrimariesESDAntiProtons->Divide(gYPrimariesMCAntiProtons);
  SetError(gYPrimariesESDAntiProtons,gYPrimariesMCAntiProtons);
  gYPrimariesESDAntiProtons->Scale(100.);
  gYPrimariesESDAntiProtons->SetMarkerStyle(kFullCircle);

  //Protons from weak decays
  TH1D *gYESDProtonsFromWeak = (TH1D *)gHistESDYPtProtonsFromWeak->ProjectionX("gYESDProtonsFromWeak",0,gHistESDYPtProtonsFromWeak->GetXaxis()->GetNbins(),"e");
  TH1D *gYMCProtonsFromWeak = (TH1D *)gHistMCYPtProtonsFromWeak->ProjectionX("gYMCProtonsFromWeak",0,gHistMCYPtProtonsFromWeak->GetXaxis()->GetNbins(),"e");
  gYESDProtonsFromWeak->Divide(gYMCProtonsFromWeak);
  SetError(gYESDProtonsFromWeak,gYMCProtonsFromWeak);
  gYESDProtonsFromWeak->Scale(100.);
  gYESDProtonsFromWeak->SetMarkerStyle(21);
  gYESDProtonsFromWeak->SetMarkerColor(2);

  //AntiProtons from weak decays
  TH1D *gYESDAntiProtonsFromWeak = (TH1D *)gHistESDYPtAntiProtonsFromWeak->ProjectionX("gYESDAntiProtonsFromWeak",0,gHistESDYPtAntiProtonsFromWeak->GetXaxis()->GetNbins(),"e");
  TH1D *gYMCAntiProtonsFromWeak = (TH1D *)gHistMCYPtAntiProtonsFromWeak->ProjectionX("gYMCAntiProtonsFromWeak",0,gHistMCYPtProtonsFromWeak->GetXaxis()->GetNbins(),"e");
  gYESDAntiProtonsFromWeak->Divide(gYMCAntiProtonsFromWeak);
  SetError(gYESDAntiProtonsFromWeak,gYMCAntiProtonsFromWeak);
  gYESDAntiProtonsFromWeak->Scale(100.);
  gYESDAntiProtonsFromWeak->SetMarkerStyle(21);
  gYESDAntiProtonsFromWeak->SetMarkerColor(2);

  //Protons from hadronic interactions
  TH1D *gYESDProtonsFromHadronic = (TH1D *)gHistESDYPtProtonsFromHadronic->ProjectionX("gYESDProtonsFromHadronic",0,gHistESDYPtProtonsFromHadronic->GetXaxis()->GetNbins(),"e");
  TH1D *gYMCProtonsFromHadronic = (TH1D *)gHistMCYPtProtonsFromHadronic->ProjectionX("gYMCProtonsFromHadronic",0,gHistMCYPtProtonsFromHadronic->GetXaxis()->GetNbins(),"e");
  gYESDProtonsFromHadronic->Divide(gYMCProtonsFromHadronic);
  SetError(gYESDProtonsFromHadronic,gYMCProtonsFromHadronic);
  gYESDProtonsFromHadronic->Scale(100.);
  gYESDProtonsFromHadronic->SetMarkerStyle(22);
  gYESDProtonsFromHadronic->SetMarkerColor(3);

  //AntiProtons from hadronic interactions
  TH1D *gYESDAntiProtonsFromHadronic = (TH1D *)gHistESDYPtAntiProtonsFromHadronic->ProjectionX("gYESDAntiProtonsFromHadronic",0,gHistESDYPtAntiProtonsFromHadronic->GetXaxis()->GetNbins(),"e");
  TH1D *gYMCAntiProtonsFromHadronic = (TH1D *)gHistMCYPtAntiProtonsFromHadronic->ProjectionX("gYMCAntiProtonsFromHadronic",0,gHistMCYPtProtonsFromHadronic->GetXaxis()->GetNbins(),"e");
  gYESDAntiProtonsFromHadronic->Divide(gYMCAntiProtonsFromHadronic);
  SetError(gYESDAntiProtonsFromHadronic,gYMCAntiProtonsFromHadronic);
  gYESDAntiProtonsFromHadronic->Scale(100.);
  gYESDAntiProtonsFromHadronic->SetMarkerStyle(22);
  gYESDAntiProtonsFromHadronic->SetMarkerColor(3);

  c14->cd(1)->SetBottomMargin(0.15); 
  c14->cd(1)->SetLeftMargin(0.15); 
  c14->cd(1)->SetGridx(); c14->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(gYPrimariesESDAntiProtons->GetXaxis()->GetXmin()-0.2,
				   gYPrimariesESDAntiProtons->GetXaxis()->GetXmax()+0.2);
  hEmpty->GetXaxis()->SetTitle(gYPrimariesESDAntiProtons->GetXaxis()->GetTitle());
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gYPrimariesESDProtons->DrawCopy("ESAME");
  gYESDProtonsFromWeak->DrawCopy("ESAME");
  gYESDProtonsFromHadronic->DrawCopy("ESAME");

  c14->cd(2)->SetBottomMargin(0.15); 
  c14->cd(2)->SetLeftMargin(0.15); 
  c14->cd(2)->SetGridx(); c14->cd(2)->SetGridy();
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gYPrimariesESDAntiProtons->DrawCopy("ESAME");
  gYESDAntiProtonsFromWeak->DrawCopy("ESAME");
  gYESDAntiProtonsFromHadronic->DrawCopy("ESAME");
  c14->SaveAs("ReconstructionEfficiency-Protons-Rapidity.gif");

  //pT dependence
  TCanvas *c15 = new TCanvas("c15",
			     "(Anti)Proton reconstruction efficiency vs pT",
			     700,700,700,400);
  c15->SetHighLightColor(10); c15->Divide(2,1);

  //Primary Protons
  TH1D *gPtPrimariesESDProtons = (TH1D *)gHistPrimariesESDYPtProtons->ProjectionY("gPtPrimariesESDProtons",0,gHistPrimariesESDYPtProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtPrimariesMCProtons = (TH1D *)gHistPrimariesMCYPtProtons->ProjectionY("gPtPrimariesMCProtons",0,gHistPrimariesMCYPtProtons->GetYaxis()->GetNbins(),"e");
  gPtPrimariesESDProtons->Divide(gPtPrimariesMCProtons);
  SetError(gPtPrimariesESDProtons,gPtPrimariesMCProtons);
  gPtPrimariesESDProtons->Scale(100.);
  gPtPrimariesESDProtons->SetMarkerStyle(kFullCircle);

  //Primary AntiProtons
  TH1D *gPtPrimariesESDAntiProtons = (TH1D *)gHistPrimariesESDYPtAntiProtons->ProjectionY("gPtPrimariesESDAntiProtons",0,gHistPrimariesESDYPtAntiProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtPrimariesMCAntiProtons = (TH1D *)gHistPrimariesMCYPtAntiProtons->ProjectionY("gPtPrimariesMCAntiProtons",0,gHistPrimariesMCYPtProtons->GetYaxis()->GetNbins(),"e");
  gPtPrimariesESDAntiProtons->Divide(gPtPrimariesMCAntiProtons);
  SetError(gPtPrimariesESDAntiProtons,gPtPrimariesMCAntiProtons);
  gPtPrimariesESDAntiProtons->Scale(100.);
  gPtPrimariesESDAntiProtons->SetMarkerStyle(kFullCircle);

  //Protons from weak decays
  TH1D *gPtESDProtonsFromWeak = (TH1D *)gHistESDYPtProtonsFromWeak->ProjectionY("gPtESDProtonsFromWeak",0,gHistESDYPtProtonsFromWeak->GetYaxis()->GetNbins(),"e");
  TH1D *gPtMCProtonsFromWeak = (TH1D *)gHistMCYPtProtonsFromWeak->ProjectionY("gPtMCProtonsFromWeak",0,gHistMCYPtProtonsFromWeak->GetYaxis()->GetNbins(),"e");
  gPtESDProtonsFromWeak->Divide(gPtMCProtonsFromWeak);
  SetError(gPtESDProtonsFromWeak,gPtMCProtonsFromWeak);
  gPtESDProtonsFromWeak->Scale(100.);
  gPtESDProtonsFromWeak->SetMarkerStyle(21);
  gPtESDProtonsFromWeak->SetMarkerColor(2);

  //AntiProtons from weak decays
  TH1D *gPtESDAntiProtonsFromWeak = (TH1D *)gHistESDYPtAntiProtonsFromWeak->ProjectionY("gPtESDAntiProtonsFromWeak",0,gHistESDYPtAntiProtonsFromWeak->GetYaxis()->GetNbins(),"e");
  TH1D *gPtMCAntiProtonsFromWeak = (TH1D *)gHistMCYPtAntiProtonsFromWeak->ProjectionY("gPtMCAntiProtonsFromWeak",0,gHistMCYPtProtonsFromWeak->GetYaxis()->GetNbins(),"e");
  gPtESDAntiProtonsFromWeak->Divide(gPtMCAntiProtonsFromWeak);
  SetError(gPtESDAntiProtonsFromWeak,gPtMCAntiProtonsFromWeak);
  gPtESDAntiProtonsFromWeak->Scale(100.);
  gPtESDAntiProtonsFromWeak->SetMarkerStyle(21);
  gPtESDAntiProtonsFromWeak->SetMarkerColor(2);

  //Protons from hadronic interactions
  TH1D *gPtESDProtonsFromHadronic = (TH1D *)gHistESDYPtProtonsFromHadronic->ProjectionY("gPtESDProtonsFromHadronic",0,gHistESDYPtProtonsFromHadronic->GetYaxis()->GetNbins(),"e");
  TH1D *gPtMCProtonsFromHadronic = (TH1D *)gHistMCYPtProtonsFromHadronic->ProjectionY("gPtMCProtonsFromHadronic",0,gHistMCYPtProtonsFromHadronic->GetYaxis()->GetNbins(),"e");
  gPtESDProtonsFromHadronic->Divide(gPtMCProtonsFromHadronic);
  SetError(gPtESDProtonsFromHadronic,gPtMCProtonsFromHadronic);
  gPtESDProtonsFromHadronic->Scale(100.);
  gPtESDProtonsFromHadronic->SetMarkerStyle(22);
  gPtESDProtonsFromHadronic->SetMarkerColor(3);

  //AntiProtons from hadronic interactions
  TH1D *gPtESDAntiProtonsFromHadronic = (TH1D *)gHistESDYPtAntiProtonsFromHadronic->ProjectionY("gPtESDAntiProtonsFromHadronic",0,gHistESDYPtAntiProtonsFromHadronic->GetYaxis()->GetNbins(),"e");
  TH1D *gPtMCAntiProtonsFromHadronic = (TH1D *)gHistMCYPtAntiProtonsFromHadronic->ProjectionY("gPtMCAntiProtonsFromHadronic",0,gHistMCYPtProtonsFromHadronic->GetYaxis()->GetNbins(),"e");
  gPtESDAntiProtonsFromHadronic->Divide(gPtMCAntiProtonsFromHadronic);
  SetError(gPtESDAntiProtonsFromHadronic,gPtMCAntiProtonsFromHadronic);
  gPtESDAntiProtonsFromHadronic->Scale(100.);
  gPtESDAntiProtonsFromHadronic->SetMarkerStyle(22);
  gPtESDAntiProtonsFromHadronic->SetMarkerColor(3);


  c15->cd(1)->SetBottomMargin(0.15); 
  c15->cd(1)->SetLeftMargin(0.15); 
  c15->cd(1)->SetGridx(); c15->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(gPtPrimariesESDAntiProtons->GetXaxis()->GetXmin()-0.2,
				   gPtPrimariesESDAntiProtons->GetXaxis()->GetXmax()+0.2);
  hEmpty->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gPtPrimariesESDProtons->DrawCopy("ESAME");
  gPtESDProtonsFromWeak->DrawCopy("ESAME");
  gPtESDProtonsFromHadronic->DrawCopy("ESAME");

  c15->cd(2)->SetBottomMargin(0.15); 
  c15->cd(2)->SetLeftMargin(0.15); 
  c15->cd(2)->SetGridx(); c15->cd(2)->SetGridy();
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gPtPrimariesESDAntiProtons->DrawCopy("ESAME");
  gPtESDAntiProtonsFromWeak->DrawCopy("ESAME");
  gPtESDAntiProtonsFromHadronic->DrawCopy("ESAME");
  c15->SaveAs("ReconstructionEfficiency-Protons-Pt.gif");

  //______________//
  //PID efficiency//
  //______________//
  TH3D *gHistESDInitYPtProtons = (TH3D *)list->At(12);
  TH3D *gHistESDIdYPtProtons = (TH3D *)list->At(13);
  TH3D *gHistESDRecIdYPtProtons = (TH3D *)list->At(14);
  TH3D *gHistESDContamYPtProtons = (TH3D *)list->At(15);

  TCanvas *c16 = new TCanvas("c16",
			     "(Anti)Proton PID efficiency vs y and pT",
			     750,750,700,400);
  c16->SetHighLightColor(10); c16->Divide(3,1);

  //rapidity dependence
  //protons pid efficiency
  //TH1D *gYESDIdProtons = (TH1D *)gHistESDIdYPtProtons->ProjectionX("gYESDIdProtons",0,gHistESDIdYPtProtons->GetXaxis()->GetNbins(),"e");
  //TH1D *gYESDInitProtons = (TH1D *)gHistESDInitYPtProtons->ProjectionX("gYESDInitProtons",0,gHistESDInitYPtProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYESDIdProtons = (TH1D *)gHistESDIdYPtProtons->Project3D("x");
  TH1D *gYESDInitProtons = (TH1D *)gHistESDInitYPtProtons->Project3D("x");
  gYESDIdProtons->Divide(gYESDInitProtons);
  SetError(gYESDIdProtons,gYESDInitProtons);
  gYESDIdProtons->Scale(100.);
  gYESDIdProtons->SetMarkerStyle(kFullCircle);

  //protons pid contamination
  //TH1D *gYESDContamProtons = (TH1D *)gHistESDContamYPtProtons->ProjectionX("gYESDContamProtons",0,gHistESDContamYPtProtons->GetXaxis()->GetNbins(),"e");
  //TH1D *gYESDRecIdProtons = (TH1D *)gHistESDRecIdYPtProtons->ProjectionX("gYESDRecIdProtons",0,gHistESDRecIdYPtProtons->GetXaxis()->GetNbins(),"e");
  TH1D *gYESDContamProtons = (TH1D *)gHistESDContamYPtProtons->Project3D("x");
  TH1D *gYESDRecIdProtons = (TH1D *)gHistESDRecIdYPtProtons->Project3D("x");
  gYESDContamProtons->Divide(gYESDRecIdProtons);
  SetError(gYESDContamProtons,gYESDRecIdProtons);
  gYESDContamProtons->Scale(100.);
  gYESDContamProtons->SetMarkerStyle(kOpenCircle);

  c16->cd(1)->SetBottomMargin(0.15); 
  c16->cd(1)->SetLeftMargin(0.15); 
  c16->cd(1)->SetGridx(); c16->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(-1.0.,1.0);
  hEmpty->GetXaxis()->SetTitle(gYESDContamProtons->GetXaxis()->GetTitle());
  hEmpty->DrawCopy();
  gYESDIdProtons->DrawCopy("ESAME");
  gYESDContamProtons->DrawCopy("ESAME");

  //pT dependence
  //protons pid efficiency
  //TH1D *gPtESDIdProtons = (TH1D *)gHistESDIdYPtProtons->ProjectionY("gPtESDIdProtons",0,gHistESDIdYPtProtons->GetYaxis()->GetNbins(),"e");
  //TH1D *gPtESDInitProtons = (TH1D *)gHistESDInitYPtProtons->ProjectionY("gPtESDInitProtons",0,gHistESDInitYPtProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtESDIdProtons = (TH1D *)gHistESDIdYPtProtons->Project3D("y");
  TH1D *gPtESDInitProtons = (TH1D *)gHistESDInitYPtProtons->Project3D("y");
  gPtESDIdProtons->Divide(gPtESDInitProtons);
  SetError(gPtESDIdProtons,gPtESDInitProtons);
  gPtESDIdProtons->Scale(100.);
  gPtESDIdProtons->SetMarkerStyle(kFullCircle);

  //protons pid contamination
  //TH1D *gPtESDContamProtons = (TH1D *)gHistESDContamYPtProtons->ProjectionY("gPtESDContamProtons",0,gHistESDContamYPtProtons->GetYaxis()->GetNbins(),"e");
  //TH1D *gPtESDRecIdProtons = (TH1D *)gHistESDRecIdYPtProtons->ProjectionY("gPtESDRecIdProtons",0,gHistESDRecIdYPtProtons->GetYaxis()->GetNbins(),"e");
  TH1D *gPtESDContamProtons = (TH1D *)gHistESDContamYPtProtons->Project3D("y");
  TH1D *gPtESDRecIdProtons = (TH1D *)gHistESDRecIdYPtProtons->Project3D("y");
  gPtESDContamProtons->Divide(gPtESDRecIdProtons);
  SetError(gPtESDContamProtons,gPtESDRecIdProtons);
  gPtESDContamProtons->Scale(100.);
  gPtESDContamProtons->SetMarkerStyle(kOpenCircle);

  c16->cd(2)->SetBottomMargin(0.15); 
  c16->cd(2)->SetLeftMargin(0.15); 
  c16->cd(2)->SetGridx(); c16->cd(2)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(0.0,1.2);
  hEmpty->GetXaxis()->SetTitle(gPtESDContamProtons->GetXaxis()->GetTitle());
  hEmpty->DrawCopy();
  gPtESDIdProtons->DrawCopy("ESAME");
  gPtESDContamProtons->DrawCopy("ESAME");

  //N_points dependence
  //protons pid efficiency
  TH1D *gNPointsESDIdProtons = (TH1D *)gHistESDIdYPtProtons->Project3D("y");
  TH1D *gNPointsESDInitProtons = (TH1D *)gHistESDInitYPtProtons->Project3D("y");
  gNPointsESDIdProtons->Divide(gNPointsESDInitProtons);
  SetError(gNPointsESDIdProtons,gNPointsESDInitProtons);
  gNPointsESDIdProtons->Scale(100.);
  gNPointsESDIdProtons->SetMarkerStyle(kFullCircle);

  //protons pid contamination
  TH1D *gNPointsESDContamProtons = (TH1D *)gHistESDContamYPtProtons->Project3D("z");
  TH1D *gNPointsESDRecIdProtons = (TH1D *)gHistESDRecIdYPtProtons->Project3D("z");
  gNPointsESDContamProtons->Divide(gNPointsESDRecIdProtons);
  SetError(gNPointsESDContamProtons,gNPointsESDRecIdProtons);
  gNPointsESDContamProtons->Scale(100.);
  gNPointsESDContamProtons->SetMarkerStyle(kOpenCircle);

  c16->cd(3)->SetBottomMargin(0.15); 
  c16->cd(3)->SetLeftMargin(0.15); 
  c16->cd(3)->SetGridx(); c16->cd(2)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(0.0,1.2);
  hEmpty->GetXaxis()->SetTitle(gNPointsESDContamProtons->GetXaxis()->GetTitle());
  hEmpty->DrawCopy();
  gNPointsESDIdProtons->DrawCopy("ESAME");
  gNPointsESDContamProtons->DrawCopy("ESAME");

  c16->SaveAs("PIDEfficiency-Protons.gif");

  TString outputFileName = "Reconstruction-PID-Efficiency.";
  outputFileName += analysisType; outputFileName += ".root";
  TFile *fout = TFile::Open(outputFileName.Data(),"recreate");
  gYPrimariesESDProtons->Write();
  gYESDProtonsFromWeak->Write();
  gYESDProtonsFromHadronic->Write();
  gPtPrimariesESDProtons->Write();
  gPtESDProtonsFromWeak->Write();
  gPtESDProtonsFromHadronic->Write();
  gYPrimariesESDAntiProtons->Write();
  gYESDAntiProtonsFromWeak->Write();
  gYESDAntiProtonsFromHadronic->Write();
  gPtPrimariesESDAntiProtons->Write();
  gPtESDAntiProtonsFromWeak->Write();
  gPtESDAntiProtonsFromHadronic->Write();
  gYESDIdProtons->Write();
  gYESDContamProtons->Write();
  gPtESDIdProtons->Write();
  gPtESDContamProtons->Write();
  gNPointsESDIdProtons->Write();
  gNPointsESDContamProtons->Write();
  fout->Close();
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

//________________________________________________//
void drawKineQA(const char *filename) {
  //Draws the QA plots for the kinematic variables for protons and antiprotons
  gStyle->SetPalette(1,0);
  gStyle->SetCanvasColor(41);
  gStyle->SetFrameFillColor(10);

  TFile *f = TFile::Open(filename);
  TList *acceptedList = (TList *)f->Get("acceptedCutList");
  TH3D *gHistEtaPhiNClustersPrimaryProtonsPass = (TH3D *)acceptedList->At(44);
  TH3D *gHistEtaPhiNClustersSecondaryProtonsPass = (TH3D *)acceptedList->At(46);
  TH3D *gHistEtaPhiNClustersPrimaryAntiProtonsPass = (TH3D *)acceptedList->At(45);
  TH3D *gHistEtaPhiNClustersSecondaryAntiProtonsPass = (TH3D *)acceptedList->At(47);

  TList *rejectedList = (TList *)f->Get("rejectedCutList");
  TH3D *gHistEtaPhiNClustersPrimaryProtonsReject = (TH3D *)rejectedList->At(0);
  gHistEtaPhiNClustersPrimaryProtonsPass->Add(gHistEtaPhiNClustersPrimaryProtonsReject);
  TH3D *gHistEtaPhiNClustersSecondaryProtonsReject = (TH3D *)rejectedList->At(2);
  gHistEtaPhiNClustersSecondaryProtonsPass->Add(gHistEtaPhiNClustersSecondaryProtonsReject);
  TH3D *gHistEtaPhiNClustersPrimaryAntiProtonsReject = (TH3D *)rejectedList->At(1);
  gHistEtaPhiNClustersPrimaryAntiProtonsPass->Add(gHistEtaPhiNClustersPrimaryAntiProtonsReject);
  TH3D *gHistEtaPhiNClustersSecondaryAntiProtonsReject = (TH3D *)rejectedList->At(3);
  gHistEtaPhiNClustersSecondaryAntiProtonsPass->Add(gHistEtaPhiNClustersSecondaryAntiProtonsReject);

  //eta-phi
  TCanvas *c21 = new TCanvas("c21",
			     "#eta-#phi",
			     0,0,600,600);
  c21->SetHighLightColor(10); c21->Divide(2,2);
  c21->cd(1)->SetLeftMargin(0.15); c21->cd(1)->SetBottomMargin(0.15);  
  c21->cd(1)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryProtonsPass->Project3D("yx")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryProtonsPass->Project3D("yx")))->DrawCopy("colz");
  c21->cd(2)->SetLeftMargin(0.15); c21->cd(2)->SetBottomMargin(0.15);  
  c21->cd(2)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryProtonsPass->Project3D("yx")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryProtonsPass->Project3D("yx")))->DrawCopy("colz");
  c21->cd(3)->SetLeftMargin(0.15); c21->cd(3)->SetBottomMargin(0.15);  
  c21->cd(3)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryAntiProtonsPass->Project3D("yx")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryAntiProtonsPass->Project3D("yx")))->DrawCopy("colz");
  c21->cd(4)->SetLeftMargin(0.15); c21->cd(4)->SetBottomMargin(0.15);  
  c21->cd(4)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryAntiProtonsPass->Project3D("yx")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryAntiProtonsPass->Project3D("yx")))->DrawCopy("colz");
  c21->SaveAs("EtaPhi.gif");

  //eta-Nclusters
  TCanvas *c22 = new TCanvas("c22",
			     "#eta-N_{clusters}",
			     100,100,600,600);
  c22->SetHighLightColor(10); c22->Divide(2,2);
  c22->cd(1)->SetLeftMargin(0.15); c22->cd(1)->SetBottomMargin(0.15);  
  c22->cd(1)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryProtonsPass->Project3D("zx")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryProtonsPass->Project3D("zx")))->DrawCopy("colz");
  c22->cd(2)->SetLeftMargin(0.15); c22->cd(2)->SetBottomMargin(0.15);  
  c22->cd(2)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryProtonsPass->Project3D("zx")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryProtonsPass->Project3D("zx")))->DrawCopy("colz");
  c22->cd(3)->SetLeftMargin(0.15); c22->cd(3)->SetBottomMargin(0.15);  
  c22->cd(3)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryAntiProtonsPass->Project3D("zx")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryAntiProtonsPass->Project3D("zx")))->DrawCopy("colz");
  c22->cd(4)->SetLeftMargin(0.15); c22->cd(4)->SetBottomMargin(0.15);  
  c22->cd(4)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryAntiProtonsPass->Project3D("zx")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryAntiProtonsPass->Project3D("zx")))->DrawCopy("colz");
  c22->SaveAs("EtaNClusters.gif");

  //phi-Nclusters
  TCanvas *c23 = new TCanvas("c23",
			     "#phi-N_{clusters}",
			     200,200,600,600);
  c23->SetHighLightColor(10); c23->Divide(2,2);
  c23->cd(1)->SetLeftMargin(0.15); c23->cd(1)->SetBottomMargin(0.15);  
  c23->cd(1)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryProtonsPass->Project3D("zy")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryProtonsPass->Project3D("zy")))->DrawCopy("colz");
  c23->cd(2)->SetLeftMargin(0.15); c23->cd(2)->SetBottomMargin(0.15);  
  c23->cd(2)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryProtonsPass->Project3D("zy")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryProtonsPass->Project3D("zy")))->DrawCopy("colz");
  c23->cd(3)->SetLeftMargin(0.15); c23->cd(3)->SetBottomMargin(0.15);  
  c23->cd(3)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryAntiProtonsPass->Project3D("zy")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersPrimaryAntiProtonsPass->Project3D("zy")))->DrawCopy("colz");
  c23->cd(4)->SetLeftMargin(0.15); c23->cd(4)->SetBottomMargin(0.15);  
  c23->cd(4)->SetRightMargin(0.2);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryAntiProtonsPass->Project3D("zy")))->SetStats(kFALSE);
  ((TH2D *)(gHistEtaPhiNClustersSecondaryAntiProtonsPass->Project3D("zy")))->DrawCopy("colz");
  c23->SaveAs("PhiNClusters.gif");

  f->Close();
}

//________________________________________________//
void drawEfficiencies(const char *filename,
		      Bool_t gShowPrimaries = kTRUE,
		      Bool_t gShowWeak = kFALSE,
		      Bool_t gShowHadronic = kFALSE) {
  //Macro to display the reconstruction, cut and pid efficiencies
  gStyle->SetPalette(1,0);
  gStyle->SetCanvasColor(41);
  gStyle->SetFrameFillColor(10);

  TH2F *hEmpty = new TH2F("hEmpty","",100,-1.5,5.0,100,0,120);
  hEmpty->SetStats(kFALSE);
  hEmpty->GetXaxis()->SetTitleColor(1);
  hEmpty->GetXaxis()->SetNdivisions(15);
  hEmpty->GetYaxis()->SetNdivisions(15);
  hEmpty->GetYaxis()->SetTitleOffset(1.4);
  hEmpty->GetYaxis()->SetTitle("#epsilon");

  TLatex *t1 = new TLatex();
  t1->SetTextSize(0.04);

  TPaveText *tpave = new TPaveText();
  tpave->SetFillColor(10);
  Double_t bottomY = 0.0;

  TFile *f = TFile::Open(filename);
  TH1D *gYPrimariesESDProtons = (TH1D *)f->Get("gYPrimariesESDProtons");
  TH1D *gYESDProtonsFromWeak = (TH1D *)f->Get("gYESDProtonsFromWeak");
  TH1D *gYESDProtonsFromHadronic = (TH1D *)f->Get("gYESDProtonsFromHadronic");
  TH1D *gPtPrimariesESDProtons = (TH1D *)f->Get("gPtPrimariesESDProtons");
  TH1D *gPtESDProtonsFromWeak = (TH1D *)f->Get("gPtESDProtonsFromWeak");
  TH1D *gPtESDProtonsFromHadronic = (TH1D *)f->Get("gPtESDProtonsFromHadronic");
  TH1D *gYPrimariesESDAntiProtons = (TH1D *)f->Get("gYPrimariesESDAntiProtons");
  TH1D *gYESDAntiProtonsFromWeak = (TH1D *)f->Get("gYESDAntiProtonsFromWeak");
  TH1D *gYESDAntiProtonsFromHadronic = (TH1D *)f->Get("gYESDAntiProtonsFromHadronic");	
  TH1D *gPtPrimariesESDAntiProtons = (TH1D *)f->Get("gPtPrimariesESDAntiProtons");	
  TH1D *gPtESDAntiProtonsFromWeak = (TH1D *)f->Get("gPtESDAntiProtonsFromWeak");	
  TH1D *gPtESDAntiProtonsFromHadronic = (TH1D *)f->Get("gPtESDAntiProtonsFromHadronic");	
  TH1D *gYESDIdProtons = (TH1D *)f->Get("gYESDIdProtons");
  TH1D *gYESDContamProtons = (TH1D *)f->Get("gYESDContamProtons");
  TH1D *gPtESDIdProtons = (TH1D *)f->Get("gPtESDIdProtons");
  TH1D *gPtESDContamProtons = (TH1D *)f->Get("gPtESDContamProtons");	

  //Reconstruction efficiencies - protons
  TCanvas *c1 = new TCanvas("c1","Reconstruction efficiencies - Protons",
			    0,0,700,400);
  c1->SetFillColor(41); c1->SetHighLightColor(41);
  c1->Divide(2,1);
  c1->cd(1)->SetLeftMargin(0.15); c1->cd(1)->SetBottomMargin(0.15);
  c1->cd(1)->SetGridx(); c1->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(gYPrimariesESDProtons->GetXaxis()->GetXmin()-0.1,
				   gYPrimariesESDProtons->GetXaxis()->GetXmax()+0.1);
  hEmpty->GetXaxis()->SetTitle(gYPrimariesESDProtons->GetXaxis()->GetTitle());
  hEmpty->DrawCopy();
  if(gShowPrimaries)
    bottomY = 111.0;
  if(gShowWeak)
    bottomY = 106.0;
  if(gShowHadronic)
    bottomY = 101.0;
  tpave->DrawPave(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.1,119,
		  gYPrimariesESDProtons->GetXaxis()->GetXmax()-0.1,bottomY);
  if(gShowPrimaries) {
    DrawMarker(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.2, 115, 20, 1.2, 1);
    t1->DrawLatex(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.3,113,"Primary p");
    gYPrimariesESDProtons->DrawCopy("ESAME");
  }
  if(gShowWeak) {
    DrawMarker(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.2, 110, 21, 1.2, 2);
    t1->DrawLatex(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.3,108,"#Lambda -> p + #pi^{-}");
    gYESDProtonsFromWeak->DrawCopy("ESAME");
  }
  if(gShowHadronic) {
    DrawMarker(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.2, 105, 22, 1.2, 3);
    t1->DrawLatex(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.3,103,"X + A -> p + B");
    gYESDProtonsFromHadronic->DrawCopy("ESAME");
  }

  c1->cd(2)->SetLeftMargin(0.15); c1->cd(2)->SetBottomMargin(0.15);
  c1->cd(2)->SetGridx(); c1->cd(2)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(gPtPrimariesESDProtons->GetXaxis()->GetXmin()-0.1,
				   gPtPrimariesESDProtons->GetXaxis()->GetXmax()+0.1);
  hEmpty->GetXaxis()->SetTitle(gPtPrimariesESDProtons->GetXaxis()->GetTitle());
  hEmpty->DrawCopy();
  if(gShowPrimaries)
    bottomY = 111.0;
  if(gShowWeak)
    bottomY = 106.0;
  if(gShowHadronic)
    bottomY = 101.0;
  tpave->DrawPave(gPtPrimariesESDProtons->GetXaxis()->GetXmin(),119,
		  gPtPrimariesESDProtons->GetXaxis()->GetXmax(),bottomY);
  if(gShowPrimaries) {
    DrawMarker(gPtPrimariesESDProtons->GetXaxis()->GetXmin()+0.05, 115, 20, 1.2, 1);
    t1->DrawLatex(gPtPrimariesESDProtons->GetXaxis()->GetXmin()+0.15,113,"Primary p");
    gPtPrimariesESDProtons->DrawCopy("ESAME");
  }
  if(gShowWeak) {
    DrawMarker(gPtPrimariesESDProtons->GetXaxis()->GetXmin()+0.05, 110, 21, 1.2, 2);
    t1->DrawLatex(gPtPrimariesESDProtons->GetXaxis()->GetXmin()+0.15,108,"#Lambda -> p + #pi^{-}");
    gPtESDProtonsFromWeak->DrawCopy("ESAME");
  }
  if(gShowHadronic) {
    DrawMarker(gPtPrimariesESDProtons->GetXaxis()->GetXmin()+0.05, 105, 22, 1.2, 3);
    t1->DrawLatex(gPtPrimariesESDProtons->GetXaxis()->GetXmin()+0.15,103,"X + A -> p + B");
    gPtESDProtonsFromHadronic->DrawCopy("ESAME");
  }
  c1->SaveAs("ReconstructionEfficiency-Protons.gif");

  //Reconstruction efficiencies - antiprotons
  TCanvas *c2 = new TCanvas("c2","Reconstruction efficiencies - Antirotons",
			    100,100,700,400);
  c2->SetFillColor(41); c2->SetHighLightColor(41);
  c2->Divide(2,1);
  c2->cd(1)->SetLeftMargin(0.15); c2->cd(1)->SetBottomMargin(0.15);
  c2->cd(1)->SetGridx(); c2->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(gYPrimariesESDAntiProtons->GetXaxis()->GetXmin()-0.1,
				   gYPrimariesESDAntiProtons->GetXaxis()->GetXmax()+0.1);
  hEmpty->GetXaxis()->SetTitle(gYPrimariesESDAntiProtons->GetXaxis()->GetTitle());
  hEmpty->DrawCopy();
  if(gShowPrimaries)
    bottomY = 111.0;
  if(gShowWeak)
    bottomY = 104.0;
  if(gShowHadronic)
    bottomY = 97.0;
  tpave->DrawPave(gYPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.1,119,
		  gYPrimariesESDAntiProtons->GetXaxis()->GetXmax()-0.1,bottomY);
  if(gShowPrimaries) {
    DrawMarker(gYPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.2, 115, 20, 1.2, 1);
    t1->DrawLatex(gYPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.3,113,"Primary #bar{p}");
    gYPrimariesESDAntiProtons->DrawCopy("ESAME");
  }
  if(gShowWeak) {
    DrawMarker(gYPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.2, 108, 21, 1.2, 2);
    t1->DrawLatex(gYPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.3,106,"#bar{#Lambda} -> #bar{p} + #pi^{+}");
    gYESDAntiProtonsFromWeak->DrawCopy("ESAME");
  }
  if(gShowHadronic) {
    DrawMarker(gYPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.2, 101, 22, 1.2, 3);
    t1->DrawLatex(gYPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.3,98,"X + A -> #bar{p} + B");
    gYESDAntiProtonsFromHadronic->DrawCopy("ESAME");
  }

  c2->cd(2)->SetLeftMargin(0.15); c2->cd(2)->SetBottomMargin(0.15);
  c2->cd(2)->SetGridx(); c2->cd(2)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(gPtPrimariesESDAntiProtons->GetXaxis()->GetXmin()-0.1,
				   gPtPrimariesESDAntiProtons->GetXaxis()->GetXmax()+0.1);
  hEmpty->GetXaxis()->SetTitle(gPtPrimariesESDAntiProtons->GetXaxis()->GetTitle());
  hEmpty->DrawCopy();
  hEmpty->DrawCopy();
  if(gShowPrimaries)
    bottomY = 111.0;
  if(gShowWeak)
    bottomY = 104.0;
  if(gShowHadronic)
    bottomY = 97.0;
  tpave->DrawPave(gPtPrimariesESDAntiProtons->GetXaxis()->GetXmin(),119,
		  gPtPrimariesESDAntiProtons->GetXaxis()->GetXmax(),bottomY);
  if(gShowPrimaries) {
    DrawMarker(gPtPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.05, 115, 20, 1.2, 1);
    t1->DrawLatex(gPtPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.15,113,"Primary #bar{p}");
    gPtPrimariesESDAntiProtons->DrawCopy("ESAME");
  }
  if(gShowWeak) {
    DrawMarker(gPtPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.05, 108, 21, 1.2, 2);
    t1->DrawLatex(gPtPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.15,106,"#bar{#Lambda} -> #bar{p} + #pi^{+}");
    gPtESDAntiProtonsFromWeak->DrawCopy("ESAME");
  }
  if(gShowHadronic) {
    DrawMarker(gPtPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.05, 101, 22, 1.2, 3);
    t1->DrawLatex(gPtPrimariesESDAntiProtons->GetXaxis()->GetXmin()+0.15,98,"X + A -> #bar{p} + B");
    gPtESDAntiProtonsFromHadronic->DrawCopy("ESAME");
  }
  c2->SaveAs("ReconstructionEfficiency-AntiProtons.gif");

  //PID efficiencies - (anti)protons
  TCanvas *c3 = new TCanvas("c3","PID efficiencies",
			    200,200,700,400);
  c3->SetFillColor(41); c3->SetHighLightColor(41);
  c3->Divide(2,1);
  c3->cd(1)->SetLeftMargin(0.15); c3->cd(1)->SetBottomMargin(0.15);
  c3->cd(1)->SetGridx(); c3->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(gYPrimariesESDProtons->GetXaxis()->GetXmin()-0.1,
				   gYPrimariesESDProtons->GetXaxis()->GetXmax()+0.1);
  hEmpty->GetXaxis()->SetTitle(gYPrimariesESDProtons->GetXaxis()->GetTitle());
  hEmpty->DrawCopy();
  tpave->DrawPave(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.1,119,
		  gYPrimariesESDProtons->GetXaxis()->GetXmax()-0.1,104);
  gYESDIdProtons->DrawCopy("ESAME");
  gYESDContamProtons->DrawCopy("ESAME");
  DrawMarker(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.2, 115, 20, 1.2, 1);
  t1->DrawLatex(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.3,113,"Efficiency");
  DrawMarker(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.2, 108, 24, 1.2, 1);
  t1->DrawLatex(gYPrimariesESDProtons->GetXaxis()->GetXmin()+0.3,106,"Contamination");
 
  c3->cd(2)->SetLeftMargin(0.15); c3->cd(2)->SetBottomMargin(0.15);
  c3->cd(2)->SetGridx(); c3->cd(2)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(gPtPrimariesESDProtons->GetXaxis()->GetXmin()-0.1,
				   gPtPrimariesESDProtons->GetXaxis()->GetXmax()+0.1);
  hEmpty->GetXaxis()->SetTitle(gPtPrimariesESDProtons->GetXaxis()->GetTitle());
  hEmpty->DrawCopy();
  tpave->DrawPave(gPtPrimariesESDProtons->GetXaxis()->GetXmin(),119,
		  gPtPrimariesESDProtons->GetXaxis()->GetXmax(),104);
  gPtESDIdProtons->DrawCopy("ESAME");
  gPtESDContamProtons->DrawCopy("ESAME");
  DrawMarker(gPtPrimariesESDProtons->GetXaxis()->GetXmin()+0.05, 115, 20, 1.2, 1);
  t1->DrawLatex(gPtPrimariesESDProtons->GetXaxis()->GetXmin()+0.15,113,"Efficiency");
  DrawMarker(gPtPrimariesESDProtons->GetXaxis()->GetXmin()+0.05, 108, 24, 1.2, 1);
  t1->DrawLatex(gPtPrimariesESDProtons->GetXaxis()->GetXmin()+0.15,106,"Contamination");
  c3->SaveAs("PIDEfficiency-Protons.gif");
}

//________________________________________________//
void compareEfficiencies(const char *filenameTPC,
			 const char *filenameHybrid,
			 Bool_t gShowPrimaries = kTRUE,
			 Bool_t gShowWeak = kFALSE,
			 Bool_t gShowHadronic = kFALSE) {
  //Function to compare the reconstruction efficiencies between two methods:
  //e.g. TPC standalone and global tracking
  gStyle->SetPalette(1,0);
  gStyle->SetCanvasColor(41);
  gStyle->SetFrameFillColor(10);

  TH2F *hEmpty = new TH2F("hEmpty","",100,-1.5,2.0,100,0,120);
  hEmpty->SetStats(kFALSE);
  hEmpty->GetXaxis()->SetTitleColor(1);
  hEmpty->GetXaxis()->SetNdivisions(15);
  hEmpty->GetYaxis()->SetNdivisions(15);
  hEmpty->GetYaxis()->SetTitleOffset(1.4);
  hEmpty->GetYaxis()->SetTitle("#epsilon");

  TLatex *t1 = new TLatex();
  t1->SetTextSize(0.04);

  TPaveText *tpave = new TPaveText();
  tpave->SetFillColor(10);
  Double_t bottomY = 0.0;

  //TPC standalone
  TFile *f1 = TFile::Open(filenameTPC);
  TH1D *g1YPrimariesESDProtons = (TH1D *)f1->Get("gYPrimariesESDProtons");
  TH1D *g1YESDProtonsFromWeak = (TH1D *)f1->Get("gYESDProtonsFromWeak");
  TH1D *g1YESDProtonsFromHadronic = (TH1D *)f1->Get("gYESDProtonsFromHadronic");
  TH1D *g1PtPrimariesESDProtons = (TH1D *)f1->Get("gPtPrimariesESDProtons");
  TH1D *g1PtESDProtonsFromWeak = (TH1D *)f1->Get("gPtESDProtonsFromWeak");
  TH1D *g1PtESDProtonsFromHadronic = (TH1D *)f1->Get("gPtESDProtonsFromHadronic");
  TH1D *g1YPrimariesESDAntiProtons = (TH1D *)f1->Get("gYPrimariesESDAntiProtons");
  TH1D *g1YESDAntiProtonsFromWeak = (TH1D *)f1->Get("gYESDAntiProtonsFromWeak");
  TH1D *g1YESDAntiProtonsFromHadronic = (TH1D *)f1->Get("gYESDAntiProtonsFromHadronic");	
  TH1D *g1PtPrimariesESDAntiProtons = (TH1D *)f1->Get("gPtPrimariesESDAntiProtons");	
  TH1D *g1PtESDAntiProtonsFromWeak = (TH1D *)f1->Get("gPtESDAntiProtonsFromWeak");	
  TH1D *g1PtESDAntiProtonsFromHadronic = (TH1D *)f1->Get("gPtESDAntiProtonsFromHadronic");	

  //Global tracking
  TFile *f2 = TFile::Open(filenameHybrid);
  TH1D *g2YPrimariesESDProtons = (TH1D *)f2->Get("gYPrimariesESDProtons");
  TH1D *g2YESDProtonsFromWeak = (TH1D *)f2->Get("gYESDProtonsFromWeak");
  TH1D *g2YESDProtonsFromHadronic = (TH1D *)f2->Get("gYESDProtonsFromHadronic");
  TH1D *g2PtPrimariesESDProtons = (TH1D *)f2->Get("gPtPrimariesESDProtons");
  TH1D *g2PtESDProtonsFromWeak = (TH1D *)f2->Get("gPtESDProtonsFromWeak");
  TH1D *g2PtESDProtonsFromHadronic = (TH1D *)f2->Get("gPtESDProtonsFromHadronic");
  TH1D *g2YPrimariesESDAntiProtons = (TH1D *)f2->Get("gYPrimariesESDAntiProtons");
  TH1D *g2YESDAntiProtonsFromWeak = (TH1D *)f2->Get("gYESDAntiProtonsFromWeak");
  TH1D *g2YESDAntiProtonsFromHadronic = (TH1D *)f2->Get("gYESDAntiProtonsFromHadronic");	
  TH1D *g2PtPrimariesESDAntiProtons = (TH1D *)f2->Get("gPtPrimariesESDAntiProtons");	
  TH1D *g2PtESDAntiProtonsFromWeak = (TH1D *)f2->Get("gPtESDAntiProtonsFromWeak");	
  TH1D *g2PtESDAntiProtonsFromHadronic = (TH1D *)f2->Get("gPtESDAntiProtonsFromHadronic");	

  //Reconstruction efficiencies - protons
  TCanvas *c1 = new TCanvas("c1","Reconstruction efficiencies - Protons",
			    0,0,700,400);
  c1->SetFillColor(41); c1->SetHighLightColor(41);
  c1->Divide(2,1);
  c1->cd(1)->SetLeftMargin(0.15); c1->cd(1)->SetBottomMargin(0.15);
  c1->cd(1)->SetGridx(); c1->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(-1.1,1.1);
  hEmpty->GetXaxis()->SetTitle("#eta");
  hEmpty->DrawCopy();
  if(gShowPrimaries)
    bottomY = 106.0;
  if(gShowWeak)
    bottomY = 99.0;
  if(gShowHadronic)
    bottomY = 92.0;
  tpave->DrawPave(-0.4,119,1.1,bottomY);
  t1->DrawLatex(-0.3,113,"TPC");
  t1->DrawLatex(0.0,113,"Global");
  if(gShowPrimaries) {
    DrawMarker(-0.15, 110, 20, 1.2, 1);
    DrawMarker(0.15, 110, 24, 1.2, 1);
    t1->DrawLatex(0.35,108,"Primary p");
    g1YPrimariesESDProtons->DrawCopy("ESAME");
    g2YPrimariesESDProtons->SetMarkerStyle(24);
    g2YPrimariesESDProtons->DrawCopy("ESAME");
  }
  if(gShowWeak) {
    DrawMarker(-0.15, 103, 21, 1.2, 2);
    DrawMarker(0.15, 103, 25, 1.2, 2);
    t1->DrawLatex(0.35,101,"#Lambda -> p + #pi^{-}");
    g1YESDProtonsFromWeak->DrawCopy("ESAME");
    g2YESDProtonsFromWeak->SetMarkerStyle(25);
    g2YESDProtonsFromWeak->DrawCopy("ESAME");
  }
  if(gShowHadronic) {
    DrawMarker(-0.15, 96, 22, 1.2, 3);
    DrawMarker(0.15, 96, 26, 1.2, 3);
    t1->DrawLatex(0.35,94,"X + A -> p + B");
    g1YESDProtonsFromHadronic->DrawCopy("ESAME");
    g2YESDProtonsFromHadronic->SetMarkerStyle(26);
    g2YESDProtonsFromHadronic->DrawCopy("ESAME");
  }

  c1->cd(2)->SetLeftMargin(0.15); c1->cd(2)->SetBottomMargin(0.15);
  c1->cd(2)->SetGridx(); c1->cd(2)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(0.0,1.7);
  hEmpty->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  hEmpty->DrawCopy();
  if(gShowPrimaries)
    bottomY = 106.0;
  if(gShowWeak)
    bottomY = 99.0;
  if(gShowHadronic)
    bottomY = 92.0;
  tpave->DrawPave(0.05,119,1.2,bottomY);
  t1->DrawLatex(0.1,113,"TPC");
  t1->DrawLatex(0.3,113,"Global");
  if(gShowPrimaries) {
    DrawMarker(0.15, 110, 20, 1.2, 1);
    DrawMarker(0.4, 110, 24, 1.2, 1);
    t1->DrawLatex(0.6,108,"Primary p");
    g1PtPrimariesESDProtons->DrawCopy("ESAME");
    g2PtPrimariesESDProtons->SetMarkerStyle(24);
    g2PtPrimariesESDProtons->DrawCopy("ESAME");
  }
  if(gShowWeak) {
    DrawMarker(0.15, 103, 21, 1.2, 2);
    DrawMarker(0.4, 103, 25, 1.2, 2);
    t1->DrawLatex(0.6,101,"#Lambda -> p + #pi^{-}");
    g1PtESDProtonsFromWeak->DrawCopy("ESAME");
    g2PtESDProtonsFromWeak->SetMarkerStyle(25);
    g2PtESDProtonsFromWeak->DrawCopy("ESAME");
  }
  if(gShowHadronic) {
    DrawMarker(0.15, 96, 22, 1.2, 3);
    DrawMarker(0.4, 96, 26, 1.2, 3);
    t1->DrawLatex(0.6,94,"X + A -> p + B");
    g1PtESDProtonsFromHadronic->DrawCopy("ESAME");
    g2PtESDProtonsFromHadronic->SetMarkerStyle(26);
    g2PtESDProtonsFromHadronic->DrawCopy("ESAME");
  }
  c1->SaveAs("ReconstructionEfficiency-Protons.gif");

  //Reconstruction efficiencies - antiprotons
  TCanvas *c2 = new TCanvas("c2","Reconstruction efficiencies - Antirotons",
			    100,100,700,400);
  c2->SetFillColor(41); c2->SetHighLightColor(41);
  c2->Divide(2,1);
  c2->cd(1)->SetLeftMargin(0.15); c2->cd(1)->SetBottomMargin(0.15);
  c2->cd(1)->SetGridx(); c2->cd(1)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(-1.1,1.1);
  hEmpty->GetXaxis()->SetTitle("#eta");
  hEmpty->DrawCopy();
  if(gShowPrimaries)
    bottomY = 106.0;
  if(gShowWeak)
    bottomY = 99.0;
  if(gShowHadronic)
    bottomY = 92.0;
  tpave->DrawPave(-0.4,119,1.1,bottomY);
  t1->DrawLatex(-0.3,113,"TPC");
  t1->DrawLatex(0.0,113,"Global");
  if(gShowPrimaries) {
    DrawMarker(-0.15, 110, 20, 1.2, 1);
    DrawMarker(0.15, 110, 24, 1.2, 1);
    t1->DrawLatex(0.35,108,"Primary #bar{p}");
    g1YPrimariesESDAntiProtons->DrawCopy("ESAME");
    g2YPrimariesESDAntiProtons->SetMarkerStyle(24);
    g2YPrimariesESDAntiProtons->DrawCopy("ESAME");
  }
  if(gShowWeak) {
    DrawMarker(-0.15, 103, 21, 1.2, 2);
    DrawMarker(0.15, 103, 25, 1.2, 2);
    t1->DrawLatex(0.35,101,"#bar{#Lambda} -> #bar{p} + #pi^{+}");
    g1YESDAntiProtonsFromWeak->DrawCopy("ESAME");
    g2YESDAntiProtonsFromWeak->SetMarkerStyle(25);
    g2YESDAntiProtonsFromWeak->DrawCopy("ESAME");
  }
  if(gShowHadronic) {
    DrawMarker(-0.15, 96, 22, 1.2, 3);
    DrawMarker(0.15, 96, 26, 1.2, 3);
    t1->DrawLatex(0.35,94,"X + A -> #bar{p} + B");
    g1YESDAntiProtonsFromHadronic->DrawCopy("ESAME");
    g2YESDAntiProtonsFromHadronic->SetMarkerStyle(26);
    g2YESDAntiProtonsFromHadronic->DrawCopy("ESAME");
  }

  c2->cd(2)->SetLeftMargin(0.15); c2->cd(2)->SetBottomMargin(0.15);
  c2->cd(2)->SetGridx(); c2->cd(2)->SetGridy();
  hEmpty->GetXaxis()->SetRangeUser(0.0,1.7);
  hEmpty->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  hEmpty->DrawCopy();
  if(gShowPrimaries)
    bottomY = 106.0;
  if(gShowWeak)
    bottomY = 99.0;
  if(gShowHadronic)
    bottomY = 92.0;
  tpave->DrawPave(0.05,119,1.2,bottomY);
  t1->DrawLatex(0.1,113,"TPC");
  t1->DrawLatex(0.3,113,"Global");
  if(gShowPrimaries) {
    DrawMarker(0.15, 110, 20, 1.2, 1);
    DrawMarker(0.4, 110, 24, 1.2, 1);
    t1->DrawLatex(0.6,108,"Primary #bar{p}");
    g1PtPrimariesESDAntiProtons->DrawCopy("ESAME");
    g2PtPrimariesESDAntiProtons->SetMarkerStyle(24);
    g2PtPrimariesESDAntiProtons->DrawCopy("ESAME");
  }
  if(gShowWeak) {
    DrawMarker(0.15, 103, 21, 1.2, 2);
    DrawMarker(0.4, 103, 25, 1.2, 2);
    t1->DrawLatex(0.6,101,"#bar{#Lambda} -> #bar{p} + #pi^{+}");
    g1PtESDAntiProtonsFromWeak->DrawCopy("ESAME");
    g2PtESDAntiProtonsFromWeak->SetMarkerStyle(25);
    g2PtESDAntiProtonsFromWeak->DrawCopy("ESAME");
  }
  if(gShowHadronic) {
    DrawMarker(0.15, 96, 22, 1.2, 3);
    DrawMarker(0.4, 96, 26, 1.2, 3);
    t1->DrawLatex(0.6,94,"X + A -> #bar{p} + B");
    g1PtESDAntiProtonsFromHadronic->DrawCopy("ESAME");
    g2PtESDAntiProtonsFromHadronic->SetMarkerStyle(26);
    g2PtESDAntiProtonsFromHadronic->DrawCopy("ESAME");
  }
  c2->SaveAs("ReconstructionEfficiency-AntiProtons.gif");

}

//________________________________________________//
void drawCutParametersDistributions(const char* filename = "Protons.QA.Histograms.root") {
  //macro that takes as an input the third file 
  //created by the proton QA analysis task 
  //and draws the DCA distributions of protons 
  //and antiprotons (both primary & secondaries)
  const Int_t nEvents = 1;

  TFile *f = TFile::Open(filename);

  //cut list
  TH1F *gCutListHistograms[100];
  TList *listCut = (TList *)f->Get("acceptedCutList");
  Int_t iCounter = 0;
  cout<<"Cut list entries: "<<listCut->GetEntries()<<endl;
  for(Int_t iEntry = 0; iEntry < listCut->GetEntries(); iEntry++) {
    if(iCounter == 4) iCounter = 0;
    iCounter += 1;
    gCutListHistograms[iEntry] = (TH1F *)listCut->At(iEntry);
    gCutListHistograms[iEntry]->Scale(1./nEvents);
    if(iCounter < 3) {
      gCutListHistograms[iEntry]->SetFillColor(4);
      gCutListHistograms[iEntry]->SetMarkerColor(4);
      gCutListHistograms[iEntry]->SetMarkerStyle(20);
    }
    else {
      gCutListHistograms[iEntry]->SetFillColor(kOrange+1);
      gCutListHistograms[iEntry]->SetMarkerColor(kOrange+1);
      gCutListHistograms[iEntry]->SetMarkerStyle(29);
    }
    /*cout<<"Entry: "<<iEntry<<
      " - Counter: "<<iCounter<<
      " - Name: "<<gCutListHistograms[iEntry]->GetName()<<endl;*/
  }

  //DCA list
  TH1F *gDCAListHistograms[20];
  TList *listDCA = (TList *)f->Get("acceptedDCAList");
  iCounter = 0;
  cout<<"DCA list entries: "<<listDCA->GetEntries()<<endl;
  for(Int_t iEntry = 0; iEntry < listDCA->GetEntries(); iEntry++) {
    if(iCounter == 4) iCounter = 0;
    iCounter += 1;
    gDCAListHistograms[iEntry] = (TH1F *)listDCA->At(iEntry);
    gDCAListHistograms[iEntry]->Scale(1./nEvents);
    if(iCounter < 3) {
      gDCAListHistograms[iEntry]->SetFillColor(4);
      gDCAListHistograms[iEntry]->SetMarkerColor(4);
      gDCAListHistograms[iEntry]->SetMarkerStyle(20);
    }
    else {
      gDCAListHistograms[iEntry]->SetFillColor(kOrange+1);
      gDCAListHistograms[iEntry]->SetMarkerColor(kOrange+1);
      gDCAListHistograms[iEntry]->SetMarkerStyle(29);
    }
    /*cout<<"Entry: "<<iEntry<<
      " - Counter: "<<iCounter<<
      " - Name: "<<gDCAListHistograms[iEntry]->GetName()<<endl;*/
  }

  //_________________________________________________________//
  TF1 *gDCA = new TF1("gDCA",
		      "[0]*TMath::Power(1+TMath::Exp((x-[1])/[2]),-1)",
		      0.1,100.0);
  gDCA->SetParameter(0,1.74221e+07);
  gDCA->SetParameter(1,-1.12221e+01);
  gDCA->SetParameter(2,1.02726);
  //_________________________________________________________//
  TH2F *hEmpty = new TH2F("hEmpty","",300,-100,200,100,1e-01,1e+06);
  hEmpty->GetYaxis()->SetTitle("Entries/Event");
  hEmpty->GetYaxis()->SetNdivisions(10);
  hEmpty->GetXaxis()->SetNdivisions(10);
  hEmpty->SetStats(kFALSE);
  //_________________________________________________________//

  //Cut parameters
  TCanvas *c1 = new TCanvas("c1","ITS Cluster map",0,0,650,350);
  c1->SetFillColor(10); c1->GetFrame()->SetFillColor(10);
  c1->SetHighLightColor(10); c1->Divide(2,1);
  c1->cd(1)->SetBottomMargin(0.2); c1->cd(1)->SetLeftMargin(0.2);
  c1->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("ITS layer");
  hEmpty->GetXaxis()->SetRangeUser(0.0,7.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[0]->Draw("ESAME");
  gCutListHistograms[2]->Draw("ESAME");
  c1->cd(2)->SetBottomMargin(0.2); c1->cd(2)->SetLeftMargin(0.2);
  c1->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("ITS layer");
  hEmpty->GetXaxis()->SetRangeUser(0.0,7.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[1]->Draw("ESAME");
  gCutListHistograms[3]->Draw("ESAME");
  c1->SaveAs("ITSClusterMap.gif");

  TCanvas *c2 = new TCanvas("c2","Number of ITS Clusters",50,50,650,350);
  c2->SetFillColor(10); c2->GetFrame()->SetFillColor(10);
  c2->SetHighLightColor(10); c2->Divide(2,1);
  c2->cd(1)->SetBottomMargin(0.2); c2->cd(1)->SetLeftMargin(0.2);
  c2->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("N_{clusters}(ITS)");
  hEmpty->GetXaxis()->SetRangeUser(0.0,7.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[4]->Draw("ESAME");
  gCutListHistograms[6]->Draw("ESAME");
  c2->cd(2)->SetBottomMargin(0.2); c2->cd(2)->SetLeftMargin(0.2);
  c2->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("N_{clusters}(ITS)");
  hEmpty->GetXaxis()->SetRangeUser(0.0,7.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[5]->Draw("ESAME");
  gCutListHistograms[7]->Draw("ESAME");
  c2->SaveAs("NITSClusters.gif");

  TCanvas *c3 = new TCanvas("c3","Chi2 per ITS cluster",100,100,650,350);
  c3->SetFillColor(10); c3->GetFrame()->SetFillColor(10);
  c3->SetHighLightColor(10); c3->Divide(2,1);
  c3->cd(1)->SetBottomMargin(0.2); c3->cd(1)->SetLeftMargin(0.2);
  c3->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#chi^{2}/N_{clusters}(ITS)");
  hEmpty->GetXaxis()->SetRangeUser(0.0,20.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[8]->Draw("ESAME");
  gCutListHistograms[10]->Draw("ESAME");
  c3->cd(2)->SetBottomMargin(0.2); c3->cd(2)->SetLeftMargin(0.2);
  c3->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#chi^{2}/N_{clusters}(ITS)");
  hEmpty->GetXaxis()->SetRangeUser(0.0,20.0);
  hEmpty->SetTitle("AntiPpotons");
  hEmpty->DrawCopy();
  gCutListHistograms[9]->Draw("ESAME");
  gCutListHistograms[11]->Draw("ESAME");
  c3->SaveAs("Chi2PerITSCluster.gif");

  TCanvas *c4 = new TCanvas("c4","Constrain chi2 - vertex",150,150,650,350);
  c4->SetFillColor(10); c4->GetFrame()->SetFillColor(10);
  c4->SetHighLightColor(10); c4->Divide(2,1);
  c4->cd(1)->SetBottomMargin(0.2); c4->cd(1)->SetLeftMargin(0.2);
  c4->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("log_{10}(#chi^{2}) (vertex)");
  hEmpty->GetXaxis()->SetRangeUser(-10.0,10.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[12]->Draw("ESAME");
  gCutListHistograms[14]->Draw("ESAME");
  c4->cd(2)->SetBottomMargin(0.2); c4->cd(2)->SetLeftMargin(0.2);
  c4->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("log_{10}(#chi^{2}) (vertex )");
  hEmpty->GetXaxis()->SetRangeUser(-10.0,10.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[13]->Draw("ESAME");
  gCutListHistograms[15]->Draw("ESAME");
  c4->SaveAs("ConstrainChi2Vertex.gif");

  TCanvas *c5 = new TCanvas("c5","Number of TPC Clusters",200,200,650,350);
  c5->SetFillColor(10); c5->GetFrame()->SetFillColor(10);
  c5->SetHighLightColor(10); c5->Divide(2,1);
  c5->cd(1)->SetBottomMargin(0.2); c5->cd(1)->SetLeftMargin(0.2);
  c5->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("N_{clusters}(TPC");
  hEmpty->GetXaxis()->SetRangeUser(0.0,200.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[16]->Draw("ESAME");
  gCutListHistograms[18]->Draw("ESAME");
  c5->cd(2)->SetBottomMargin(0.2); c5->cd(2)->SetLeftMargin(0.2);
  c5->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("N_{clusters}(TPC");
  hEmpty->GetXaxis()->SetRangeUser(0.0,200.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[17]->Draw("ESAME");
  gCutListHistograms[19]->Draw("ESAME");
  c5->SaveAs("NTPCClusters.gif");

  TCanvas *c6 = new TCanvas("c6","Chi2 per TPC cluster",250,250,650,350);
  c6->SetFillColor(10); c6->GetFrame()->SetFillColor(10);
  c6->SetHighLightColor(10); c6->Divide(2,1);
  c6->cd(1)->SetBottomMargin(0.2); c6->cd(1)->SetLeftMargin(0.2);
  c6->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#chi^{2}/N_{clusters}(TPC)");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[20]->Draw("ESAME");
  gCutListHistograms[22]->Draw("ESAME");
  c6->cd(2)->SetBottomMargin(0.2); c6->cd(2)->SetLeftMargin(0.2);
  c6->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#chi^{2}/N_{clusters}(TPC)");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[21]->Draw("ESAME");
  gCutListHistograms[23]->Draw("ESAME");
  c6->SaveAs("Chi2PerTPCCluster.gif");

  TCanvas *c7 = new TCanvas("c7","Covariance matrix 11",300,300,650,350);
  c7->SetFillColor(10); c7->GetFrame()->SetFillColor(10);
  c7->SetHighLightColor(10); c7->Divide(2,1);
  c7->cd(1)->SetBottomMargin(0.2); c7->cd(1)->SetLeftMargin(0.2);
  c7->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#sigma_{y} [cm]");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[24]->Draw("ESAME");
  gCutListHistograms[26]->Draw("ESAME");
  c7->cd(2)->SetBottomMargin(0.2); c7->cd(2)->SetLeftMargin(0.2);
  c7->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#sigma_{y} [cm]");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[25]->Draw("ESAME");
  gCutListHistograms[27]->Draw("ESAME");
  c7->SaveAs("Cov11.gif");

  TCanvas *c8 = new TCanvas("c8","Covariance matrix 22",350,350,650,350);
  c8->SetFillColor(10); c8->GetFrame()->SetFillColor(10);
  c8->SetHighLightColor(10); c8->Divide(2,1);
  c8->cd(1)->SetBottomMargin(0.2); c8->cd(1)->SetLeftMargin(0.2);
  c8->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#sigma_{z} [cm]");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[28]->Draw("ESAME");
  gCutListHistograms[30]->Draw("ESAME");
  c8->cd(2)->SetBottomMargin(0.2); c8->cd(2)->SetLeftMargin(0.2);
  c8->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#sigma_{z} [cm]");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[29]->Draw("ESAME");
  gCutListHistograms[31]->Draw("ESAME");
  c8->SaveAs("Cov22.gif");

  TCanvas *c9 = new TCanvas("c9","Covariance matrix 33",400,400,650,350);
  c9->SetFillColor(10); c9->GetFrame()->SetFillColor(10);
  c9->SetHighLightColor(10); c9->Divide(2,1);
  c9->cd(1)->SetBottomMargin(0.2); c9->cd(1)->SetLeftMargin(0.2);
  c9->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#sigma_{sin(#phi)}");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[32]->Draw("ESAME");
  gCutListHistograms[34]->Draw("ESAME");
  c9->cd(2)->SetBottomMargin(0.2); c9->cd(2)->SetLeftMargin(0.2);
  c9->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#sigma_{sin(#phi)}");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[33]->Draw("ESAME");
  gCutListHistograms[35]->Draw("ESAME");
  c9->SaveAs("Cov33.gif");

  TCanvas *c10 = new TCanvas("c10","Covariance matrix 44",450,450,650,350);
  c10->SetFillColor(10); c10->GetFrame()->SetFillColor(10);
  c10->SetHighLightColor(10); c10->Divide(2,1);
  c10->cd(1)->SetBottomMargin(0.2); c10->cd(1)->SetLeftMargin(0.2);
  c10->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#sigma_{tan(#lambda)}");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[36]->Draw("ESAME");
  gCutListHistograms[38]->Draw("ESAME");
  c10->cd(2)->SetBottomMargin(0.2); c10->cd(2)->SetLeftMargin(0.2);
  c10->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#sigma_{tan(#lambda)}");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[37]->Draw("ESAME");
  gCutListHistograms[39]->Draw("ESAME");
  c10->SaveAs("Cov44.gif");

  TCanvas *c11 = new TCanvas("c11","Covariance matrix 55",500,500,650,350);
  c11->SetFillColor(10); c11->GetFrame()->SetFillColor(10);
  c11->SetHighLightColor(10); c11->Divide(2,1);
  c11->cd(1)->SetBottomMargin(0.2); c11->cd(1)->SetLeftMargin(0.2);
  c11->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#sigma_{1/P_{T}} [GeV/c]^{-1}");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[40]->Draw("ESAME");
  gCutListHistograms[42]->Draw("ESAME");
  c11->cd(2)->SetBottomMargin(0.2); c11->cd(2)->SetLeftMargin(0.2);
  c11->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("#sigma_{1/P_{T}} [GeV/c]^{-1}");
  hEmpty->GetXaxis()->SetRangeUser(0.0,4.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[41]->Draw("ESAME");
  gCutListHistograms[43]->Draw("ESAME");
  c11->SaveAs("Cov55.gif");

  TCanvas *c12 = new TCanvas("c12","Number of TPC points (dE/dx)",550,550,650,350);
  c12->SetFillColor(10); c12->GetFrame()->SetFillColor(10);
  c12->SetHighLightColor(10); c12->Divide(2,1);
  c12->cd(1)->SetBottomMargin(0.2); c12->cd(1)->SetLeftMargin(0.2);
  c12->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("N_{points} (TPC-dE/dx");
  hEmpty->GetXaxis()->SetRangeUser(0.0,200.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gCutListHistograms[56]->Draw("ESAME");
  gCutListHistograms[58]->Draw("ESAME");
  c12->cd(2)->SetBottomMargin(0.2); c12->cd(2)->SetLeftMargin(0.2);
  c12->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("N_{points} (TPC-dE/dx");
  hEmpty->GetXaxis()->SetRangeUser(0.0,200.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gCutListHistograms[57]->Draw("ESAME");
  gCutListHistograms[59]->Draw("ESAME");
  c12->SaveAs("Npoints-TPCdEdx.gif");


  //DCA cut parameters  
  TCanvas *c13 = new TCanvas("c13","DCA xy",600,600,650,350);
  c13->SetFillColor(10); c13->GetFrame()->SetFillColor(10);
  c13->SetHighLightColor(10); c13->Divide(2,1);
  c13->cd(1)->SetBottomMargin(0.2); c13->cd(1)->SetLeftMargin(0.2);
  c13->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("DCA_{xy} [cm]");
  hEmpty->GetXaxis()->SetRangeUser(0.0,20.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gDCAListHistograms[2]->Draw("ESAME");
  //gDCAListHistograms[2]->Fit("gDCA","","esame",0.1,12);
  gDCAListHistograms[0]->Draw("ESAME");
  c13->cd(2)->SetBottomMargin(0.15); c13->cd(2)->SetLeftMargin(0.15);
  c13->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("DCA_{xy} [cm]");
  hEmpty->GetXaxis()->SetRangeUser(0.0,20.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gDCAListHistograms[1]->Draw("ESAME");
  gDCAListHistograms[3]->Draw("ESAME");
  c13->SaveAs("DCAxy.gif");

  TCanvas *c14 = new TCanvas("c14","DCA z",650,650,650,350);
  c14->SetFillColor(10); c14->GetFrame()->SetFillColor(10);
  c14->SetHighLightColor(10); c14->Divide(2,1);
  c14->cd(1)->SetBottomMargin(0.2); c14->cd(1)->SetLeftMargin(0.2);
  c14->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("DCA_{z} [cm]");
  hEmpty->GetXaxis()->SetRangeUser(0.0,20.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gDCAListHistograms[4]->Draw("ESAME");
  gDCAListHistograms[6]->Draw("ESAME");
  //gDCAListHistograms[6]->Fit("gDCA","","esame",0.1,12);
  c14->cd(2)->SetBottomMargin(0.15); c14->cd(2)->SetLeftMargin(0.15);
  c14->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("DCA_{z} [cm]");
  hEmpty->GetXaxis()->SetRangeUser(0.0,20.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gDCAListHistograms[5]->Draw("ESAME");
  gDCAListHistograms[7]->Draw("ESAME");
  c14->SaveAs("DCAz.gif");

  TCanvas *c15 = new TCanvas("c15","Sigma to vertex",700,700,650,350);
  c15->SetFillColor(10); c15->GetFrame()->SetFillColor(10);
  c15->SetHighLightColor(10); c15->Divide(2,1);
  c15->cd(1)->SetBottomMargin(0.2); c15->cd(1)->SetLeftMargin(0.2);
  c15->cd(1)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("N_{#sigma}(Vertex)");
  hEmpty->GetXaxis()->SetRangeUser(0.0,7.0);
  hEmpty->SetTitle("Protons");
  hEmpty->DrawCopy();
  gDCAListHistograms[8]->DrawCopy("ESAME");
  gDCAListHistograms[10]->DrawCopy("ESAME");
  c15->cd(2)->SetBottomMargin(0.15); c15->cd(2)->SetLeftMargin(0.15);
  c15->cd(2)->SetLogy();
  hEmpty->GetXaxis()->SetTitle("N_{#sigma}(Vertex)");
  hEmpty->GetXaxis()->SetRangeUser(0.0,7.0);
  hEmpty->SetTitle("Antiprotons");
  hEmpty->DrawCopy();
  gDCAListHistograms[9]->DrawCopy("ESAME");
  gDCAListHistograms[11]->DrawCopy("ESAME");
  c15->SaveAs("NSigmaToVertex.gif");


  f->Close();
}
