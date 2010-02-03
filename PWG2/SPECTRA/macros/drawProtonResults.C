void drawProtonResults(const char* analysisOutput = 0x0,
		       Bool_t kShowResults = kTRUE,
		       Bool_t kShowQAPlots = kFALSE,
		       Bool_t kMC = kFALSE) {
  //Macro to visualize the proton ratio results
  //It also visualizes the QA plots
  gStyle->SetPalette(1,0);
  if(!analysisOutput)
    Error("drawProtonResults::The analysis output was not defined!!!");
  if(kShowResults) drawResults(analysisOutput);
  if(kShowQAPlots) drawQAPlots(analysisOutput, kMC);
}

//___________________________________________________//
void drawResults(const char* analysisOutput) {
  //Draws the main results from the ratio analysis
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG2spectra.so");

  //Create the AliProtonAnalysis object
  AliProtonAnalysis *analysis = new AliProtonAnalysis();
  analysis->ReadFromFile(analysisOutput);
  TH1F *gHistEventStats = dynamic_cast<TH1F *>(analysis->GetEventStatistics());
  TH2D *gHistYPtProtons = dynamic_cast<TH2D *>(analysis->GetProtonYPtHistogram());
  TH2D *gHistYPtAntiProtons = dynamic_cast<TH2D *>(analysis->GetAntiProtonYPtHistogram());
  TH1D *gHistYProtons = dynamic_cast<TH1D *>(analysis->GetProtonYHistogram());
  TH1D *gHistYAntiProtons = dynamic_cast<TH1D *>(analysis->GetAntiProtonYHistogram());
  TH1D *gHistPtProtons = dynamic_cast<TH1D *>(analysis->GetProtonPtHistogram());
  TH1D *gHistPtAntiProtons = dynamic_cast<TH1D *>(analysis->GetAntiProtonPtHistogram());
  TH1D *gHistYRatio = dynamic_cast<TH1D *>(analysis->GetYRatioHistogram());
  TH1D *gHistPtRatio = dynamic_cast<TH1D *>(analysis->GetPtRatioHistogram());

  //==================================================================//
  TH2F *hEmptyRatio = new TH2F("hEmptyRatio",";;#bar{p}/p",100,-1.1,1.1,100,0.1,1.1);
  hEmptyRatio->SetStats(kFALSE);
  hEmptyRatio->GetXaxis()->SetNdivisions(10);
  hEmptyRatio->GetYaxis()->SetNdivisions(10);

  TLatex *latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->SetTextColor(2);

  TF1 *fFitFunction = new TF1("fFitFunction","[0]",-0.5,0.5);
  fFitFunction->SetParameter(0,0.7);

  TCanvas *c2D = new TCanvas("c2D","eta-pT (anti)protons",0,0,700,400);
  c2D->SetFillColor(10); c2D->SetHighLightColor(10); c2D->Divide(2,1);
  c2D->cd(1); gHistYPtProtons->Draw("col");
  c2D->cd(2); gHistYPtAntiProtons->Draw("col");

  TCanvas *cEventStats = new TCanvas("cEventStats","Event statistics",
				     0,0,500,500);
  cEventStats->SetFillColor(10); cEventStats->SetHighLightColor(10);
  gHistEventStats->Draw();

  TCanvas *cEta = new TCanvas("cEta","Eta",100,0,600,400);
  cEta->SetFillColor(10); cEta->SetHighLightColor(10); cEta->Divide(2,1);
  cEta->cd(1); gHistYProtons->Draw("E"); PrintYields(gHistYProtons);
  cEta->cd(2); gHistYAntiProtons->Draw("E"); PrintYields(gHistYAntiProtons);

  TCanvas *cPt = new TCanvas("cPt","Pt",100,200,600,400);
  cPt->SetFillColor(10); cPt->SetHighLightColor(10); cPt->Divide(2,1);
  cPt->cd(1)->SetLogy(); gHistPtProtons->Draw("E");
  cPt->cd(2)->SetLogy(); gHistPtAntiProtons->Draw("E");

  TCanvas *cRatio = new TCanvas("cRatio","Ratio",300,0,600,400);
  cRatio->SetFillColor(10); cRatio->SetHighLightColor(10); cRatio->Divide(2,1);
  cRatio->cd(1); hEmptyRatio->GetXaxis()->SetTitle("eta"); 
  hEmptyRatio->GetXaxis()->SetRangeUser(-1.0,1.0); 
  hEmptyRatio->DrawCopy(); gHistYRatio->Draw("ESAME");
  gHistYRatio->Fit("fFitFunction","N");
  latex->DrawLatex(-0.1,0.45,"ALICE PRELIMINARY");
  latex->DrawLatex(-0.1,0.4,"p-p: #sqrt{s} = 900 GeV");
  cRatio->cd(2);  hEmptyRatio->GetXaxis()->SetTitle("P_{T} [GeV/c]"); 
  hEmptyRatio->GetXaxis()->SetRangeUser(0.3,1.1); 
  hEmptyRatio->DrawCopy(); gHistPtRatio->Draw("ESAME");
  latex->DrawLatex(0.6,0.45,"ALICE PRELIMINARY");
  latex->DrawLatex(0.6,0.4,"p-p: #sqrt{s} = 900 GeV");


  Printf("==========================================");
  for(Int_t iBin = 1; iBin <= gHistYRatio->GetNbinsX(); iBin++)
    Printf("Eta: %lf - Ratio: %lf - Error: %lf",
	   gHistYRatio->GetBinCenter(iBin),
	   gHistYRatio->GetBinContent(iBin),
	   gHistYRatio->GetBinError(iBin));
  Printf("==========================================");

  Printf("Fit result: %lf - %lf",fFitFunction->GetParameter(0),fFitFunction->GetParError(0));
  analysis->PrintMean(gHistYRatio,0.5);
  Printf("==========================================");
}

//___________________________________________________//
void drawQAPlots(const char* analysisOutput,
		 Bool_t kMC) {
  //Draws the QA plots from the output of the analysis
  //=========================================================//
  //List of cuts
  TFile *fCutFile = TFile::Open("ListOfCuts.root");
  TCanvas *cListOfCuts = dynamic_cast<TCanvas *>(fCutFile->Get("cListOfCuts"));
  cListOfCuts->Draw();

  //=========================================================//
  //QA plots
  TFile *f = TFile::Open(analysisOutput);
  TList *listQA = dynamic_cast<TList *>(f->Get("outputQAList"));
  TList *gListGlobalQA = dynamic_cast<TList *>(listQA->At(0));

  //================QA plots================//
  TList *fQA2DList = dynamic_cast<TList *>(gListGlobalQA->At(0));
  //2D de/dx vs P
  TH2F *gHistdEdxP = dynamic_cast<TH2F *>(fQA2DList->At(0));
  gHistdEdxP->SetStats(kFALSE);
  drawdEdx(gHistdEdxP);
  TH2F *gHistProtonsdEdxP = dynamic_cast<TH2F *>(fQA2DList->At(1));
  gHistProtonsdEdxP->SetStats(kFALSE);

  //Theoretical Bethe-Bloch
  Double_t fAlephParameters[5];
  if(kMC) {
    fAlephParameters[0] = 2.15898e+00/50.;
    fAlephParameters[1] = 1.75295e+01;
    fAlephParameters[2] = 3.40030e-09;
    fAlephParameters[3] = 1.96178e+00;
    fAlephParameters[4] = 3.91720e+00;
  }
  else {
    fAlephParameters[0] = 0.0283086;
    fAlephParameters[1] = 2.63394e+01;
    fAlephParameters[2] = 5.04114e-11;
    fAlephParameters[3] = 2.12543e+00;
    fAlephParameters[4] = 4.88663e+00;
  }

  AliTPCPIDResponse *tpcResponse = new AliTPCPIDResponse();
  tpcResponse->SetBetheBlochParameters(fAlephParameters[0],
				       fAlephParameters[1],
				       fAlephParameters[2],
				       fAlephParameters[3],
				       fAlephParameters[4]);
  const Int_t nEntries = 10000;
  Double_t mom[nEntries];
  Double_t dEdxElectrons[nEntries];
  Double_t dEdxMuons[nEntries];
  Double_t dEdxPions[nEntries];
  Double_t dEdxKaons[nEntries];
  Double_t dEdxProtons[nEntries];
  for(Int_t i = 0; i < nEntries; i++) {
    mom[i] = 0.01 + 0.01*i;
    dEdxElectrons[i] = tpcResponse->GetExpectedSignal(mom[i],0);
    dEdxMuons[i] = tpcResponse->GetExpectedSignal(mom[i],1);
    dEdxPions[i] = tpcResponse->GetExpectedSignal(mom[i],2);
    dEdxKaons[i] = tpcResponse->GetExpectedSignal(mom[i],3);
    dEdxProtons[i] = tpcResponse->GetExpectedSignal(mom[i],4);
  }

  TGraph *grElectrons = new TGraph(nEntries,mom,dEdxElectrons);
  grElectrons->SetName("grElectrons");
  grElectrons->SetLineColor(6); grElectrons->SetLineWidth(2);
  TGraph *grMuons = new TGraph(nEntries,mom,dEdxMuons);
  grMuons->SetLineColor(3); grMuons->SetLineWidth(2);
  grMuons->SetName("grMuons");
  TGraph *grPions = new TGraph(nEntries,mom,dEdxPions);
  grPions->SetLineColor(1); grPions->SetLineWidth(2);
  grPions->SetName("grPions");
  TGraph *grKaons = new TGraph(nEntries,mom,dEdxKaons);
  grKaons->SetLineColor(2); grKaons->SetLineWidth(2);
  grKaons->SetName("grKaons");
  TGraph *grProtons = new TGraph(nEntries,mom,dEdxProtons);
  grProtons->SetLineColor(4); grProtons->SetLineWidth(2);
  grProtons->SetName("grProtons");

  //3D eta-phi-NPoints(dEdx)
  TH3F *gHistEtaPhiTPCdEdxNPoints = dynamic_cast<TH3F *>(fQA2DList->At(2));
  TH2D *gHistEtaPhi = dynamic_cast<TH2D *>gHistEtaPhiTPCdEdxNPoints->Project3D("yx");
  gHistEtaPhi->SetStats(kFALSE);
  TH2D *gHistEtaTPCdEdxNPoints = dynamic_cast<TH2D *>gHistEtaPhiTPCdEdxNPoints->Project3D("zx");
  gHistEtaTPCdEdxNPoints->SetStats(kFALSE);
  TH2D *gHistPhiTPCdEdxNPoints = dynamic_cast<TH2D *>gHistEtaPhiTPCdEdxNPoints->Project3D("zy");
  gHistPhiTPCdEdxNPoints->SetStats(kFALSE);

  //3D eta-phi-NPoints(dEdx): protons
  TH3F *gHistProtonsEtaPhiTPCdEdxNPoints = dynamic_cast<TH3F *>(fQA2DList->At(3));
  TH2D *gHistProtonsEtaPhi = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCdEdxNPoints->Project3D("yx");
  gHistProtonsEtaPhi->SetStats(kFALSE);
  TH2D *gHistProtonsEtaTPCdEdxNPoints = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCdEdxNPoints->Project3D("zx");
  gHistProtonsEtaTPCdEdxNPoints->SetStats(kFALSE);
  TH2D *gHistProtonsPhiTPCdEdxNPoints = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCdEdxNPoints->Project3D("zy");
  gHistProtonsPhiTPCdEdxNPoints->SetStats(kFALSE);

  //3D eta-phi-NPoints
  TH3F *gHistEtaPhiTPCNPoints = dynamic_cast<TH3F *>(fQA2DList->At(4));
  TH2D *gHistEtaPhi = dynamic_cast<TH2D *>gHistEtaPhiTPCNPoints->Project3D("yx");
  gHistEtaPhi->SetStats(kFALSE);
  TH2D *gHistEtaTPCNPoints = dynamic_cast<TH2D *>gHistEtaPhiTPCNPoints->Project3D("zx");
  gHistEtaTPCNPoints->SetStats(kFALSE);
  TH2D *gHistPhiTPCNPoints = dynamic_cast<TH2D *>gHistEtaPhiTPCNPoints->Project3D("zy");
  gHistPhiTPCNPoints->SetStats(kFALSE);

  //3D eta-phi-NPoints: protons
  TH3F *gHistProtonsEtaPhiTPCNPoints = dynamic_cast<TH3F *>(fQA2DList->At(5));
  TH2D *gHistProtonsEtaPhi = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCNPoints->Project3D("yx");
  gHistProtonsEtaPhi->SetStats(kFALSE);
  TH2D *gHistProtonsEtaTPCNPoints = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCNPoints->Project3D("zx");
  gHistProtonsEtaTPCNPoints->SetStats(kFALSE);
  TH2D *gHistProtonsPhiTPCNPoints = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCNPoints->Project3D("zy");
  gHistProtonsPhiTPCNPoints->SetStats(kFALSE);

  //3D pt-phi-NPoints(dEdx)
  TH3F *gHistPtPhiTPCdEdxNPoints = dynamic_cast<TH3F *>(fQA2DList->At(6));
  TH2D *gHistPtPhi = dynamic_cast<TH2D *>gHistPtPhiTPCdEdxNPoints->Project3D("yx");
  gHistPtPhi->SetStats(kFALSE);
  TH2D *gHistPtTPCdEdxNPoints = dynamic_cast<TH2D *>gHistPtPhiTPCdEdxNPoints->Project3D("zx");
  gHistPtTPCdEdxNPoints->SetStats(kFALSE);
  TH2D *gHistPhiTPCdEdxNPoints = dynamic_cast<TH2D *>gHistPtPhiTPCdEdxNPoints->Project3D("zy");
  gHistPhiTPCdEdxNPoints->SetStats(kFALSE);

  //3D pt-phi-NPoints(dEdx): protons
  TH3F *gHistProtonsPtPhiTPCdEdxNPoints = dynamic_cast<TH3F *>(fQA2DList->At(7));
  TH2D *gHistProtonsPtPhi = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCdEdxNPoints->Project3D("yx");
  gHistProtonsPtPhi->SetStats(kFALSE);
  TH2D *gHistProtonsPtTPCdEdxNPoints = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCdEdxNPoints->Project3D("zx");
  gHistProtonsPtTPCdEdxNPoints->SetStats(kFALSE);
  TH2D *gHistProtonsPhiTPCdEdxNPoints = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCdEdxNPoints->Project3D("zy");
  gHistProtonsPhiTPCdEdxNPoints->SetStats(kFALSE);

  //3D pt-phi-NPoints
  TH3F *gHistPtPhiTPCNPoints = dynamic_cast<TH3F *>(fQA2DList->At(8));
  TH2D *gHistPtPhi = dynamic_cast<TH2D *>gHistPtPhiTPCNPoints->Project3D("yx");
  gHistPtPhi->SetStats(kFALSE);
  TH2D *gHistPtTPCNPoints = dynamic_cast<TH2D *>gHistPtPhiTPCNPoints->Project3D("zx");
  gHistPtTPCNPoints->SetStats(kFALSE);
  TH2D *gHistPhiTPCNPoints = dynamic_cast<TH2D *>gHistPtPhiTPCNPoints->Project3D("zy");
  gHistPhiTPCNPoints->SetStats(kFALSE);

  //3D pt-phi-NPoints: protons
  TH3F *gHistProtonsPtPhiTPCNPoints = dynamic_cast<TH3F *>(fQA2DList->At(9));
  TH2D *gHistProtonsPtPhi = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCNPoints->Project3D("yx");
  gHistProtonsPtPhi->SetStats(kFALSE);
  TH2D *gHistProtonsPtTPCNPoints = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCNPoints->Project3D("zx");
  gHistProtonsPtTPCNPoints->SetStats(kFALSE);
  TH2D *gHistProtonsPhiTPCNPoints = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCNPoints->Project3D("zy");
  gHistProtonsPhiTPCNPoints->SetStats(kFALSE);

  //2D eta-phi- accepted protons & antiprotons
  TH2F *gHistEtaPhiProtons = dynamic_cast<TH2F *>(fQA2DList->At(10));
  gHistEtaPhiProtons->SetStats(kFALSE);
  TH2F *gHistEtaPhiAntiProtons = dynamic_cast<TH2F *>(fQA2DList->At(11));
  gHistEtaPhiAntiProtons->SetStats(kFALSE);

  //2D dca vs pT - accepted protons & antiprotons
  TH2F *gHistDCAxyPtProtons = dynamic_cast<TH2F *>(fQA2DList->At(12));
  gHistDCAxyPtProtons->SetStats(kFALSE);
  TH2F *gHistDCAzPtProtons = dynamic_cast<TH2F *>(fQA2DList->At(13));
  gHistDCAzPtProtons->SetStats(kFALSE);
  TH2F *gHistDCAxyPtAntiProtons = dynamic_cast<TH2F *>(fQA2DList->At(14));
  gHistDCAxyPtAntiProtons->SetStats(kFALSE);
  TH2F *gHistDCAzPtAntiProtons = dynamic_cast<TH2F *>(fQA2DList->At(15));
  gHistDCAzPtAntiProtons->SetStats(kFALSE);

  //__________________________________________________//
  TH2F *hEmptydEdx = new TH2F("hEmptydEdx",
			      "TPC dE/dx parametrization;P[GeV/c];dE/dx [a.u]",
			      100,0.01,110.,100,30,1000);
  hEmptydEdx->SetStats(kFALSE);

  TLatex *latex = new TLatex();
  latex->SetTextSize(0.035);

  TCanvas *cdEdx = new TCanvas("cdEdx","dE/dx (TPC)",0,0,700,400);
  cdEdx->SetFillColor(10); cdEdx->SetHighLightColor(10); cdEdx->Divide(2,1);
  cdEdx->cd(1)->SetLogx(); hEmptydEdx->DrawCopy();
  gHistdEdxP->Draw("colsame");
  grElectrons->Draw("LSAME"); latex->SetTextColor(6); latex->DrawLatex(0.02,55,"e");
  grMuons->Draw("LSAME"); latex->SetTextColor(3); latex->DrawLatex(0.02,400,"#mu");
  grPions->Draw("LSAME"); latex->SetTextColor(1); latex->DrawLatex(0.05,400,"#pi");
  grKaons->Draw("LSAME"); latex->SetTextColor(2); latex->DrawLatex(0.17,400,"K");
  grProtons->Draw("LSAME"); latex->SetTextColor(4); latex->DrawLatex(0.35,400,"p");
  
  cdEdx->cd(2)->SetLogx(); gHistProtonsdEdxP->Draw("col");

  TCanvas *cEtaPhi = new TCanvas("cEtaPhi",
				 "eta-phi",
				 0,0,700,400);
  cEtaPhi->SetFillColor(10); 
  cEtaPhi->SetHighLightColor(10); cEtaPhi->Divide(2,1);
  cEtaPhi->cd(1); 
  gHistEtaPhiProtons->SetTitle("Accepted protons - eta vs phi");
  gHistEtaPhiProtons->Draw("colz");
  cEtaPhi->cd(2); 
  gHistEtaPhiAntiProtons->SetTitle("Accepted antiprotons - eta vs phi");
  gHistEtaPhiAntiProtons->Draw("colz");

  TCanvas *cDCAPt = new TCanvas("cDCAPt","pT-dca",0,0,700,700);
  cDCAPt->SetFillColor(10); 
  cDCAPt->SetHighLightColor(10); cDCAPt->Divide(2,2);
  cDCAPt->cd(1); 
  gHistDCAxyPtProtons->SetTitle("Accepted protons - dca(xy) vs Pt");
  gHistDCAxyPtProtons->Draw("colz");
  cDCAPt->cd(2); 
  gHistDCAzPtProtons->SetTitle("Accepted protons - dca(z) vs Pt");
  gHistDCAzPtProtons->Draw("colz");
  cDCAPt->cd(3); 
  gHistDCAxyPtAntiProtons->SetTitle("Accepted antiprotons - dca(xy) vs Pt");
  gHistDCAxyPtAntiProtons->Draw("colz");
  cDCAPt->cd(4); 
  gHistDCAzPtAntiProtons->SetTitle("Accepted antiprotons - dca(z) vs Pt");
  gHistDCAzPtAntiProtons->Draw("colz");

  /*TCanvas *cEtaPhiNPointsdEdx = new TCanvas("cEtaPhiNPointsdEdx",
					    "eta-phi-NPoints(dE/dx)",
					    0,0,900,600);
  cEtaPhiNPointsdEdx->SetFillColor(10); 
  cEtaPhiNPointsdEdx->SetHighLightColor(10); cEtaPhiNPointsdEdx->Divide(3,2);
  cEtaPhiNPointsdEdx->cd(1); gHistEtaPhi->Draw("col");
  cEtaPhiNPointsdEdx->cd(2); gHistEtaTPCdEdxNPoints->Draw("col");
  cEtaPhiNPointsdEdx->cd(3); gHistPhiTPCdEdxNPoints->Draw("col");
  cEtaPhiNPointsdEdx->cd(4); gHistProtonsEtaPhi->Draw("col");
  cEtaPhiNPointsdEdx->cd(5); gHistProtonsEtaTPCdEdxNPoints->Draw("col");
  cEtaPhiNPointsdEdx->cd(6); gHistProtonsPhiTPCdEdxNPoints->Draw("col");

  TCanvas *cEtaPhiNPoints = new TCanvas("cEtaPhiNPoints",
					"eta-phi-NPoints",
					0,0,900,600);
  cEtaPhiNPoints->SetFillColor(10); 
  cEtaPhiNPoints->SetHighLightColor(10); cEtaPhiNPoints->Divide(3,2);
  cEtaPhiNPoints->cd(1); gHistEtaPhi->Draw("col");
  cEtaPhiNPoints->cd(2); gHistEtaTPCNPoints->Draw("col");
  cEtaPhiNPoints->cd(3); gHistPhiTPCNPoints->Draw("col");
  cEtaPhiNPoints->cd(4); gHistProtonsEtaPhi->Draw("col");
  cEtaPhiNPoints->cd(5); gHistProtonsEtaTPCNPoints->Draw("col");
  cEtaPhiNPoints->cd(6); gHistProtonsPhiTPCNPoints->Draw("col");

  TCanvas *cPtPhiNPointsdEdx = new TCanvas("cPtPhiNPointsdEdx",
					   "pt-phi-NPoints(dE/dx)",
					   0,0,900,600);
  cPtPhiNPointsdEdx->SetFillColor(10); 
  cPtPhiNPointsdEdx->SetHighLightColor(10); cPtPhiNPointsdEdx->Divide(3,2);
  cPtPhiNPointsdEdx->cd(1); gHistPtPhi->Draw("col");
  cPtPhiNPointsdEdx->cd(2); gHistPtTPCdEdxNPoints->Draw("col");
  cPtPhiNPointsdEdx->cd(3); gHistPhiTPCdEdxNPoints->Draw("col");
  cPtPhiNPointsdEdx->cd(4); gHistProtonsPtPhi->Draw("col");
  cPtPhiNPointsdEdx->cd(5); gHistProtonsPtTPCdEdxNPoints->Draw("col");
  cPtPhiNPointsdEdx->cd(6); gHistProtonsPhiTPCdEdxNPoints->Draw("col");

  TCanvas *cPtPhiNPoints = new TCanvas("cPtPhiNPoints",
				       "pt-phi-NPoints",
				       0,0,900,600);
  cPtPhiNPoints->SetFillColor(10); 
  cPtPhiNPoints->SetHighLightColor(10); cPtPhiNPoints->Divide(3,2);
  cPtPhiNPoints->cd(1); gHistPtPhi->Draw("col");
  cPtPhiNPoints->cd(2); gHistPtTPCNPoints->Draw("col");
  cPtPhiNPoints->cd(3); gHistPhiTPCNPoints->Draw("col");
  cPtPhiNPoints->cd(4); gHistProtonsPtPhi->Draw("col");
  cPtPhiNPoints->cd(5); gHistProtonsPtTPCNPoints->Draw("col");
  cPtPhiNPoints->cd(6); gHistProtonsPhiTPCNPoints->Draw("col");*/

  //Accepted protons
  TList *fQAProtonsAcceptedList = dynamic_cast<TList *>(gListGlobalQA->At(1));
  TH1F *gProtonsITSClustersPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(0));
  TH1F *gProtonsChi2PerClusterITSPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(1));
  TH1F *gProtonsTPCClustersPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(2));
  TH1F *gProtonsChi2PerClusterTPCPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(3));
  TH1F *gProtonsExtCov11Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(4));
  TH1F *gProtonsExtCov22Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(5));
  TH1F *gProtonsExtCov33Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(6));
  TH1F *gProtonsExtCov44Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(7));
  TH1F *gProtonsExtCov55Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(8));
  TH1F *gProtonsSigmaToVertexPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(9));
  TH1F *gProtonsSigmaToVertexTPCPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(10));
  TH1F *gProtonsDCAXYPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(11));
  TH1F *gProtonsDCAXYTPCPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(12));
  TH1F *gProtonsDCAZPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(13));
  TH1F *gProtonsDCAZTPCPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(14));
  TH1F *gProtonsConstrainChi2Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(15));
  TH1F *gProtonsITSRefitPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(16));
  TH1F *gProtonsTPCRefitPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(17));
  TH1F *gProtonsESDpidPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(18));
  TH1F *gProtonsTPCpidPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(19));
  TH1F *gProtonsPointOnITSLayer1Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(20));
  TH1F *gProtonsPointOnITSLayer2Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(21));
  TH1F *gProtonsPointOnITSLayer3Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(22));
  TH1F *gProtonsPointOnITSLayer4Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(23));
  TH1F *gProtonsPointOnITSLayer5Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(24));
  TH1F *gProtonsPointOnITSLayer6Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(25));
  TH1F *gProtonsNumberOfTPCdEdxPointsPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(26));

  //Rejected protons
  TList *fQAProtonsRejectedList = dynamic_cast<TList *>(gListGlobalQA->At(2));
  TH1F *gProtonsITSClustersReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(0));
  TH1F *gProtonsChi2PerClusterITSReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(1));
  TH1F *gProtonsTPCClustersReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(2));
  TH1F *gProtonsChi2PerClusterTPCReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(3));
  TH1F *gProtonsExtCov11Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(4));
  TH1F *gProtonsExtCov22Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(5));
  TH1F *gProtonsExtCov33Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(6));
  TH1F *gProtonsExtCov44Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(7));
  TH1F *gProtonsExtCov55Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(8));
  TH1F *gProtonsSigmaToVertexReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(9));
  TH1F *gProtonsSigmaToVertexTPCReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(10));
  TH1F *gProtonsDCAXYReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(11));
  TH1F *gProtonsDCAXYTPCReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(12));
  TH1F *gProtonsDCAZReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(13));
  TH1F *gProtonsDCAZTPCReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(14));
  TH1F *gProtonsConstrainChi2Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(15));
  TH1F *gProtonsITSRefitReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(16));
  TH1F *gProtonsTPCRefitReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(17));
  TH1F *gProtonsESDpidReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(18));
  TH1F *gProtonsTPCpidReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(19));
  TH1F *gProtonsPointOnITSLayer1Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(20));
  TH1F *gProtonsPointOnITSLayer2Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(21));
  TH1F *gProtonsPointOnITSLayer3Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(22));
  TH1F *gProtonsPointOnITSLayer4Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(23));
  TH1F *gProtonsPointOnITSLayer5Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(24));
  TH1F *gProtonsPointOnITSLayer6Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(25));
  TH1F *gProtonsNumberOfTPCdEdxPointsReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(26));

  //Accepted antiprotons
  TList *fQAAntiProtonsAcceptedList = dynamic_cast<TList *>(gListGlobalQA->At(3));
  TH1F *gAntiProtonsITSClustersPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(0));
  TH1F *gAntiProtonsChi2PerClusterITSPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(1));
  TH1F *gAntiProtonsTPCClustersPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(2));
  TH1F *gAntiProtonsChi2PerClusterTPCPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(3));
  TH1F *gAntiProtonsExtCov11Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(4));
  TH1F *gAntiProtonsExtCov22Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(5));
  TH1F *gAntiProtonsExtCov33Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(6));
  TH1F *gAntiProtonsExtCov44Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(7));
  TH1F *gAntiProtonsExtCov55Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(8));
  TH1F *gAntiProtonsSigmaToVertexPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(9));
  TH1F *gAntiProtonsSigmaToVertexTPCPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(10));
  TH1F *gAntiProtonsDCAXYPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(11));
  TH1F *gAntiProtonsDCAXYTPCPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(12));
  TH1F *gAntiProtonsDCAZPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(13));
  TH1F *gAntiProtonsDCAZTPCPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(14));
  TH1F *gAntiProtonsConstrainChi2Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(15));
  TH1F *gAntiProtonsITSRefitPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(16));
  TH1F *gAntiProtonsTPCRefitPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(17));
  TH1F *gAntiProtonsESDpidPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(18));
  TH1F *gAntiProtonsTPCpidPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(19));
  TH1F *gAntiProtonsPointOnITSLayer1Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(20));
  TH1F *gAntiProtonsPointOnITSLayer2Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(21));
  TH1F *gAntiProtonsPointOnITSLayer3Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(22));
  TH1F *gAntiProtonsPointOnITSLayer4Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(23));
  TH1F *gAntiProtonsPointOnITSLayer5Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(24));
  TH1F *gAntiProtonsPointOnITSLayer6Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(25));
  TH1F *gAntiProtonsNumberOfTPCdEdxPointsPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(26));

  //Rejected antiprotons
  TList *fQAAntiProtonsRejectedList = dynamic_cast<TList *>(gListGlobalQA->At(4));
  TH1F *gAntiProtonsITSClustersReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(0));
  TH1F *gAntiProtonsChi2PerClusterITSReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(1));
  TH1F *gAntiProtonsTPCClustersReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(2));
  TH1F *gAntiProtonsChi2PerClusterTPCReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(3));
  TH1F *gAntiProtonsExtCov11Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(4));
  TH1F *gAntiProtonsExtCov22Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(5));
  TH1F *gAntiProtonsExtCov33Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(6));
  TH1F *gAntiProtonsExtCov44Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(7));
  TH1F *gAntiProtonsExtCov55Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(8));
  TH1F *gAntiProtonsSigmaToVertexReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(9));
  TH1F *gAntiProtonsSigmaToVertexTPCReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(10));
  TH1F *gAntiProtonsDCAXYReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(11));
  TH1F *gAntiProtonsDCAXYTPCReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(12));
  TH1F *gAntiProtonsDCAZReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(13));
  TH1F *gAntiProtonsDCAZTPCReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(14));
  TH1F *gAntiProtonsConstrainChi2Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(15));
  TH1F *gAntiProtonsITSRefitReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(16));
  TH1F *gAntiProtonsTPCRefitReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(17));
  TH1F *gAntiProtonsESDpidReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(18));
  TH1F *gAntiProtonsTPCpidReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(19));
  TH1F *gAntiProtonsPointOnITSLayer1Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(20));
  TH1F *gAntiProtonsPointOnITSLayer2Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(21));
  TH1F *gAntiProtonsPointOnITSLayer3Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(22));
  TH1F *gAntiProtonsPointOnITSLayer4Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(23));
  TH1F *gAntiProtonsPointOnITSLayer5Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(24));
  TH1F *gAntiProtonsPointOnITSLayer6Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(25));
  TH1F *gAntiProtonsNumberOfTPCdEdxPointsReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(26));

  //__________________________________________________//
  TCanvas *c1 = new TCanvas("c1","ITS clusters",0,0,600,400);
  c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->Divide(2,1);
  c1->cd(1); gProtonsITSClustersPass->Draw(); 
  gProtonsITSClustersReject->Draw("same");
  c1->cd(2); gAntiProtonsITSClustersPass->Draw(); 
  gAntiProtonsITSClustersReject->Draw("same");

  TCanvas *c2 = new TCanvas("c2","chi^2 per ITS cluster",0,100,600,400);
  c2->SetFillColor(10); c2->SetHighLightColor(10);
  c2->Divide(2,1);
  c2->cd(1); gProtonsChi2PerClusterITSPass->Draw(); 
  gProtonsChi2PerClusterITSReject->Draw("same");
  c2->cd(2); gAntiProtonsChi2PerClusterITSPass->Draw(); 
  gAntiProtonsChi2PerClusterITSReject->Draw("same");

  TCanvas *c3 = new TCanvas("c3","TPC clusters",0,200,600,400);
  c3->SetFillColor(10); c3->SetHighLightColor(10);
  c3->Divide(2,1);
  c3->cd(1); gProtonsTPCClustersPass->Draw();
  gProtonsTPCClustersReject->Draw("same");
  c3->cd(2); gAntiProtonsTPCClustersPass->Draw();
  gAntiProtonsTPCClustersReject->Draw("same");

  TCanvas *c4 = new TCanvas("c4","chi^2 per TPC cluster",0,300,600,400);
  c4->SetFillColor(10); c4->SetHighLightColor(10);
  c4->Divide(2,1);
  c4->cd(1); gProtonsChi2PerClusterTPCPass->Draw(); 
  gProtonsChi2PerClusterTPCReject->Draw("same");
  c4->cd(2); gAntiProtonsChi2PerClusterTPCPass->Draw(); 
  gAntiProtonsChi2PerClusterTPCReject->Draw("same");

  if(gProtonsExtCov11Pass->GetEntries() != 0) {
    TCanvas *c5 = new TCanvas("c5","Cov11",0,400,600,400);
    c5->SetFillColor(10); c5->SetHighLightColor(10);
    c5->Divide(2,1);
    c5->cd(1)->SetLogy(); gProtonsExtCov11Pass->Draw(); 
    gProtonsExtCov11Reject->Draw("same");
    c5->cd(2)->SetLogy(); gAntiProtonsExtCov11Pass->Draw(); 
    gAntiProtonsExtCov11Reject->Draw("same");
  }

  if(gProtonsExtCov11Pass->GetEntries() != 0) {
    TCanvas *c6 = new TCanvas("c6","Cov22",0,500,600,400);
    c6->SetFillColor(10); c6->SetHighLightColor(10);
    c6->Divide(2,1);
    c6->cd(1)->SetLogy(); gProtonsExtCov22Pass->Draw(); 
    gProtonsExtCov22Reject->Draw("same");
    c6->cd(2)->SetLogy(); gAntiProtonsExtCov22Pass->Draw(); 
    gAntiProtonsExtCov22Reject->Draw("same");
  }

  if(gProtonsExtCov11Pass->GetEntries() != 0) {
    TCanvas *c7 = new TCanvas("c7","Cov33",600,0,600,400);
    c7->SetFillColor(10); c7->SetHighLightColor(10);
    c7->Divide(2,1);
    c7->cd(1)->SetLogy(); gProtonsExtCov33Pass->Draw(); 
    gProtonsExtCov33Reject->Draw("same");
    c7->cd(2)->SetLogy(); gAntiProtonsExtCov33Pass->Draw(); 
    gAntiProtonsExtCov33Reject->Draw("same");
  }

  if(gProtonsExtCov11Pass->GetEntries() != 0) {
    TCanvas *c8 = new TCanvas("c8","Cov44",600,100,600,400);
    c8->SetFillColor(10); c8->SetHighLightColor(10);
    c8->Divide(2,1);
    c8->cd(1)->SetLogy(); gProtonsExtCov44Pass->Draw(); 
    gProtonsExtCov44Reject->Draw("same");
    c8->cd(2)->SetLogy(); gAntiProtonsExtCov44Pass->Draw(); 
    gAntiProtonsExtCov44Reject->Draw("same");
  }

  if(gProtonsExtCov11Pass->GetEntries() != 0) {
    TCanvas *c9 = new TCanvas("c9","Cov55",600,200,600,400);
    c9->SetFillColor(10); c9->SetHighLightColor(10);
    c9->Divide(2,1);
    c9->cd(1)->SetLogy(); gProtonsExtCov55Pass->Draw(); 
    gProtonsExtCov55Reject->Draw("same");
    c9->cd(2)->SetLogy(); gAntiProtonsExtCov55Pass->Draw(); 
    gAntiProtonsExtCov55Reject->Draw("same");
  }

  if(gProtonsSigmaToVertexPass->GetEntries() != 0) {
    TCanvas *c10 = new TCanvas("c10","N-sigma to Vertex",600,300,600,400);
    c10->SetFillColor(10); c10->SetHighLightColor(10);
    c10->Divide(2,1);
    c10->cd(1)->SetLogy(); gProtonsSigmaToVertexPass->Draw(); 
    gProtonsSigmaToVertexReject->Draw("same");
    c10->cd(2)->SetLogy(); gAntiProtonsSigmaToVertexPass->Draw(); 
    gAntiProtonsSigmaToVertexReject->Draw("same");
  }

  if(gProtonsSigmaToVertexTPCPass->GetEntries() != 0) {
    TCanvas *c11 = new TCanvas("c11","N-sigma to Vertex (TPC)",600,400,600,400);
    c11->SetFillColor(10); c11->SetHighLightColor(10);
    c11->Divide(2,1);
    c11->cd(1)->SetLogy(); gProtonsSigmaToVertexTPCPass->Draw(); 
    gProtonsSigmaToVertexTPCReject->Draw("same");
    c11->cd(2)->SetLogy(); gAntiProtonsSigmaToVertexTPCPass->Draw(); 
    gAntiProtonsSigmaToVertexTPCReject->Draw("same");
  }

  if(gProtonsDCAXYPass->GetEntries() != 0) {
    TCanvas *c12 = new TCanvas("c12","dca(xy)",600,500,600,400);
    c12->SetFillColor(10); c12->SetHighLightColor(10);
    c12->Divide(2,1);
    c12->cd(1)->SetLogy(); gProtonsDCAXYPass->Draw(); 
    gProtonsDCAXYReject->Draw("same");
    c12->cd(2)->SetLogy(); gAntiProtonsDCAXYPass->Draw(); 
    gAntiProtonsDCAXYReject->Draw("same");
  }

  if(gProtonsDCAXYTPCPass->GetEntries() != 0) {
    TCanvas *c13 = new TCanvas("c13","dca(xy - TPC)",1200,0,600,400);
    c13->SetFillColor(10); c13->SetHighLightColor(10);
    c13->Divide(2,1);
    c13->cd(1)->SetLogy(); gProtonsDCAXYTPCPass->Draw(); 
    gProtonsDCAXYTPCReject->Draw("same");
    c13->cd(2)->SetLogy(); gAntiProtonsDCAXYTPCPass->Draw(); 
    gAntiProtonsDCAXYTPCReject->Draw("same");
  }

  if(gProtonsDCAZPass->GetEntries() != 0) {
    TCanvas *c14 = new TCanvas("c14","dca(z)",1200,100,600,400);
    c14->SetFillColor(10); c14->SetHighLightColor(10);
    c14->Divide(2,1);
    c14->cd(1)->SetLogy(); gProtonsDCAZPass->Draw(); 
    gProtonsDCAZReject->Draw("same");
    c14->cd(2)->SetLogy(); gAntiProtonsDCAZPass->Draw(); 
    gAntiProtonsDCAZReject->Draw("same");
  }

  if(gProtonsDCAZTPCPass->GetEntries() != 0) {
    TCanvas *c15 = new TCanvas("c15","dca(z - TPC)",1200,200,600,400);
    c15->SetFillColor(10); c15->SetHighLightColor(10);
    c15->Divide(2,1);
    c15->cd(1)->SetLogy(); gProtonsDCAZTPCPass->Draw(); 
    gProtonsDCAZTPCReject->Draw("same");
    c15->cd(2)->SetLogy(); gAntiProtonsDCAZTPCPass->Draw(); 
    gAntiProtonsDCAZTPCReject->Draw("same");
  }

  TCanvas *c16 = new TCanvas("c16","TPC clusters (dE/dx)",1200,300,600,400);
  c16->SetFillColor(10); c16->SetHighLightColor(10);
  c16->Divide(2,1);
  c16->cd(1); gProtonsNumberOfTPCdEdxPointsPass->Draw(); 
  gProtonsNumberOfTPCdEdxPointsReject->Draw("same");
  c16->cd(2); gAntiProtonsNumberOfTPCdEdxPointsPass->Draw(); 
  gAntiProtonsNumberOfTPCdEdxPointsReject->Draw("same");

  //================Vertex QA================//
  TList *gListVertexQA = dynamic_cast<TList *>(listQA->At(1));
  TH1F *gHistVx = dynamic_cast<TH1F *>(gListVertexQA->At(0));
  TH1F *gHistVxAccepted = dynamic_cast<TH1F *>(gListVertexQA->At(1));
  gHistVxAccepted->SetFillColor(10);
  TH1F *gHistVy = dynamic_cast<TH1F *>(gListVertexQA->At(2));
  TH1F *gHistVyAccepted = dynamic_cast<TH1F *>(gListVertexQA->At(3));
  gHistVyAccepted->SetFillColor(10);
  TH1F *gHistVz = dynamic_cast<TH1F *>(gListVertexQA->At(4));
  TH1F *gHistVzAccepted = dynamic_cast<TH1F *>(gListVertexQA->At(5));
  gHistVzAccepted->SetFillColor(10);
  TH1F *gHistNumberOfContributors = dynamic_cast<TH1F *>(gListVertexQA->At(6));
  gHistNumberOfContributors->SetFillColor(10);

  TCanvas *cVertex = new TCanvas("cVertex","Vertex QA",0,0,900,400);
  cVertex->SetFillColor(10); cVertex->SetHighLightColor(10);
  cVertex->Divide(3,1);
  cVertex->cd(1)->SetLogy(); gHistVx->Draw(); gHistVxAccepted->Draw("same");
  cVertex->cd(2)->SetLogy(); gHistVy->Draw(); gHistVyAccepted->Draw("same");
  cVertex->cd(3)->SetLogy(); gHistVz->Draw(); gHistVzAccepted->Draw("same");

  TCanvas *cVertexNContributors = new TCanvas("cVertexNContributors",
					      "Vertex QA",0,0,400,400);
  gHistNumberOfContributors->Draw();
}

//________________________________________________________//
void PrintYields(TH1 *h) {
  Double_t sum = 0.0, error = 0.0;
  for(Int_t iBin = 1; iBin <= h->GetNbinsX(); iBin++) {
    sum += h->GetBinContent(iBin);
    error += TMath::Power(h->GetBinError(iBin),2);
  }
  error = TMath::Sqrt(error);

  Printf("==================================");
  Printf("Histogram: %s",h->GetName());
  Printf("Yields: %lf - %lf",sum,error);
  Printf("==================================");
}

//___________________________________________________//
void drawdEdx(TH2F *gHistdEdxP) {
  //Draws the dE/dx distributions for the different momentum bins
  TString title;
  TH1D *gHist[100];
  Int_t iCounter = 0;
  Double_t binMin = gHistdEdxP->GetXaxis()->GetXmin();
  Double_t binMax = gHistdEdxP->GetXaxis()->GetXmin() + 
    (gHistdEdxP->GetXaxis()->GetXmax() - gHistdEdxP->GetXaxis()->GetXmin())/gHistdEdxP->GetNbinsX();

  TCanvas *c[100];

  for(Int_t iBin = 1; iBin <= gHistdEdxP->GetNbinsX(); iBin++) {
    if((binMax > 0.41)&&(binMin < 0.91)) {
      title = "P: "; title += binMin; title += " - "; 
      title += binMax; title += "GeV/c";
      c[iCounter] = new TCanvas(title.Data(),title.Data(),0,0,500,500);
      c[iCounter]->SetFillColor(10); c[iCounter]->SetHighLightColor(10); 
      gHist[iCounter] = gHistdEdxP->ProjectionY(title.Data(),iBin,iBin);
      gHist[iCounter]->SetTitle(title.Data());
      gHist[iCounter]->SetStats(kFALSE);
      if(gHist[iCounter]->GetEntries() != 0)
	c[iCounter]->SetLogy();
      gHist[iCounter]->Draw();
      Printf("Bin: %d - Pmin: %lf - Pmax: %lf : %s",iBin,binMin,binMax,title.Data());
      iCounter += 1;
    }
    binMin += (gHistdEdxP->GetXaxis()->GetXmax() - gHistdEdxP->GetXaxis()->GetXmin())/gHistdEdxP->GetNbinsX();
    binMax += (gHistdEdxP->GetXaxis()->GetXmax() - gHistdEdxP->GetXaxis()->GetXmin())/gHistdEdxP->GetNbinsX();
  }
  
}
