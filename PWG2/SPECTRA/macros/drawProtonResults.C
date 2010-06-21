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

  TList *gYRatioInPtBinsList = dynamic_cast<TList *>(analysis->GetYRatioHistogramsInPtBins());
  drawRatioInPtBins(gYRatioInPtBinsList);

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

  TCanvas *c2D = new TCanvas("c2D","eta-pT (anti)protons",0,0,900,400);
  c2D->SetFillColor(10); c2D->SetHighLightColor(10); c2D->Divide(3,1);
  c2D->cd(1); gHistYPtProtons->Draw("col");
  c2D->cd(2); gHistYPtAntiProtons->DrawCopy("col");
  c2D->cd(3); gHistYPtAntiProtons->Divide(gHistYPtProtons);
  gHistYPtAntiProtons->SetStats(kFALSE); gHistYPtAntiProtons->DrawCopy("colz");
  TFile *fout = TFile::Open("test.root","recreate");
  gHistYPtAntiProtons->SetName("gHistYPtAntiProtons");
  gHistYPtAntiProtons->Write();
  fout->Close();
  TCanvas *cEventStats = new TCanvas("cEventStats","Event statistics",
				     0,0,500,500);
  cEventStats->SetFillColor(10); cEventStats->SetHighLightColor(10);
  cEventStats->SetLeftMargin(0.15);
  gHistEventStats->GetYaxis()->SetTitleOffset(1.4);
  gHistEventStats->SetStats(kFALSE); gHistEventStats->Draw();

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
  cRatio->cd(1); hEmptyRatio->GetXaxis()->SetTitle("#eta"); 
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

  TFile *fout = TFile::Open("RawRatioPlots.root","recreate");
  gHistYRatio->Write();
  gHistPtRatio->Write();
  fout->Close();

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
void drawRatioInPtBins(TList *gYRatioInPtBinsList) {
  Printf("drawRatioInPtBins:: %d entries",gYRatioInPtBinsList->GetEntries());
  static const Int_t nEntries = gYRatioInPtBinsList->GetEntries();
  TCanvas *cRatioInPtBins[100];
  TH1D *gHistRatioInPtBins[100];
  TString title;
  for(Int_t iEntry = 0; iEntry < gYRatioInPtBinsList->GetEntries(); iEntry++) {
    title = "ratioPtBin"; title += iEntry+1;
    cRatioInPtBins[iEntry] = new TCanvas(title.Data(),
					 title.Data(),
					 0,0,400,400);
    cRatioInPtBins[iEntry]->SetFillColor(10);
    cRatioInPtBins[iEntry]->SetHighLightColor(10);
    gHistRatioInPtBins[iEntry] = dynamic_cast<TH1D *>(gYRatioInPtBinsList->At(iEntry));
    gHistRatioInPtBins[iEntry]->SetStats(kFALSE);
    gHistRatioInPtBins[iEntry]->SetMarkerStyle(20);
    gHistRatioInPtBins[iEntry]->SetMarkerColor(4);
    gHistRatioInPtBins[iEntry]->GetYaxis()->SetRangeUser(0.0,1.4);
    gHistRatioInPtBins[iEntry]->DrawCopy("E");
  }
}

//___________________________________________________//
void drawQAPlots(const char* analysisOutput,
		 Bool_t kMC) {
  //Draws the QA plots from the output of the analysis
  //=========================================================//
  //QA plots
  TFile *f = TFile::Open(analysisOutput);
  //List of cuts
  TCanvas *cListOfCuts = dynamic_cast<TCanvas *>(f->Get("cListOfCuts"));
  if(cListOfCuts)
    cListOfCuts->Draw();

  TList *listQA = dynamic_cast<TList *>(f->Get("outputQAList"));
  TList *gListGlobalQA = dynamic_cast<TList *>(listQA->At(0));

  //================QA plots================//
  TList *fQA2DList = dynamic_cast<TList *>(gListGlobalQA->At(0));
  //2D de/dx vs P
  TH2F *gHistdEdxP = dynamic_cast<TH2F *>(fQA2DList->At(0));
  gHistdEdxP->SetStats(kFALSE);
  drawdEdx(gHistdEdxP,0);
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

  //2D de/dx vs P
  TH2F *gHistZP = dynamic_cast<TH2F *>(fQA2DList->At(2));
  gHistZP->SetStats(kFALSE);
  drawdEdx(gHistZP,1);
  TH2F *gHistProtonsZP = dynamic_cast<TH2F *>(fQA2DList->At(3));
  gHistProtonsZP->SetStats(kFALSE);

  //3D eta-phi-NPoints(dEdx)
  TH3F *gHistEtaPhiTPCdEdxNPoints = dynamic_cast<TH3F *>(fQA2DList->At(4));
  TH2D *gHistEtaPhi = dynamic_cast<TH2D *>gHistEtaPhiTPCdEdxNPoints->Project3D("yx");
  gHistEtaPhi->SetStats(kFALSE);
  TH2D *gHistEtaTPCdEdxNPoints = dynamic_cast<TH2D *>gHistEtaPhiTPCdEdxNPoints->Project3D("zx");
  gHistEtaTPCdEdxNPoints->SetStats(kFALSE);
  TH2D *gHistPhiTPCdEdxNPoints = dynamic_cast<TH2D *>gHistEtaPhiTPCdEdxNPoints->Project3D("zy");
  gHistPhiTPCdEdxNPoints->SetStats(kFALSE);

  //3D eta-phi-NPoints(dEdx): protons
  TH3F *gHistProtonsEtaPhiTPCdEdxNPoints = dynamic_cast<TH3F *>(fQA2DList->At(5));
  TH2D *gHistProtonsEtaPhi = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCdEdxNPoints->Project3D("yx");
  gHistProtonsEtaPhi->SetStats(kFALSE);
  TH2D *gHistProtonsEtaTPCdEdxNPoints = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCdEdxNPoints->Project3D("zx");
  gHistProtonsEtaTPCdEdxNPoints->SetStats(kFALSE);
  TH2D *gHistProtonsPhiTPCdEdxNPoints = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCdEdxNPoints->Project3D("zy");
  gHistProtonsPhiTPCdEdxNPoints->SetStats(kFALSE);

  //3D eta-phi-NPoints
  TH3F *gHistEtaPhiTPCNPoints = dynamic_cast<TH3F *>(fQA2DList->At(6));
  TH2D *gHistEtaPhi = dynamic_cast<TH2D *>gHistEtaPhiTPCNPoints->Project3D("yx");
  gHistEtaPhi->SetStats(kFALSE);
  TH2D *gHistEtaTPCNPoints = dynamic_cast<TH2D *>gHistEtaPhiTPCNPoints->Project3D("zx");
  gHistEtaTPCNPoints->SetStats(kFALSE);
  TH2D *gHistPhiTPCNPoints = dynamic_cast<TH2D *>gHistEtaPhiTPCNPoints->Project3D("zy");
  gHistPhiTPCNPoints->SetStats(kFALSE);

  //3D eta-phi-NPoints: protons
  TH3F *gHistProtonsEtaPhiTPCNPoints = dynamic_cast<TH3F *>(fQA2DList->At(7));
  TH2D *gHistProtonsEtaPhi = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCNPoints->Project3D("yx");
  gHistProtonsEtaPhi->SetStats(kFALSE);
  TH2D *gHistProtonsEtaTPCNPoints = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCNPoints->Project3D("zx");
  gHistProtonsEtaTPCNPoints->SetStats(kFALSE);
  TH2D *gHistProtonsPhiTPCNPoints = dynamic_cast<TH2D *>gHistProtonsEtaPhiTPCNPoints->Project3D("zy");
  gHistProtonsPhiTPCNPoints->SetStats(kFALSE);

  //3D pt-phi-NPoints(dEdx)
  TH3F *gHistPtPhiTPCdEdxNPoints = dynamic_cast<TH3F *>(fQA2DList->At(8));
  TH2D *gHistPtPhi = dynamic_cast<TH2D *>gHistPtPhiTPCdEdxNPoints->Project3D("yx");
  gHistPtPhi->SetStats(kFALSE);
  TH2D *gHistPtTPCdEdxNPoints = dynamic_cast<TH2D *>gHistPtPhiTPCdEdxNPoints->Project3D("zx");
  gHistPtTPCdEdxNPoints->SetStats(kFALSE);
  TH2D *gHistPhiTPCdEdxNPoints = dynamic_cast<TH2D *>gHistPtPhiTPCdEdxNPoints->Project3D("zy");
  gHistPhiTPCdEdxNPoints->SetStats(kFALSE);

  //3D pt-phi-NPoints(dEdx): protons
  TH3F *gHistProtonsPtPhiTPCdEdxNPoints = dynamic_cast<TH3F *>(fQA2DList->At(9));
  TH2D *gHistProtonsPtPhi = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCdEdxNPoints->Project3D("yx");
  gHistProtonsPtPhi->SetStats(kFALSE);
  TH2D *gHistProtonsPtTPCdEdxNPoints = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCdEdxNPoints->Project3D("zx");
  gHistProtonsPtTPCdEdxNPoints->SetStats(kFALSE);
  TH2D *gHistProtonsPhiTPCdEdxNPoints = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCdEdxNPoints->Project3D("zy");
  gHistProtonsPhiTPCdEdxNPoints->SetStats(kFALSE);

  //3D pt-phi-NPoints
  TH3F *gHistPtPhiTPCNPoints = dynamic_cast<TH3F *>(fQA2DList->At(10));
  TH2D *gHistPtPhi = dynamic_cast<TH2D *>gHistPtPhiTPCNPoints->Project3D("yx");
  gHistPtPhi->SetStats(kFALSE);
  TH2D *gHistPtTPCNPoints = dynamic_cast<TH2D *>gHistPtPhiTPCNPoints->Project3D("zx");
  gHistPtTPCNPoints->SetStats(kFALSE);
  TH2D *gHistPhiTPCNPoints = dynamic_cast<TH2D *>gHistPtPhiTPCNPoints->Project3D("zy");
  gHistPhiTPCNPoints->SetStats(kFALSE);

  //3D pt-phi-NPoints: protons
  TH3F *gHistProtonsPtPhiTPCNPoints = dynamic_cast<TH3F *>(fQA2DList->At(11));
  TH2D *gHistProtonsPtPhi = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCNPoints->Project3D("yx");
  gHistProtonsPtPhi->SetStats(kFALSE);
  TH2D *gHistProtonsPtTPCNPoints = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCNPoints->Project3D("zx");
  gHistProtonsPtTPCNPoints->SetStats(kFALSE);
  TH2D *gHistProtonsPhiTPCNPoints = dynamic_cast<TH2D *>gHistProtonsPtPhiTPCNPoints->Project3D("zy");
  gHistProtonsPhiTPCNPoints->SetStats(kFALSE);

  //2D eta-phi- accepted protons & antiprotons
  TH2F *gHistEtaPhiProtons = dynamic_cast<TH2F *>(fQA2DList->At(12));
  gHistEtaPhiProtons->SetStats(kFALSE);
  TH2F *gHistEtaPhiAntiProtons = dynamic_cast<TH2F *>(fQA2DList->At(13));
  gHistEtaPhiAntiProtons->SetStats(kFALSE);

  //2D dca vs pT - accepted protons & antiprotons
  TH3F *gHistDCAxyEtaPtProtons = dynamic_cast<TH3F *>(fQA2DList->At(14));
  gHistDCAxyEtaPtProtons->SetStats(kFALSE);
  TH2D *gHistDCAxyPtProtons = gHistDCAxyEtaPtProtons->Project3D("zy");
  gHistDCAxyPtProtons->SetStats(kFALSE);
  TH3F *gHistDCAzEtaPtProtons = dynamic_cast<TH3F *>(fQA2DList->At(15));
  gHistDCAzEtaPtProtons->SetStats(kFALSE);
  TH2D *gHistDCAzPtProtons = gHistDCAzEtaPtProtons->Project3D("zy");
  gHistDCAzPtProtons->SetStats(kFALSE);
  TH3F *gHistDCAxyEtaPtAntiProtons = dynamic_cast<TH3F *>(fQA2DList->At(16));
  gHistDCAxyEtaPtAntiProtons->SetStats(kFALSE);
  TH2D *gHistDCAxyPtAntiProtons = gHistDCAxyEtaPtAntiProtons->Project3D("zy");
  gHistDCAxyPtAntiProtons->SetStats(kFALSE);
  TH3F *gHistDCAzEtaPtAntiProtons = dynamic_cast<TH3F *>(fQA2DList->At(17));
  gHistDCAzEtaPtAntiProtons->SetStats(kFALSE);
  TH2D *gHistDCAzPtAntiProtons = gHistDCAzEtaPtAntiProtons->Project3D("zy");
  gHistDCAzPtAntiProtons->SetStats(kFALSE);

  //__________________________________________________//
  TH2F *hEmptydEdx = new TH2F("hEmptydEdx",
			      ";p [GeV/c];dE/dx [arb. units]",
			      100,0.08,20.,100,30,1000);
  hEmptydEdx->SetStats(kFALSE);

  TLatex *latex = new TLatex();
  latex->SetTextSize(0.035);

  TCanvas *cdEdx = new TCanvas("cdEdx","dE/dx (TPC)",0,0,700,400);
  cdEdx->SetFillColor(10); cdEdx->SetHighLightColor(10); //cdEdx->Divide(2,1);
  cdEdx->cd(1)->SetLogx(); cdEdx->cd(1)->SetLogy(); hEmptydEdx->DrawCopy();
  gHistdEdxP->Draw("colsame");
  latex->DrawLatex(7.,400,"ALICE");
  latex->DrawLatex(6.,400,"p+p @ #sqrt{s} = 900 GeV");
  grElectrons->Draw("LSAME"); latex->SetTextColor(6); latex->DrawLatex(0.09,55,"e");
  grMuons->Draw("LSAME"); latex->SetTextColor(3); latex->DrawLatex(0.09,400,"#mu");
  grPions->Draw("LSAME"); latex->SetTextColor(1); latex->DrawLatex(0.1,200,"#pi");
  grKaons->Draw("LSAME"); latex->SetTextColor(2); latex->DrawLatex(0.17,400,"K");
  grProtons->Draw("LSAME"); latex->SetTextColor(4); latex->DrawLatex(0.35,400,"p");
  
  //cdEdx->cd(2)->SetLogx(); gHistProtonsdEdxP->Draw("col");

  TCanvas *cZdEdx = new TCanvas("cZdEdx","Normalized dE/dx (TPC)",500,0,700,400);
  cZdEdx->SetFillColor(10); cZdEdx->SetHighLightColor(10); cZdEdx->Divide(2,1);
  cZdEdx->cd(1); gHistZP->Draw("col");
  cZdEdx->cd(2); gHistProtonsZP->Draw("col");

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
  TH1F *gProtonsITSClusterMapPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(27));
  TH1F *gProtonsDCA3DPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(28));

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
  TH1F *gProtonsITSClusterMapReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(27));
  TH1F *gProtonsDCA3DReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(28));

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
  TH1F *gAntiProtonsITSClusterMapPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(27));
  TH1F *gAntiProtonsDCA3DPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(28));

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
  TH1F *gAntiProtonsITSClusterMapReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(27));
  TH1F *gAntiProtonsDCA3DReject = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(28));

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

  TCanvas *c17 = new TCanvas("c17","ITS cluster map",1200,400,600,400);
  c17->SetFillColor(10); c17->SetHighLightColor(10);
  c17->Divide(2,1);
  c17->cd(1); gProtonsITSClusterMapPass->GetYaxis()->SetRangeUser(0,18000);
  gProtonsITSClusterMapPass->Draw(); 
  gProtonsITSClusterMapReject->Draw("same"); 
  c17->cd(2); gAntiProtonsITSClusterMapPass->GetYaxis()->SetRangeUser(0,18000);
  gAntiProtonsITSClusterMapPass->Draw(); 
  gAntiProtonsITSClusterMapReject->Draw("same"); 

  if(gProtonsDCA3DPass->GetEntries() != 0) {  
    TCanvas *c18 = new TCanvas("c18","DCA 3D",1200,500,600,400);
    c18->SetFillColor(10); c18->SetHighLightColor(10);
    c18->Divide(2,1);
    c18->cd(1)->SetLogy(); gProtonsDCA3DPass->Draw();
    gProtonsDCA3DReject->Draw("SAME");
    c18->cd(2)->SetLogy(); gAntiProtonsDCA3DPass->Draw();
    gAntiProtonsDCA3DReject->Draw("SAME");
  }

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
void drawdEdx(TH2F *gHistdEdxP, Int_t iMode) {
  //Draws the dE/dx distributions for the different momentum bins
  //iMode == 0: dEdx vs P
  //iMode == 1: normalized(dEdx) vs P
  TString title;
  TH1D *gHist[100];
  Int_t iCounter = 0;
  Double_t binMin = gHistdEdxP->GetXaxis()->GetXmin();
  Double_t binMax = gHistdEdxP->GetXaxis()->GetXmin() + 
    (gHistdEdxP->GetXaxis()->GetXmax() - gHistdEdxP->GetXaxis()->GetXmin())/gHistdEdxP->GetNbinsX();

  TCanvas *c[100];
  TString canvasTitle;

  for(Int_t iBin = 1; iBin <= gHistdEdxP->GetNbinsX(); iBin++) {
    if((binMax > 0.45)&&(binMin < 1.05)) {
      if(iMode == 0) {
	title = "(dEdx)P: "; title += binMin; title += " - "; 
	title += binMax; title += "GeV/c";
	canvasTitle = "dedxMomentumSlice."; canvasTitle += iCounter;
	canvasTitle += ".gif";
	c[iCounter] = new TCanvas(title.Data(),title.Data(),0,0,500,500);
      }
      if(iMode == 1) {
	title = "(normdEdx)P: "; title += binMin; title += " - "; 
	title += binMax; title += "GeV/c";
	canvasTitle = "normdedxMomentumSlice."; canvasTitle += iCounter;
	canvasTitle += ".gif";
	c[iCounter] = new TCanvas(title.Data(),title.Data(),500,0,500,500);
      }
      c[iCounter]->SetFillColor(10); c[iCounter]->SetHighLightColor(10); 
      gHist[iCounter] = gHistdEdxP->ProjectionY(title.Data(),iBin,iBin+1);
      gHist[iCounter]->SetTitle(title.Data());
      gHist[iCounter]->SetStats(kFALSE);
      if(iMode == 0)
	gHist[iCounter]->GetXaxis()->SetRangeUser(0.0,300.);
      if(iMode == 1)
	gHist[iCounter]->GetXaxis()->SetRangeUser(-2.0,2.0);
      if(gHist[iCounter]->GetEntries() != 0)
	c[iCounter]->SetLogy();
      gHist[iCounter]->Draw("E");
      c[iCounter]->SaveAs(canvasTitle.Data());
      //Printf("Bin: %d - Pmin: %lf - Pmax: %lf : %s",iBin,binMin,binMax,title.Data());
      iCounter += 1;
    }
    binMin += (gHistdEdxP->GetXaxis()->GetXmax() - gHistdEdxP->GetXaxis()->GetXmin())/gHistdEdxP->GetNbinsX();
    binMax += (gHistdEdxP->GetXaxis()->GetXmax() - gHistdEdxP->GetXaxis()->GetXmin())/gHistdEdxP->GetNbinsX();
  }
  
}

//___________________________________________________//
void drawMCvsData(const char* analysisOutputMC,
		  const char* analysisOutputData) {
  //Draws the QA plots from the output of the analysis
  //=========================================================//
  gStyle->SetPalette(1,0);

  //=========================================================//
  //QA MC plots
  TFile *fMC = TFile::Open(analysisOutputMC);
  TList *listQAMC = dynamic_cast<TList *>(fMC->Get("outputQAList"));
  TList *gListGlobalQAMC = dynamic_cast<TList *>(listQAMC->At(0));

  //================QA plots================//
  TList *fQAMC2DList = dynamic_cast<TList *>(gListGlobalQAMC->At(0));

  //Rejected protons
  TList *fQAMCProtonsRejectedList = dynamic_cast<TList *>(gListGlobalQAMC->At(2));
  TH1F *gMCProtonsITSClustersReject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(0));
  TH1F *gMCProtonsTPCClustersReject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(2));
  TH1F *gMCProtonsChi2PerClusterTPCReject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(3));
  TH1F *gMCProtonsExtCov11Reject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(4));
  TH1F *gMCProtonsExtCov22Reject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(5));
  TH1F *gMCProtonsExtCov33Reject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(6));
  TH1F *gMCProtonsExtCov44Reject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(7));
  TH1F *gMCProtonsExtCov55Reject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(8));
  TH1F *gMCProtonsDCAXYReject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(11));
  TH1F *gMCProtonsDCAZReject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(13));
  TH1F *gMCProtonsNumberOfTPCdEdxPointsReject = dynamic_cast<TH1F *>(fQAMCProtonsRejectedList->At(26));

  //Accepted protons
  TList *fQAMCProtonsAcceptedList = dynamic_cast<TList *>(gListGlobalQAMC->At(1));
  TH1F *gMCProtonsITSClustersPass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(0));
  gMCProtonsITSClustersPass->SetStats(kFALSE);
  gMCProtonsITSClustersPass->Add(gMCProtonsITSClustersReject);
  gMCProtonsITSClustersPass->Sumw2();
  gMCProtonsITSClustersPass->Scale(1./gMCProtonsITSClustersPass->Integral(1,gMCProtonsITSClustersPass->GetNbinsX()));
  gMCProtonsITSClustersPass->SetMarkerStyle(20);
  TH1F *gMCProtonsTPCClustersPass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(2));
  gMCProtonsTPCClustersPass->SetStats(kFALSE);
  gMCProtonsTPCClustersPass->Add(gMCProtonsTPCClustersReject);
  gMCProtonsTPCClustersPass->Sumw2();
  gMCProtonsTPCClustersPass->Scale(1./gMCProtonsTPCClustersPass->Integral(1,gMCProtonsTPCClustersPass->GetNbinsX()));
  gMCProtonsTPCClustersPass->SetMarkerStyle(20);
  TH1F *gMCProtonsChi2PerClusterTPCPass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(3));
  gMCProtonsChi2PerClusterTPCPass->SetStats(kFALSE);
  gMCProtonsChi2PerClusterTPCPass->Add(gMCProtonsChi2PerClusterTPCReject);
  gMCProtonsChi2PerClusterTPCPass->Sumw2();
  gMCProtonsChi2PerClusterTPCPass->Scale(1./gMCProtonsChi2PerClusterTPCPass->Integral(1,gMCProtonsChi2PerClusterTPCPass->GetNbinsX()));
  gMCProtonsChi2PerClusterTPCPass->SetMarkerStyle(20);
  TH1F *gMCProtonsExtCov11Pass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(4));
  gMCProtonsExtCov11Pass->SetStats(kFALSE);
  gMCProtonsExtCov11Pass->Add(gMCProtonsExtCov11Reject);
  gMCProtonsExtCov11Pass->Sumw2();
  gMCProtonsExtCov11Pass->Scale(1./gMCProtonsExtCov11Pass->Integral(1,gMCProtonsExtCov11Pass->GetNbinsX()));
  gMCProtonsExtCov11Pass->SetMarkerStyle(20);
  TH1F *gMCProtonsExtCov22Pass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(5));
  gMCProtonsExtCov22Pass->SetStats(kFALSE);
  gMCProtonsExtCov22Pass->Add(gMCProtonsExtCov22Reject);
  gMCProtonsExtCov22Pass->Sumw2();
  gMCProtonsExtCov22Pass->Scale(1./gMCProtonsExtCov22Pass->Integral(1,gMCProtonsExtCov22Pass->GetNbinsX()));
  gMCProtonsExtCov22Pass->SetMarkerStyle(20);
  TH1F *gMCProtonsExtCov33Pass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(6));
  gMCProtonsExtCov33Pass->SetStats(kFALSE);
  gMCProtonsExtCov33Pass->Add(gMCProtonsExtCov33Reject);
  gMCProtonsExtCov33Pass->Sumw2();
  gMCProtonsExtCov33Pass->Scale(1./gMCProtonsExtCov33Pass->Integral(1,gMCProtonsExtCov33Pass->GetNbinsX()));
  gMCProtonsExtCov33Pass->SetMarkerStyle(20);
  TH1F *gMCProtonsExtCov44Pass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(7));
  gMCProtonsExtCov44Pass->SetStats(kFALSE);
  gMCProtonsExtCov44Pass->Add(gMCProtonsExtCov44Reject);
  gMCProtonsExtCov44Pass->Sumw2();
  gMCProtonsExtCov44Pass->Scale(1./gMCProtonsExtCov44Pass->Integral(1,gMCProtonsExtCov44Pass->GetNbinsX()));
  gMCProtonsExtCov44Pass->SetMarkerStyle(20);
  TH1F *gMCProtonsExtCov55Pass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(8));
  gMCProtonsExtCov55Pass->SetStats(kFALSE);
  gMCProtonsExtCov55Pass->Add(gMCProtonsExtCov55Reject);
  gMCProtonsExtCov55Pass->Sumw2();
  gMCProtonsExtCov55Pass->Scale(1./gMCProtonsExtCov55Pass->Integral(1,gMCProtonsExtCov55Pass->GetNbinsX()));
  gMCProtonsExtCov55Pass->SetMarkerStyle(20);
  TH1F *gMCProtonsDCAXYPass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(11));
  gMCProtonsDCAXYPass->SetStats(kFALSE);  
  gMCProtonsDCAXYPass->Add(gMCProtonsDCAXYReject);
  gMCProtonsDCAXYPass->Sumw2();
  gMCProtonsDCAXYPass->Scale(1./gMCProtonsDCAXYPass->Integral(1,gMCProtonsDCAXYPass->GetNbinsX()));
  gMCProtonsDCAXYPass->SetMarkerStyle(20);
  TH1F *gMCProtonsDCAZPass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(13));
  gMCProtonsDCAZPass->SetStats(kFALSE);  
  gMCProtonsDCAZPass->Add(gMCProtonsDCAZReject);
  gMCProtonsDCAZPass->Sumw2();
  gMCProtonsDCAZPass->Scale(1./gMCProtonsDCAZPass->Integral(1,gMCProtonsDCAZPass->GetNbinsX()));
  gMCProtonsDCAZPass->SetMarkerStyle(20);
  TH1F *gMCProtonsNumberOfTPCdEdxPointsPass = dynamic_cast<TH1F *>(fQAMCProtonsAcceptedList->At(26));
  gMCProtonsNumberOfTPCdEdxPointsPass->SetStats(kFALSE);  
  gMCProtonsNumberOfTPCdEdxPointsPass->Add(gMCProtonsNumberOfTPCdEdxPointsReject);
  gMCProtonsNumberOfTPCdEdxPointsPass->Sumw2();
  gMCProtonsNumberOfTPCdEdxPointsPass->Scale(1./gMCProtonsNumberOfTPCdEdxPointsPass->Integral(1,gMCProtonsNumberOfTPCdEdxPointsPass->GetNbinsX()));
  gMCProtonsNumberOfTPCdEdxPointsPass->SetMarkerStyle(20);

  //Rejected antiprotons
  TList *fQAMCAntiProtonsRejectedList = dynamic_cast<TList *>(gListGlobalQAMC->At(4));
  TH1F *gMCAntiProtonsITSClustersReject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(0));
  TH1F *gMCAntiProtonsTPCClustersReject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(2));
  TH1F *gMCAntiProtonsChi2PerClusterTPCReject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(3));
  TH1F *gMCAntiProtonsExtCov11Reject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(4));
  TH1F *gMCAntiProtonsExtCov22Reject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(5));
  TH1F *gMCAntiProtonsExtCov33Reject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(6));
  TH1F *gMCAntiProtonsExtCov44Reject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(7));
  TH1F *gMCAntiProtonsExtCov55Reject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(8));
  TH1F *gMCAntiProtonsDCAXYReject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(11));
  TH1F *gMCAntiProtonsDCAZReject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(13));
  TH1F *gMCAntiProtonsNumberOfTPCdEdxPointsReject = dynamic_cast<TH1F *>(fQAMCAntiProtonsRejectedList->At(26));

  //Accepted protons
  TList *fQAMCAntiProtonsAcceptedList = dynamic_cast<TList *>(gListGlobalQAMC->At(3));
  TH1F *gMCAntiProtonsITSClustersPass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(0));
  gMCAntiProtonsITSClustersPass->SetStats(kFALSE);
  gMCAntiProtonsITSClustersPass->Add(gMCAntiProtonsITSClustersReject);
  gMCAntiProtonsITSClustersPass->Sumw2();
  gMCAntiProtonsITSClustersPass->Scale(1./gMCAntiProtonsITSClustersPass->Integral(1,gMCAntiProtonsITSClustersPass->GetNbinsX()));
  gMCAntiProtonsITSClustersPass->SetMarkerStyle(20);
  TH1F *gMCAntiProtonsTPCClustersPass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(2));
  gMCAntiProtonsTPCClustersPass->SetStats(kFALSE);
  gMCAntiProtonsTPCClustersPass->Add(gMCAntiProtonsTPCClustersReject);
  gMCAntiProtonsTPCClustersPass->Sumw2();
  gMCAntiProtonsTPCClustersPass->Scale(1./gMCAntiProtonsTPCClustersPass->Integral(1,gMCAntiProtonsTPCClustersPass->GetNbinsX()));
  gMCAntiProtonsTPCClustersPass->SetMarkerStyle(20);
  TH1F *gMCAntiProtonsChi2PerClusterTPCPass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(3));
  gMCAntiProtonsChi2PerClusterTPCPass->SetStats(kFALSE);
  gMCAntiProtonsChi2PerClusterTPCPass->Add(gMCAntiProtonsChi2PerClusterTPCReject);
  gMCAntiProtonsChi2PerClusterTPCPass->Sumw2();
  gMCAntiProtonsChi2PerClusterTPCPass->Scale(1./gMCAntiProtonsChi2PerClusterTPCPass->Integral(1,gMCAntiProtonsChi2PerClusterTPCPass->GetNbinsX()));
  gMCAntiProtonsChi2PerClusterTPCPass->SetMarkerStyle(20);
  TH1F *gMCAntiProtonsExtCov11Pass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(4));
  gMCAntiProtonsExtCov11Pass->SetStats(kFALSE);
  gMCAntiProtonsExtCov11Pass->Add(gMCAntiProtonsExtCov11Reject);
  gMCAntiProtonsExtCov11Pass->Sumw2();
  gMCAntiProtonsExtCov11Pass->Scale(1./gMCAntiProtonsExtCov11Pass->Integral(1,gMCAntiProtonsExtCov11Pass->GetNbinsX()));
  gMCAntiProtonsExtCov11Pass->SetMarkerStyle(20);
  TH1F *gMCAntiProtonsExtCov22Pass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(5));
  gMCAntiProtonsExtCov22Pass->SetStats(kFALSE);
  gMCAntiProtonsExtCov22Pass->Add(gMCAntiProtonsExtCov22Reject);
  gMCAntiProtonsExtCov22Pass->Sumw2();
  gMCAntiProtonsExtCov22Pass->Scale(1./gMCAntiProtonsExtCov22Pass->Integral(1,gMCAntiProtonsExtCov22Pass->GetNbinsX()));
  gMCAntiProtonsExtCov22Pass->SetMarkerStyle(20);
  TH1F *gMCAntiProtonsExtCov33Pass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(6));
  gMCAntiProtonsExtCov33Pass->SetStats(kFALSE);
  gMCAntiProtonsExtCov33Pass->Add(gMCAntiProtonsExtCov33Reject);
  gMCAntiProtonsExtCov33Pass->Sumw2();
  gMCAntiProtonsExtCov33Pass->Scale(1./gMCAntiProtonsExtCov33Pass->Integral(1,gMCAntiProtonsExtCov33Pass->GetNbinsX()));
  gMCAntiProtonsExtCov33Pass->SetMarkerStyle(20);
  TH1F *gMCAntiProtonsExtCov44Pass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(7));
  gMCAntiProtonsExtCov44Pass->SetStats(kFALSE);
  gMCAntiProtonsExtCov44Pass->Add(gMCAntiProtonsExtCov44Reject);
  gMCAntiProtonsExtCov44Pass->Sumw2();
  gMCAntiProtonsExtCov44Pass->Scale(1./gMCAntiProtonsExtCov44Pass->Integral(1,gMCAntiProtonsExtCov44Pass->GetNbinsX()));
  gMCAntiProtonsExtCov44Pass->SetMarkerStyle(20);
  TH1F *gMCAntiProtonsExtCov55Pass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(8));
  gMCAntiProtonsExtCov55Pass->SetStats(kFALSE);
  gMCAntiProtonsExtCov55Pass->Add(gMCAntiProtonsExtCov55Reject);
  gMCAntiProtonsExtCov55Pass->Sumw2();
  gMCAntiProtonsExtCov55Pass->Scale(1./gMCAntiProtonsExtCov55Pass->Integral(1,gMCAntiProtonsExtCov55Pass->GetNbinsX()));
  gMCAntiProtonsExtCov55Pass->SetMarkerStyle(20);
  TH1F *gMCAntiProtonsDCAXYPass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(11));
  gMCAntiProtonsDCAXYPass->SetStats(kFALSE);  
  gMCAntiProtonsDCAXYPass->Add(gMCAntiProtonsDCAXYReject);
  gMCAntiProtonsDCAXYPass->Sumw2();
  gMCAntiProtonsDCAXYPass->Scale(1./gMCAntiProtonsDCAXYPass->Integral(1,gMCAntiProtonsDCAXYPass->GetNbinsX()));
  gMCAntiProtonsDCAXYPass->SetMarkerStyle(20);
  TH1F *gMCAntiProtonsDCAZPass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(13));
  gMCAntiProtonsDCAZPass->SetStats(kFALSE);  
  gMCAntiProtonsDCAZPass->Add(gMCAntiProtonsDCAZReject);
  gMCAntiProtonsDCAZPass->Sumw2();
  gMCAntiProtonsDCAZPass->Scale(1./gMCAntiProtonsDCAZPass->Integral(1,gMCAntiProtonsDCAZPass->GetNbinsX()));
  gMCAntiProtonsDCAZPass->SetMarkerStyle(20);
  TH1F *gMCAntiProtonsNumberOfTPCdEdxPointsPass = dynamic_cast<TH1F *>(fQAMCAntiProtonsAcceptedList->At(26));
  gMCAntiProtonsNumberOfTPCdEdxPointsPass->SetStats(kFALSE);  
  gMCAntiProtonsNumberOfTPCdEdxPointsPass->Add(gMCAntiProtonsNumberOfTPCdEdxPointsReject);
  gMCAntiProtonsNumberOfTPCdEdxPointsPass->Sumw2();
  gMCAntiProtonsNumberOfTPCdEdxPointsPass->Scale(1./gMCAntiProtonsNumberOfTPCdEdxPointsPass->Integral(1,gMCAntiProtonsNumberOfTPCdEdxPointsPass->GetNbinsX()));
  gMCAntiProtonsNumberOfTPCdEdxPointsPass->SetMarkerStyle(20);

  //=========================================================//
  //QA data plots
  TFile *fData = TFile::Open(analysisOutputData);
  TList *listQAData = dynamic_cast<TList *>(fData->Get("outputQAList"));
  TList *gListGlobalQAData = dynamic_cast<TList *>(listQAData->At(0));

  //================QA plots================//
  TList *fQAData2DList = dynamic_cast<TList *>(gListGlobalQAData->At(0));

  //Rejected protons
  TList *fQADataProtonsRejectedList = dynamic_cast<TList *>(gListGlobalQAData->At(2));
  TH1F *gDataProtonsITSClustersReject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(0));
  TH1F *gDataProtonsTPCClustersReject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(2));
  TH1F *gDataProtonsChi2PerClusterTPCReject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(3));
  TH1F *gDataProtonsExtCov11Reject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(4));
  TH1F *gDataProtonsExtCov22Reject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(5));
  TH1F *gDataProtonsExtCov33Reject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(6));
  TH1F *gDataProtonsExtCov44Reject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(7));
  TH1F *gDataProtonsExtCov55Reject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(8));
  TH1F *gDataProtonsDCAXYReject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(11));
  TH1F *gDataProtonsDCAZReject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(13));
  TH1F *gDataProtonsNumberOfTPCdEdxPointsReject = dynamic_cast<TH1F *>(fQADataProtonsRejectedList->At(26));

  //Accepted protons
  TList *fQADataProtonsAcceptedList = dynamic_cast<TList *>(gListGlobalQAData->At(1));
  TH1F *gDataProtonsITSClustersPass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(0));
  gDataProtonsITSClustersPass->SetStats(kFALSE);
  gDataProtonsITSClustersPass->Add(gDataProtonsITSClustersReject);
  gDataProtonsITSClustersPass->Sumw2();
  gDataProtonsITSClustersPass->Scale(1./gDataProtonsITSClustersPass->Integral(1,gDataProtonsITSClustersPass->GetNbinsX()));
  gDataProtonsITSClustersPass->SetMarkerStyle(1);
  TH1F *gDataProtonsTPCClustersPass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(2));
  gDataProtonsTPCClustersPass->SetStats(kFALSE);
  gDataProtonsTPCClustersPass->Add(gDataProtonsTPCClustersReject);
  gDataProtonsTPCClustersPass->Sumw2();
  gDataProtonsTPCClustersPass->Scale(1./gDataProtonsTPCClustersPass->Integral(1,gDataProtonsTPCClustersPass->GetNbinsX()));
  gDataProtonsTPCClustersPass->SetMarkerStyle(1);
  TH1F *gDataProtonsChi2PerClusterTPCPass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(3));
  gDataProtonsChi2PerClusterTPCPass->SetStats(kFALSE);
  gDataProtonsChi2PerClusterTPCPass->Add(gDataProtonsChi2PerClusterTPCReject);
  gDataProtonsChi2PerClusterTPCPass->Sumw2();
  gDataProtonsChi2PerClusterTPCPass->Scale(1./gDataProtonsChi2PerClusterTPCPass->Integral(1,gDataProtonsChi2PerClusterTPCPass->GetNbinsX()));
  gDataProtonsChi2PerClusterTPCPass->SetMarkerStyle(1);
  TH1F *gDataProtonsExtCov11Pass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(4));
  gDataProtonsExtCov11Pass->SetStats(kFALSE);
  gDataProtonsExtCov11Pass->Add(gDataProtonsExtCov11Reject);
  gDataProtonsExtCov11Pass->Sumw2();
  gDataProtonsExtCov11Pass->Scale(1./gDataProtonsExtCov11Pass->Integral(1,gDataProtonsExtCov11Pass->GetNbinsX()));
  gDataProtonsExtCov11Pass->SetMarkerStyle(1);
  TH1F *gDataProtonsExtCov22Pass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(5));
  gDataProtonsExtCov22Pass->SetStats(kFALSE);
  gDataProtonsExtCov22Pass->Add(gDataProtonsExtCov22Reject);
  gDataProtonsExtCov22Pass->Sumw2();
  gDataProtonsExtCov22Pass->Scale(1./gDataProtonsExtCov22Pass->Integral(1,gDataProtonsExtCov22Pass->GetNbinsX()));
  gDataProtonsExtCov22Pass->SetMarkerStyle(1);
  TH1F *gDataProtonsExtCov33Pass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(6));
  gDataProtonsExtCov33Pass->SetStats(kFALSE);
  gDataProtonsExtCov33Pass->Add(gDataProtonsExtCov33Reject);
  gDataProtonsExtCov33Pass->Sumw2();
  gDataProtonsExtCov33Pass->Scale(1./gDataProtonsExtCov33Pass->Integral(1,gDataProtonsExtCov33Pass->GetNbinsX()));
  gDataProtonsExtCov33Pass->SetMarkerStyle(1);
  TH1F *gDataProtonsExtCov44Pass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(7));
  gDataProtonsExtCov44Pass->SetStats(kFALSE);
  gDataProtonsExtCov44Pass->Add(gDataProtonsExtCov44Reject);
  gDataProtonsExtCov44Pass->Sumw2();
  gDataProtonsExtCov44Pass->Scale(1./gDataProtonsExtCov44Pass->Integral(1,gDataProtonsExtCov44Pass->GetNbinsX()));
  gDataProtonsExtCov44Pass->SetMarkerStyle(1);
  TH1F *gDataProtonsExtCov55Pass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(8));
  gDataProtonsExtCov55Pass->SetStats(kFALSE);
  gDataProtonsExtCov55Pass->Add(gDataProtonsExtCov55Reject);
  gDataProtonsExtCov55Pass->Sumw2();
  gDataProtonsExtCov55Pass->Scale(1./gDataProtonsExtCov55Pass->Integral(1,gDataProtonsExtCov55Pass->GetNbinsX()));
  gDataProtonsExtCov55Pass->SetMarkerStyle(1);
  TH1F *gDataProtonsDCAXYPass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(11));
  gDataProtonsDCAXYPass->SetStats(kFALSE);  
  gDataProtonsDCAXYPass->Add(gDataProtonsDCAXYReject);
  gDataProtonsDCAXYPass->Sumw2();
  gDataProtonsDCAXYPass->Scale(1./gDataProtonsDCAXYPass->Integral(1,gDataProtonsDCAXYPass->GetNbinsX()));
  gDataProtonsDCAXYPass->SetMarkerStyle(1);
  TH1F *gDataProtonsDCAZPass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(13));
  gDataProtonsDCAZPass->SetStats(kFALSE);  
  gDataProtonsDCAZPass->Add(gDataProtonsDCAZReject);
  gDataProtonsDCAZPass->Sumw2();
  gDataProtonsDCAZPass->Scale(1./gDataProtonsDCAZPass->Integral(1,gDataProtonsDCAZPass->GetNbinsX()));
  gDataProtonsDCAZPass->SetMarkerStyle(1);
  TH1F *gDataProtonsNumberOfTPCdEdxPointsPass = dynamic_cast<TH1F *>(fQADataProtonsAcceptedList->At(26));
  gDataProtonsNumberOfTPCdEdxPointsPass->SetStats(kFALSE);  
  gDataProtonsNumberOfTPCdEdxPointsPass->Add(gDataProtonsNumberOfTPCdEdxPointsReject);
  gDataProtonsNumberOfTPCdEdxPointsPass->Sumw2();
  gDataProtonsNumberOfTPCdEdxPointsPass->Scale(1./gDataProtonsNumberOfTPCdEdxPointsPass->Integral(1,gDataProtonsNumberOfTPCdEdxPointsPass->GetNbinsX()));
  gDataProtonsNumberOfTPCdEdxPointsPass->SetMarkerStyle(1);

  //Rejected antiprotons
  TList *fQADataAntiProtonsRejectedList = dynamic_cast<TList *>(gListGlobalQAData->At(4));
  TH1F *gDataAntiProtonsITSClustersReject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(0));
  TH1F *gDataAntiProtonsTPCClustersReject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(2));
  TH1F *gDataAntiProtonsChi2PerClusterTPCReject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(3));
  TH1F *gDataAntiProtonsExtCov11Reject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(4));
  TH1F *gDataAntiProtonsExtCov22Reject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(5));
  TH1F *gDataAntiProtonsExtCov33Reject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(6));
  TH1F *gDataAntiProtonsExtCov44Reject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(7));
  TH1F *gDataAntiProtonsExtCov55Reject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(8));
  TH1F *gDataAntiProtonsDCAXYReject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(11));
  TH1F *gDataAntiProtonsDCAZReject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(13));
  TH1F *gDataAntiProtonsNumberOfTPCdEdxPointsReject = dynamic_cast<TH1F *>(fQADataAntiProtonsRejectedList->At(26));

  //Accepted protons
  TList *fQADataAntiProtonsAcceptedList = dynamic_cast<TList *>(gListGlobalQAData->At(3));
  TH1F *gDataAntiProtonsITSClustersPass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(0));
  gDataAntiProtonsITSClustersPass->SetStats(kFALSE);
  gDataAntiProtonsITSClustersPass->Add(gDataAntiProtonsITSClustersReject);
  gDataAntiProtonsITSClustersPass->Sumw2();
  gDataAntiProtonsITSClustersPass->Scale(1./gDataAntiProtonsITSClustersPass->Integral(1,gDataAntiProtonsITSClustersPass->GetNbinsX()));
  gDataAntiProtonsITSClustersPass->SetMarkerStyle(1);
  TH1F *gDataAntiProtonsTPCClustersPass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(2));
  gDataAntiProtonsTPCClustersPass->SetStats(kFALSE);
  gDataAntiProtonsTPCClustersPass->Add(gDataAntiProtonsTPCClustersReject);
  gDataAntiProtonsTPCClustersPass->Sumw2();
  gDataAntiProtonsTPCClustersPass->Scale(1./gDataAntiProtonsTPCClustersPass->Integral(1,gDataAntiProtonsTPCClustersPass->GetNbinsX()));
  gDataAntiProtonsTPCClustersPass->SetMarkerStyle(1);
  TH1F *gDataAntiProtonsChi2PerClusterTPCPass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(3));
  gDataAntiProtonsChi2PerClusterTPCPass->SetStats(kFALSE);
  gDataAntiProtonsChi2PerClusterTPCPass->Add(gDataAntiProtonsChi2PerClusterTPCReject);
  gDataAntiProtonsChi2PerClusterTPCPass->Sumw2();
  gDataAntiProtonsChi2PerClusterTPCPass->Scale(1./gDataAntiProtonsChi2PerClusterTPCPass->Integral(1,gDataAntiProtonsChi2PerClusterTPCPass->GetNbinsX()));
  gDataAntiProtonsChi2PerClusterTPCPass->SetMarkerStyle(1);
  TH1F *gDataAntiProtonsExtCov11Pass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(4));
  gDataAntiProtonsExtCov11Pass->SetStats(kFALSE);
  gDataAntiProtonsExtCov11Pass->Add(gDataAntiProtonsExtCov11Reject);
  gDataAntiProtonsExtCov11Pass->Sumw2();
  gDataAntiProtonsExtCov11Pass->Scale(1./gDataAntiProtonsExtCov11Pass->Integral(1,gDataAntiProtonsExtCov11Pass->GetNbinsX()));
  gDataAntiProtonsExtCov11Pass->SetMarkerStyle(1);
  TH1F *gDataAntiProtonsExtCov22Pass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(5));
  gDataAntiProtonsExtCov22Pass->SetStats(kFALSE);
  gDataAntiProtonsExtCov22Pass->Add(gDataAntiProtonsExtCov22Reject);
  gDataAntiProtonsExtCov22Pass->Sumw2();
  gDataAntiProtonsExtCov22Pass->Scale(1./gDataAntiProtonsExtCov22Pass->Integral(1,gDataAntiProtonsExtCov22Pass->GetNbinsX()));
  gDataAntiProtonsExtCov22Pass->SetMarkerStyle(1);
  TH1F *gDataAntiProtonsExtCov33Pass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(6));
  gDataAntiProtonsExtCov33Pass->SetStats(kFALSE);
  gDataAntiProtonsExtCov33Pass->Add(gDataAntiProtonsExtCov33Reject);
  gDataAntiProtonsExtCov33Pass->Sumw2();
  gDataAntiProtonsExtCov33Pass->Scale(1./gDataAntiProtonsExtCov33Pass->Integral(1,gDataAntiProtonsExtCov33Pass->GetNbinsX()));
  gDataAntiProtonsExtCov33Pass->SetMarkerStyle(1);
  TH1F *gDataAntiProtonsExtCov44Pass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(7));
  gDataAntiProtonsExtCov44Pass->SetStats(kFALSE);
  gDataAntiProtonsExtCov44Pass->Add(gDataAntiProtonsExtCov44Reject);
  gDataAntiProtonsExtCov44Pass->Sumw2();
  gDataAntiProtonsExtCov44Pass->Scale(1./gDataAntiProtonsExtCov44Pass->Integral(1,gDataAntiProtonsExtCov44Pass->GetNbinsX()));
  gDataAntiProtonsExtCov44Pass->SetMarkerStyle(1);
  TH1F *gDataAntiProtonsExtCov55Pass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(8));
  gDataAntiProtonsExtCov55Pass->SetStats(kFALSE);
  gDataAntiProtonsExtCov55Pass->Add(gDataAntiProtonsExtCov55Reject);
  gDataAntiProtonsExtCov55Pass->Sumw2();
  gDataAntiProtonsExtCov55Pass->Scale(1./gDataAntiProtonsExtCov55Pass->Integral(1,gDataAntiProtonsExtCov55Pass->GetNbinsX()));
  gDataAntiProtonsExtCov55Pass->SetMarkerStyle(1);
  TH1F *gDataAntiProtonsDCAXYPass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(11));
  gDataAntiProtonsDCAXYPass->SetStats(kFALSE);  
  gDataAntiProtonsDCAXYPass->Add(gDataAntiProtonsDCAXYReject);
  gDataAntiProtonsDCAXYPass->Sumw2();
  gDataAntiProtonsDCAXYPass->Scale(1./gDataAntiProtonsDCAXYPass->Integral(1,gDataAntiProtonsDCAXYPass->GetNbinsX()));
  gDataAntiProtonsDCAXYPass->SetMarkerStyle(1);
  TH1F *gDataAntiProtonsDCAZPass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(13));
  gDataAntiProtonsDCAZPass->SetStats(kFALSE);  
  gDataAntiProtonsDCAZPass->Add(gDataAntiProtonsDCAZReject);
  gDataAntiProtonsDCAZPass->Sumw2();
  gDataAntiProtonsDCAZPass->Scale(1./gDataAntiProtonsDCAZPass->Integral(1,gDataAntiProtonsDCAZPass->GetNbinsX()));
  gDataAntiProtonsDCAZPass->SetMarkerStyle(1);
  TH1F *gDataAntiProtonsNumberOfTPCdEdxPointsPass = dynamic_cast<TH1F *>(fQADataAntiProtonsAcceptedList->At(26));
  gDataAntiProtonsNumberOfTPCdEdxPointsPass->SetStats(kFALSE);  
  gDataAntiProtonsNumberOfTPCdEdxPointsPass->Add(gDataAntiProtonsNumberOfTPCdEdxPointsReject);
  gDataAntiProtonsNumberOfTPCdEdxPointsPass->Sumw2();
  gDataAntiProtonsNumberOfTPCdEdxPointsPass->Scale(1./gDataAntiProtonsNumberOfTPCdEdxPointsPass->Integral(1,gDataAntiProtonsNumberOfTPCdEdxPointsPass->GetNbinsX()));
  gDataAntiProtonsNumberOfTPCdEdxPointsPass->SetMarkerStyle(1);


  //__________________________________________________//
  TCanvas *c1 = new TCanvas("c1","ITS clusters",0,0,600,400);
  c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->Divide(2,1);
  c1->cd(1); gDataProtonsITSClustersPass->Divide(gMCProtonsITSClustersPass);
  gDataProtonsITSClustersPass->SetTitle("Protons");
  gDataProtonsITSClustersPass->GetYaxis()->SetTitle("Data/MC");
  gDataProtonsITSClustersPass->Draw("E"); 
  c1->cd(2); 
  gDataAntiProtonsITSClustersPass->Divide(gMCAntiProtonsITSClustersPass);
  gDataAntiProtonsITSClustersPass->SetTitle("Antiprotons");
  gDataAntiProtonsITSClustersPass->GetYaxis()->SetTitle("Data/MC");
  gDataAntiProtonsITSClustersPass->Draw("E"); 
  c1->SaveAs("NClustersITS.gif");

  TCanvas *c3 = new TCanvas("c3","TPC clusters",0,200,600,400);
  c3->SetFillColor(10); c3->SetHighLightColor(10);
  c3->Divide(2,1);
  c3->cd(1); gDataProtonsTPCClustersPass->Divide(gMCProtonsTPCClustersPass);
  gDataProtonsTPCClustersPass->SetTitle("Protons");
  gDataProtonsTPCClustersPass->GetYaxis()->SetTitle("Data/MC");
  gDataProtonsTPCClustersPass->Draw("E");
  c3->cd(2); 
  gDataAntiProtonsTPCClustersPass->Divide(gMCAntiProtonsTPCClustersPass);
  gDataAntiProtonsTPCClustersPass->SetTitle("Antirotons");
  gDataAntiProtonsTPCClustersPass->GetYaxis()->SetTitle("Data/MC");
  gDataAntiProtonsTPCClustersPass->Draw("E");
  c3->SaveAs("NClustersTPC.gif");

  TCanvas *c4 = new TCanvas("c4","chi^2 per TPC cluster",0,300,600,400);
  c4->SetFillColor(10); c4->SetHighLightColor(10);
  c4->Divide(2,1);
  c4->cd(1); 
  gDataProtonsChi2PerClusterTPCPass->Divide(gMCProtonsChi2PerClusterTPCPass);
  gDataProtonsChi2PerClusterTPCPass->SetTitle("Protons");
  gDataProtonsChi2PerClusterTPCPass->GetYaxis()->SetTitle("Data/MC");
  gDataProtonsChi2PerClusterTPCPass->Draw("E"); 
  c4->cd(2); 
  gDataAntiProtonsChi2PerClusterTPCPass->Divide(gMCAntiProtonsChi2PerClusterTPCPass);
  gDataAntiProtonsChi2PerClusterTPCPass->SetTitle("Antirotons");
  gDataAntiProtonsChi2PerClusterTPCPass->GetYaxis()->SetTitle("Data/MC");
  gDataAntiProtonsChi2PerClusterTPCPass->Draw("E"); 
  c4->SaveAs("Chi2PerTPCCluster.gif");

  if(gMCProtonsExtCov11Pass->GetEntries() != 0) {
    TCanvas *c5 = new TCanvas("c5","Cov11",0,400,600,400);
    c5->SetFillColor(10); c5->SetHighLightColor(10);
    c5->Divide(2,1);
    c5->cd(1); 
    gDataProtonsExtCov11Pass->Divide(gMCProtonsExtCov11Pass);
    gDataProtonsExtCov11Pass->SetTitle("Protons");
    gDataProtonsExtCov11Pass->GetYaxis()->SetTitle("Data/MC");
    gDataProtonsExtCov11Pass->Draw("E"); 
    c5->cd(2); 
    gDataAntiProtonsExtCov11Pass->Divide(gMCAntiProtonsExtCov11Pass);
    gDataAntiProtonsExtCov11Pass->SetTitle("Antiprotons");
    gDataAntiProtonsExtCov11Pass->GetYaxis()->SetTitle("Data/MC");
    gDataAntiProtonsExtCov11Pass->Draw("E"); 
    c5->SaveAs("cov11.gif");
  }

  if(gMCProtonsExtCov22Pass->GetEntries() != 0) {
    TCanvas *c6 = new TCanvas("c6","Cov22",0,500,600,400);
    c6->SetFillColor(10); c6->SetHighLightColor(10);
    c6->Divide(2,1);
    c6->cd(1); 
    gDataProtonsExtCov22Pass->Divide(gMCProtonsExtCov11Pass);
    gDataProtonsExtCov22Pass->SetTitle("Protons");
    gDataProtonsExtCov22Pass->GetYaxis()->SetTitle("Data/MC");
    gDataProtonsExtCov22Pass->Draw("E"); 
    c6->cd(2); 
    gDataAntiProtonsExtCov22Pass->Divide(gMCAntiProtonsExtCov22Pass);
    gDataAntiProtonsExtCov22Pass->SetTitle("Antiprotons");
    gDataAntiProtonsExtCov22Pass->GetYaxis()->SetTitle("Data/MC");
    gDataAntiProtonsExtCov22Pass->Draw("E"); 
    c6->SaveAs("cov22.gif");
  }

  if(gMCProtonsExtCov33Pass->GetEntries() != 0) {
    TCanvas *c7 = new TCanvas("c7","Cov33",600,0,600,400);
    c7->SetFillColor(10); c7->SetHighLightColor(10);
    c7->Divide(2,1);
    c7->cd(1); 
    gDataProtonsExtCov33Pass->Divide(gMCProtonsExtCov33Pass);
    gDataProtonsExtCov33Pass->SetTitle("Protons");
    gDataProtonsExtCov33Pass->GetYaxis()->SetTitle("Data/MC");
    gDataProtonsExtCov33Pass->Draw("E"); 
    c7->cd(2); 
    gDataAntiProtonsExtCov33Pass->Divide(gMCAntiProtonsExtCov33Pass);
    gDataAntiProtonsExtCov33Pass->SetTitle("Antiprotons");
    gDataAntiProtonsExtCov33Pass->GetYaxis()->SetTitle("Data/MC");
    gDataAntiProtonsExtCov33Pass->Draw("E"); 
    c7->SaveAs("cov33.gif");
  }

  if(gMCProtonsExtCov44Pass->GetEntries() != 0) {
    TCanvas *c8 = new TCanvas("c8","Cov44",600,100,600,400);
    c8->SetFillColor(10); c8->SetHighLightColor(10);
    c8->Divide(2,1);
    c8->cd(1); 
    gDataProtonsExtCov44Pass->Divide(gMCProtonsExtCov44Pass);
    gDataProtonsExtCov44Pass->SetTitle("Protons");
    gDataProtonsExtCov44Pass->GetYaxis()->SetTitle("Data/MC");
    gDataProtonsExtCov44Pass->Draw("E"); 
    c8->cd(2); 
    gDataAntiProtonsExtCov44Pass->Divide(gMCAntiProtonsExtCov44Pass);
    gDataAntiProtonsExtCov44Pass->SetTitle("Antiprotons");
    gDataAntiProtonsExtCov44Pass->GetYaxis()->SetTitle("Data/MC");
    gDataAntiProtonsExtCov44Pass->Draw("E"); 
    c8->SaveAs("cov44.gif");
  }

  if(gMCProtonsExtCov55Pass->GetEntries() != 0) {
    TCanvas *c9 = new TCanvas("c9","Cov55",600,200,600,400);
    c9->SetFillColor(10); c9->SetHighLightColor(10);
    c9->Divide(2,1);
    c9->cd(1); 
    gDataProtonsExtCov55Pass->Divide(gMCProtonsExtCov55Pass);
    gDataProtonsExtCov55Pass->SetTitle("Protons");
    gDataProtonsExtCov55Pass->GetYaxis()->SetTitle("Data/MC");
    gDataProtonsExtCov55Pass->Draw("E"); 
    c9->cd(2); 
    gDataAntiProtonsExtCov55Pass->Divide(gMCAntiProtonsExtCov55Pass);
    gDataAntiProtonsExtCov55Pass->SetTitle("Antiprotons");
    gDataAntiProtonsExtCov55Pass->GetYaxis()->SetTitle("Data/MC");
    gDataAntiProtonsExtCov55Pass->Draw("E"); 
    c9->SaveAs("cov55.gif");
  }

  if(gMCProtonsDCAXYPass->GetEntries() != 0) {
    TCanvas *c12 = new TCanvas("c12","dca(xy)",600,500,600,400);
    c12->SetFillColor(10); c12->SetHighLightColor(10);
    c12->Divide(2,1);
    c12->cd(1); 
    gDataProtonsDCAXYPass->Divide(gMCProtonsDCAXYPass);
    gDataProtonsDCAXYPass->SetTitle("Protons");
    gDataProtonsDCAXYPass->GetYaxis()->SetTitle("Data/MC");
    gDataProtonsDCAXYPass->Draw("E"); 
    c12->cd(2); 
    gDataAntiProtonsDCAXYPass->Divide(gMCAntiProtonsDCAXYPass);
    gDataAntiProtonsDCAXYPass->SetTitle("Antiprotons");
    gDataAntiProtonsDCAXYPass->GetYaxis()->SetTitle("Data/MC");
    gDataAntiProtonsDCAXYPass->Draw("E"); 
    c12->SaveAs("dcaXY.gif");
  }

  if(gMCProtonsDCAZPass->GetEntries() != 0) {
    TCanvas *c14 = new TCanvas("c14","dca(z)",1200,100,600,400);
    c14->SetFillColor(10); c14->SetHighLightColor(10);
    c14->Divide(2,1);
    c14->cd(1); 
    gDataProtonsDCAZPass->Divide(gMCProtonsDCAZPass);
    gDataProtonsDCAZPass->SetTitle("Protons");
    gDataProtonsDCAZPass->GetYaxis()->SetTitle("Data/MC");
    gDataProtonsDCAZPass->Draw("E"); 
    c14->cd(2); 
    gDataAntiProtonsDCAZPass->Divide(gMCAntiProtonsDCAZPass);
    gDataAntiProtonsDCAZPass->SetTitle("Antiprotons");
    gDataAntiProtonsDCAZPass->GetYaxis()->SetTitle("Data/MC");
    gDataAntiProtonsDCAZPass->Draw("E"); 
    c14->SaveAs("dcaZ.gif");
  }

  TCanvas *c16 = new TCanvas("c16","TPC clusters (dE/dx)",1200,300,600,400);
  c16->SetFillColor(10); c16->SetHighLightColor(10);
  c16->Divide(2,1);
  c16->cd(1); 
  gDataProtonsNumberOfTPCdEdxPointsPass->Divide(gMCProtonsNumberOfTPCdEdxPointsPass);
  gDataProtonsNumberOfTPCdEdxPointsPass->SetTitle("Protons");
  gDataProtonsNumberOfTPCdEdxPointsPass->GetYaxis()->SetTitle("Data/MC");
  gDataProtonsNumberOfTPCdEdxPointsPass->Draw("E"); 
  c16->cd(2); 
  gDataAntiProtonsNumberOfTPCdEdxPointsPass->Divide(gMCAntiProtonsNumberOfTPCdEdxPointsPass);
  gDataAntiProtonsNumberOfTPCdEdxPointsPass->SetTitle("Antiprotons");
  gDataAntiProtonsNumberOfTPCdEdxPointsPass->GetYaxis()->SetTitle("Data/MC");
  gDataAntiProtonsNumberOfTPCdEdxPointsPass->Draw("E"); 
  c16->SaveAs("NClustersTPCdEdx.gif");
}

void drawDCAPlots(const char* fileName) {
  //Function to draw the DCA plots for protons and antiprotons
  TString histoTitle, newTitle;
  Int_t iCounter = 0;
  Double_t etaMin = 0.0, etaMax = 0.0, etaStep = 0.0;
  Double_t ptMin = 0.45, ptMax = 0.0, ptStep = 0.1;

  TFile *f = TFile::Open(fileName);
  TList *listQA = dynamic_cast<TList *>(f->Get("outputQAList"));
  TList *gListGlobalQA = dynamic_cast<TList *>(listQA->At(0));
  TList *fQA2DList = dynamic_cast<TList *>(gListGlobalQA->At(0));

  //2D dca vs pT - accepted protons & antiprotons
  TH3F *gHistDCAxyEtaPtProtons = dynamic_cast<TH3F *>(fQA2DList->At(14));
  gHistDCAxyEtaPtProtons->SetStats(kFALSE);
  TH3F *gHistDCAzEtaPtProtons = dynamic_cast<TH3F *>(fQA2DList->At(15));
  gHistDCAzEtaPtProtons->SetStats(kFALSE);
  TH3F *gHistDCAxyEtaPtAntiProtons = dynamic_cast<TH3F *>(fQA2DList->At(16));
  gHistDCAxyEtaPtAntiProtons->SetStats(kFALSE);
  TH3F *gHistDCAzEtaPtAntiProtons = dynamic_cast<TH3F *>(fQA2DList->At(17));
  gHistDCAzEtaPtAntiProtons->SetStats(kFALSE);

  //Protons dcaXY
  ptMin = 0.45, ptMax = 0.0, ptStep = 0.1;
  TH1D *gHistProtonsDCAxy[200];
  //ptStep = (gHistDCAxyEtaPtProtons->GetYaxis()->GetXmax() - 
  //gHistDCAxyEtaPtProtons->GetYaxis()->GetXmin())/gHistDCAxyEtaPtProtons->GetNbinsY();
  for(Int_t iBinY = 1; iBinY <= gHistDCAxyEtaPtProtons->GetNbinsY(); iBinY++) {
    //ptMin = gHistDCAxyEtaPtProtons->GetYaxis()->GetBinCenter(iBinY) - ptStep/2.;
    ptMax = ptMin + ptStep;
    //Printf("Pt: %lf - %lf",ptMin,ptMax);
    histoTitle = "gHistProtonsDCAxy_etaBin"; //histoTitle += iBinX;
    histoTitle += "_ptBin"; histoTitle += iBinY;
    newTitle = "(Protons) "; newTitle += ptMin; newTitle += " < Pt < ";
    newTitle += ptMax; newTitle += "GeV/c";
    ptMin += ptStep;

    //gHistPrimaryProtonsDCAxy[iCounter] = gDCAListHistograms3D[0]->ProjectionZ(histoTitle.Data(),iBinX,iBinX,iBinY,iBinY,"e");
    gHistProtonsDCAxy[iCounter] = gHistDCAxyEtaPtProtons->ProjectionZ(histoTitle.Data(),1,gHistDCAxyEtaPtProtons->GetNbinsX(),iBinY,iBinY,"e");
    gHistProtonsDCAxy[iCounter]->SetStats(kFALSE);
    gHistProtonsDCAxy[iCounter]->SetMarkerStyle(24);
    gHistProtonsDCAxy[iCounter]->SetTitle(newTitle.Data());
    gHistProtonsDCAxy[iCounter]->GetXaxis()->SetTitle("dca_{(xy)} [cm]");
    gHistProtonsDCAxy[iCounter]->GetYaxis()->SetTitle("Entries");
    iCounter += 1;
  }//loop over y axis

  //Antiprotons dcaXY
  ptMin = 0.45, ptMax = 0.0, ptStep = 0.1;
  iCounter = 0;
  TH1D *gHistAntiProtonsDCAxy[200];
  //ptStep = (gHistDCAxyEtaPtAntiProtons->GetYaxis()->GetXmax() - 
  //gHistDCAxyEtaPtAntiProtons->GetYaxis()->GetXmin())/gHistDCAxyEtaPtAntiProtons->GetNbinsY();
  for(Int_t iBinY = 1; iBinY <= gHistDCAxyEtaPtAntiProtons->GetNbinsY(); iBinY++) {
    //ptMin = gHistDCAxyEtaPtAntiProtons->GetYaxis()->GetBinCenter(iBinY) - ptStep/2.;
    ptMax = ptMin + ptStep;
    //Printf("Pt: %lf - %lf",ptMin,ptMax);
    histoTitle = "gHistAntiProtonsDCAxy_etaBin"; //histoTitle += iBinX;
    histoTitle += "_ptBin"; histoTitle += iBinY;
    newTitle = "(AntiProtons) "; newTitle += ptMin; newTitle += " < Pt < ";
    newTitle += ptMax; newTitle += "GeV/c";
    ptMin += ptStep;

    //gHistPrimaryAntiProtonsDCAxy[iCounter] = gDCAListHistograms3D[0]->ProjectionZ(histoTitle.Data(),iBinX,iBinX,iBinY,iBinY,"e");
    gHistAntiProtonsDCAxy[iCounter] = gHistDCAxyEtaPtAntiProtons->ProjectionZ(histoTitle.Data(),1,gHistDCAxyEtaPtAntiProtons->GetNbinsX(),iBinY,iBinY,"e");
    gHistAntiProtonsDCAxy[iCounter]->SetStats(kFALSE);
    gHistAntiProtonsDCAxy[iCounter]->SetMarkerStyle(20);
    gHistAntiProtonsDCAxy[iCounter]->SetTitle(newTitle.Data());
    gHistAntiProtonsDCAxy[iCounter]->GetXaxis()->SetTitle("dca_{(xy)} [cm]");
    gHistAntiProtonsDCAxy[iCounter]->GetYaxis()->SetTitle("Entries");
    iCounter += 1;
  }//loop over y axis


  //Protons dcaZ
  ptMin = 0.45, ptMax = 0.0, ptStep = 0.1;
  iCounter = 0;
  TH1D *gHistProtonsDCAz[200];
  //ptStep = (gHistDCAzEtaPtProtons->GetYaxis()->GetXmax() - 
  //gHistDCAzEtaPtProtons->GetYaxis()->GetXmin())/gHistDCAzEtaPtProtons->GetNbinsY();
  for(Int_t iBinY = 1; iBinY <= gHistDCAzEtaPtProtons->GetNbinsY(); iBinY++) {
    //ptMin = gHistDCAzEtaPtProtons->GetYaxis()->GetBinCenter(iBinY) - ptStep/2.;
    ptMax = ptMin + ptStep;
    //Printf("Pt: %lf - %lf",ptMin,ptMax);
    histoTitle = "gHistProtonsDCAz_etaBin"; //histoTitle += iBinX;
    histoTitle += "_ptBin"; histoTitle += iBinY;
    newTitle = "(Protons) "; newTitle += ptMin; newTitle += " < Pt < ";
    newTitle += ptMax; newTitle += "GeV/c";
    ptMin += ptStep;

    //gHistPrimaryProtonsDCAz[iCounter] = gDCAListHistograms3D[0]->ProjectionZ(histoTitle.Data(),iBinX,iBinX,iBinY,iBinY,"e");
    gHistProtonsDCAz[iCounter] = gHistDCAzEtaPtProtons->ProjectionZ(histoTitle.Data(),1,gHistDCAzEtaPtProtons->GetNbinsX(),iBinY,iBinY,"e");
    gHistProtonsDCAz[iCounter]->SetStats(kFALSE);
    gHistProtonsDCAz[iCounter]->SetMarkerStyle(24);
    gHistProtonsDCAz[iCounter]->SetTitle(newTitle.Data());
    gHistProtonsDCAz[iCounter]->GetXaxis()->SetTitle("dca_{(z)} [cm]");
    gHistProtonsDCAz[iCounter]->GetYaxis()->SetTitle("Entries");
    iCounter += 1;
  }//loop over y axis

  //Antiprotons dcaZ
  ptMin = 0.45, ptMax = 0.0, ptStep = 0.1;
  iCounter = 0;
  TH1D *gHistAntiProtonsDCAz[200];
  //ptStep = (gHistDCAzEtaPtAntiProtons->GetYaxis()->GetXmax() - 
  //gHistDCAzEtaPtAntiProtons->GetYaxis()->GetXmin())/gHistDCAzEtaPtAntiProtons->GetNbinsY();
  for(Int_t iBinY = 1; iBinY <= gHistDCAzEtaPtAntiProtons->GetNbinsY(); iBinY++) {
    //ptMin = gHistDCAzEtaPtAntiProtons->GetYaxis()->GetBinCenter(iBinY) - ptStep/2.;
    ptMax = ptMin + ptStep;
    //Printf("Pt: %lf - %lf",ptMin,ptMax);
    histoTitle = "gHistAntiProtonsDCAz_etaBin"; //histoTitle += iBinX;
    histoTitle += "_ptBin"; histoTitle += iBinY;
    newTitle = "(AntiProtons) "; newTitle += ptMin; newTitle += " < Pt < ";
    newTitle += ptMax; newTitle += "GeV/c";
    ptMin += ptStep;

    //gHistPrimaryAntiProtonsDCAz[iCounter] = gDCAListHistograms3D[0]->ProjectionZ(histoTitle.Data(),iBinX,iBinX,iBinY,iBinY,"e");
    gHistAntiProtonsDCAz[iCounter] = gHistDCAzEtaPtAntiProtons->ProjectionZ(histoTitle.Data(),1,gHistDCAzEtaPtAntiProtons->GetNbinsX(),iBinY,iBinY,"e");
    gHistAntiProtonsDCAz[iCounter]->SetStats(kFALSE);
    gHistAntiProtonsDCAz[iCounter]->SetMarkerStyle(20);
    gHistAntiProtonsDCAz[iCounter]->SetTitle(newTitle.Data());
    gHistAntiProtonsDCAz[iCounter]->GetXaxis()->SetTitle("dca_{(z)} [cm]");
    gHistAntiProtonsDCAz[iCounter]->GetYaxis()->SetTitle("Entries");
    iCounter += 1;
  }//loop over y axis



  //==========================================================//
  TCanvas *c1 = new TCanvas("c1","DCA(3D)",0,0,600,600);
  c1->SetFillColor(10); c1->SetHighLightColor(10); c1->Divide(2,2);
  c1->cd(1); gHistDCAxyEtaPtProtons->Draw("box");
  c1->cd(2); gHistDCAzEtaPtProtons->Draw("box");
  c1->cd(3); gHistDCAxyEtaPtAntiProtons->Draw("box");
  c1->cd(4); gHistDCAzEtaPtAntiProtons->Draw("box");

  //==========================================================//
  TCanvas *c2 = new TCanvas("c2","DCA(xy)",100,100,900,600);
  c2->SetFillColor(10); c2->SetHighLightColor(10); c2->Divide(3,2);
  c2->cd(1)->SetLogy(); gHistProtonsDCAxy[0]->Draw("E");
  gHistAntiProtonsDCAxy[0]->Draw("ESAME");
  c2->cd(2)->SetLogy(); gHistProtonsDCAxy[1]->Draw("E");
  gHistAntiProtonsDCAxy[1]->Draw("ESAME");
  c2->cd(3)->SetLogy(); gHistProtonsDCAxy[2]->Draw("E");
  gHistAntiProtonsDCAxy[2]->Draw("ESAME");
  c2->cd(4)->SetLogy(); gHistProtonsDCAxy[3]->Draw("E");
  gHistAntiProtonsDCAxy[3]->Draw("ESAME");
  c2->cd(5)->SetLogy(); gHistProtonsDCAxy[4]->Draw("E");
  gHistAntiProtonsDCAxy[4]->Draw("ESAME");
  c2->cd(6)->SetLogy(); gHistProtonsDCAxy[5]->Draw("E");
  gHistAntiProtonsDCAxy[5]->Draw("ESAME");

  TFile *fout = TFile::Open("test.root","recreate");
  gHistProtonsDCAxy[0]->Write();
  gHistProtonsDCAxy[1]->Write();
  gHistProtonsDCAxy[2]->Write();
  gHistProtonsDCAxy[3]->Write();
  gHistProtonsDCAxy[4]->Write();
  gHistProtonsDCAxy[5]->Write();
  gHistAntiProtonsDCAxy[0]->Write();
  gHistAntiProtonsDCAxy[1]->Write();
  gHistAntiProtonsDCAxy[2]->Write();
  gHistAntiProtonsDCAxy[3]->Write();
  gHistAntiProtonsDCAxy[4]->Write();
  gHistAntiProtonsDCAxy[5]->Write();
  gHistProtonsDCAz[0]->Write();
  gHistProtonsDCAz[1]->Write();
  gHistProtonsDCAz[2]->Write();
  gHistProtonsDCAz[3]->Write();
  gHistProtonsDCAz[4]->Write();
  gHistProtonsDCAz[5]->Write();
  gHistAntiProtonsDCAz[0]->Write();
  gHistAntiProtonsDCAz[1]->Write();
  gHistAntiProtonsDCAz[2]->Write();
  gHistAntiProtonsDCAz[3]->Write();
  gHistAntiProtonsDCAz[4]->Write();
  gHistAntiProtonsDCAz[5]->Write();
  fout->Close();
  //==========================================================//
  TCanvas *c3 = new TCanvas("c3","DCA(z)",200,200,900,600);
  c3->SetFillColor(10); c3->SetHighLightColor(10); c3->Divide(3,2);
  c3->cd(1)->SetLogy(); gHistProtonsDCAz[0]->Draw("E");
  gHistAntiProtonsDCAz[0]->Draw("ESAME");
  c3->cd(2)->SetLogy(); gHistProtonsDCAz[1]->Draw("E");
  gHistAntiProtonsDCAz[1]->Draw("ESAME");
  c3->cd(3)->SetLogy(); gHistProtonsDCAz[2]->Draw("E");
  gHistAntiProtonsDCAz[2]->Draw("ESAME");
  c3->cd(4)->SetLogy(); gHistProtonsDCAz[3]->Draw("E");
  gHistAntiProtonsDCAz[3]->Draw("ESAME");
  c3->cd(5)->SetLogy(); gHistProtonsDCAz[4]->Draw("E");
  gHistAntiProtonsDCAz[4]->Draw("ESAME");
  c3->cd(6)->SetLogy(); gHistProtonsDCAz[5]->Draw("E");
  gHistAntiProtonsDCAz[5]->Draw("ESAME");

}

