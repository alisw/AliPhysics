///////////////////////////////////////////////////////////
// PlotAndSave.C (called by AODQAChecks.C)               //
//                                                       //
// Written by John Groh                                  //
///////////////////////////////////////////////////////////

void PlotAndSave(Int_t runs[],
		 const Int_t nRuns,
		 Bool_t useMC,
		 TFile*& fout,

		 TH1F*& TPCnsigMeanTrendPion,
		 TH1F*& TPCnsigMeanTrendKaon,
		 TH1F*& TPCnsigMeanTrendProton,
		 TH1F*& TPCnsigSigmaTrendPion,
		 TH1F*& TPCnsigSigmaTrendKaon,
		 TH1F*& TPCnsigSigmaTrendProton,
		 TH1F*& TOFnsigMeanTrendPion,
		 TH1F*& TOFnsigMeanTrendKaon,
		 TH1F*& TOFnsigMeanTrendProton,
		 TH1F*& TOFnsigSigmaTrendPion,
		 TH1F*& TOFnsigSigmaTrendKaon,
		 TH1F*& TOFnsigSigmaTrendProton,

		 TH1F*& IntegRawYieldAll,
		 TH1F*& IntegRawYieldPiPlus,
		 TH1F*& IntegRawYieldKPlus,
		 TH1F*& IntegRawYieldProton,
		 TH1F*& IntegRawYieldPiMinus,
		 TH1F*& IntegRawYieldKMinus,
		 TH1F*& IntegRawYieldAntiproton,

		 TH1F*& EfficiencyPiPlus,
		 TH1F*& EfficiencyKPlus,
		 TH1F*& EfficiencyProton,
		 TH1F*& EfficiencyPiMinus,
		 TH1F*& EfficiencyKMinus,
		 TH1F*& EfficiencyAntiproton,

		 TH1F*& MatchEffPos,
		 TH1F*& MatchEffNeg)
{
  Printf("\n\n\n--- Entering funtion PlotAndSave() ---\n\n\n");

  // this gets rid of scientific notation for run numbers on the axes below
  TGaxis::SetMaxDigits(7);

  //-----------------------------------------------------------------------------------------
  // Draw the trends in the means and sigmas of the fits to the peaks of the projections at -
  // fixed Pt as a function of the run number                                               -
  //-----------------------------------------------------------------------------------------
  TCanvas * cNSigmaStability = new TCanvas("cNSigmaStability","cNSigmaStability",100,50,700,500);
  cNSigmaStability->Divide(2,2);
  // TPC means
  cNSigmaStability->cd(1);
  TH2F * hAxesTPCnsigMeans = new TH2F("hAxesTPCnsigMeans","",nRuns,0,nRuns,100,-2,2);
  hAxesTPCnsigMeans->GetYaxis()->SetTitle("TPC nsigma particle peak means");
  hAxesTPCnsigMeans->SetStats(kFALSE);
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesTPCnsigMeans->GetXaxis()->SetBinLabel(irun+1, Form("%i",runs[irun]));
  hAxesTPCnsigMeans->DrawCopy();
  TLegend * lTPCMeans = new TLegend(0.7,0.15,0.85,0.35);
  TPCnsigMeanTrendPion->SetMarkerStyle(Marker[0]);
  TPCnsigMeanTrendPion->SetMarkerColor(Color[0]);
  TPCnsigMeanTrendPion->SetLineColor(Color[0]);
  TPCnsigMeanTrendPion->SetStats(kFALSE);
  TPCnsigMeanTrendPion->DrawCopy("E1X0Psame");
  lTPCMeans->AddEntry(TPCnsigMeanTrendPion,"#pi^{+}, #pi^{-}","lpe");
  TPCnsigMeanTrendKaon->SetMarkerStyle(Marker[1]);
  TPCnsigMeanTrendKaon->SetMarkerColor(Color[1]);
  TPCnsigMeanTrendKaon->SetLineColor(Color[1]);
  TPCnsigMeanTrendKaon->SetStats(kFALSE);
  TPCnsigMeanTrendKaon->DrawCopy("E1X0Psame");
  lTPCMeans->AddEntry(TPCnsigMeanTrendKaon,"K^{+}, K^{-}","lpe");
  TPCnsigMeanTrendProton->SetMarkerStyle(Marker[2]);
  TPCnsigMeanTrendProton->SetMarkerColor(Color[2]);
  TPCnsigMeanTrendProton->SetLineColor(Color[2]);  
  TPCnsigMeanTrendProton->SetStats(kFALSE);
  TPCnsigMeanTrendProton->DrawCopy("E1X0Psame");
  lTPCMeans->AddEntry(TPCnsigMeanTrendProton,"p, #bar{p}","lpe");
  lTPCMeans->SetFillColor(0);
  lTPCMeans->DrawClone();
  // TPC sigmas
  cNSigmaStability->cd(3);
  TH2F * hAxesTPCnsigSigmas = new TH2F("hAxesTPCnsigSigmas","",nRuns,0,nRuns,100,-1,3);
  hAxesTPCnsigSigmas->GetYaxis()->SetTitle("TPC nsigma particle peak sigmas");
  hAxesTPCnsigSigmas->SetStats(kFALSE);
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesTPCnsigSigmas->GetXaxis()->SetBinLabel(irun+1, Form("%i",runs[irun]));
  hAxesTPCnsigSigmas->DrawCopy();
  TPCnsigSigmaTrendPion->SetMarkerStyle(Marker[0]);
  TPCnsigSigmaTrendPion->SetMarkerColor(Color[0]);
  TPCnsigSigmaTrendPion->SetLineColor(Color[0]);
  TPCnsigSigmaTrendPion->SetStats(kFALSE);
  TPCnsigSigmaTrendPion->DrawCopy("X0E1Psame");
  TPCnsigSigmaTrendKaon->SetMarkerStyle(Marker[1]);
  TPCnsigSigmaTrendKaon->SetMarkerColor(Color[1]);
  TPCnsigSigmaTrendKaon->SetLineColor(Color[1]);
  TPCnsigSigmaTrendKaon->SetStats(kFALSE);
  TPCnsigSigmaTrendKaon->DrawCopy("X0E1Psame");
  TPCnsigSigmaTrendProton->SetMarkerStyle(Marker[2]);
  TPCnsigSigmaTrendProton->SetMarkerColor(Color[2]);
  TPCnsigSigmaTrendProton->SetLineColor(Color[2]);
  TPCnsigSigmaTrendProton->SetStats(kFALSE);
  TPCnsigSigmaTrendProton->DrawCopy("X0E1Psame");
  TLegend * lTPCSigmas = new TLegend(0.7,0.15,0.85,0.35);
  lTPCSigmas->AddEntry(TPCnsigSigmaTrendPion,"#pi^{+}, #pi^{-}","lpe");
  lTPCSigmas->AddEntry(TPCnsigSigmaTrendKaon,"K^{+}, K^{-}","lpe");
  lTPCSigmas->AddEntry(TPCnsigSigmaTrendProton,"p, #bar{p}","lpe");
  lTPCSigmas->SetFillColor(0);
  lTPCSigmas->DrawClone();
  // TOF means
  cNSigmaStability->cd(2);
  TH2F * hAxesTOFnsigMeans = new TH2F("hAxesTOFnsigMeans","",nRuns,0,nRuns,100,-2,2);
  hAxesTOFnsigMeans->GetYaxis()->SetTitle("TOF nsigma particle peak means");
  hAxesTOFnsigMeans->SetStats(kFALSE);
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesTOFnsigMeans->GetXaxis()->SetBinLabel(irun+1, Form("%i",runs[irun]));
  hAxesTOFnsigMeans->DrawCopy();
  TOFnsigMeanTrendPion->SetMarkerStyle(Marker[0]);
  TOFnsigMeanTrendPion->SetMarkerColor(Color[0]);
  TOFnsigMeanTrendPion->SetLineColor(Color[0]);
  TOFnsigMeanTrendPion->SetStats(kFALSE);
  TOFnsigMeanTrendPion->DrawCopy("X0E1Psame");
  TOFnsigMeanTrendKaon->SetMarkerStyle(Marker[1]);
  TOFnsigMeanTrendKaon->SetMarkerColor(Color[1]);
  TOFnsigMeanTrendKaon->SetLineColor(Color[1]);
  TOFnsigMeanTrendKaon->SetStats(kFALSE);
  TOFnsigMeanTrendKaon->DrawCopy("X0E1Psame");
  TOFnsigMeanTrendProton->SetMarkerStyle(Marker[2]);
  TOFnsigMeanTrendProton->SetMarkerColor(Color[2]);
  TOFnsigMeanTrendProton->SetLineColor(Color[2]);  
  TOFnsigMeanTrendProton->SetStats(kFALSE);
  TOFnsigMeanTrendProton->DrawCopy("X0E1Psame");
  TLegend * lTOFMeans = new TLegend(0.7,0.15,0.85,0.35);
  lTOFMeans->AddEntry(TOFnsigMeanTrendPion,"#pi^{+}, #pi^{-}","lpe");
  lTOFMeans->AddEntry(TOFnsigMeanTrendKaon,"K^{+}, K^{-}","lpe");
  lTOFMeans->AddEntry(TOFnsigMeanTrendProton,"p, #bar{p}","lpe");
  lTOFMeans->SetFillColor(0);
  lTOFMeans->DrawClone();
  // TOF sigmas
  cNSigmaStability->cd(4);
  TH2F * hAxesTOFnsigSigmas = new TH2F("hAxesTOFnsigSigmas","",nRuns,0,nRuns,100,-1,2);
  hAxesTOFnsigSigmas->GetYaxis()->SetTitle("TOF nsigma particle peak sigmas");
  hAxesTOFnsigSigmas->SetStats(kFALSE);
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesTOFnsigSigmas->GetXaxis()->SetBinLabel(irun+1, Form("%i",runs[irun]));
  hAxesTOFnsigSigmas->DrawCopy();
  TOFnsigSigmaTrendPion->SetMarkerStyle(Marker[0]);
  TOFnsigSigmaTrendPion->SetMarkerColor(Color[0]);
  TOFnsigSigmaTrendPion->SetLineColor(Color[0]);
  TOFnsigSigmaTrendPion->SetStats(kFALSE);
  TOFnsigSigmaTrendPion->DrawCopy("X0E1Psame");
  TOFnsigSigmaTrendKaon->SetMarkerStyle(Marker[1]);
  TOFnsigSigmaTrendKaon->SetMarkerColor(Color[1]);
  TOFnsigSigmaTrendKaon->SetLineColor(Color[1]);
  TOFnsigSigmaTrendKaon->SetStats(kFALSE);
  TOFnsigSigmaTrendKaon->DrawCopy("X0E1Psame");
  TOFnsigSigmaTrendProton->SetMarkerStyle(Marker[2]);
  TOFnsigSigmaTrendProton->SetMarkerColor(Color[2]);
  TOFnsigSigmaTrendProton->SetLineColor(Color[2]);
  TOFnsigSigmaTrendProton->SetStats(kFALSE);
  TOFnsigSigmaTrendProton->DrawCopy("X0E1Psame");
  TLegend * lTOFSigmas = new TLegend(0.7,0.15,0.85,0.35);
  lTOFSigmas->AddEntry(TOFnsigSigmaTrendPion,"#pi^{+}, #pi^{-}","lpe");
  lTOFSigmas->AddEntry(TOFnsigSigmaTrendKaon,"K^{+}, K^{-}","lpe");
  lTOFSigmas->AddEntry(TOFnsigSigmaTrendProton,"p, #bar{p}","lpe");
  lTOFSigmas->SetFillColor(0);
  lTOFSigmas->DrawClone();
  // write the canvas to the file
  fout->cd();
  cNSigmaStability->Write();

  //------------------------------------------------------------------------------
  // Draw the trends in the integrated raw yield as a function of the run number -
  //------------------------------------------------------------------------------
  TCanvas * cIntegRawYieldStability = new TCanvas("cIntegRawYieldStability","cIntegRawYieldStability",150,75,700,500);
  gPad->SetLogy();
  TLegend * lIntegRawYield = new TLegend(0.75,0.4,0.85,0.6);
  lIntegRawYield->SetFillColor(0);
  // all particles
  for (Int_t irun=0; irun<nRuns; irun++)
    IntegRawYieldAll->GetXaxis()->SetBinLabel(irun+1, Form("%i",runs[irun]));
  IntegRawYieldAll->SetTitle("Raw Yield, Integrated over all p_{T}");
  IntegRawYieldAll->SetStats(kFALSE);
  IntegRawYieldAll->SetMarkerStyle(34);
  IntegRawYieldAll->SetMarkerColor(kGreen);
  IntegRawYieldAll->SetLineColor(kGreen);
  IntegRawYieldAll->GetYaxis()->SetRangeUser(50,10000);
  IntegRawYieldAll->DrawCopy("E1X0P");
  lIntegRawYield->AddEntry(IntegRawYieldAll,"All Particles","p");
  // pi+
  IntegRawYieldPiPlus->SetMarkerStyle(Marker[0]);
  IntegRawYieldPiPlus->SetMarkerColor(Color[0]);
  IntegRawYieldPiPlus->SetLineColor(Color[0]);
  IntegRawYieldPiPlus->DrawCopy("E1X0Psame");
  lIntegRawYield->AddEntry(IntegRawYieldPiPlus,"#pi^{+}","p");
  // K+
  IntegRawYieldKPlus->SetMarkerStyle(Marker[1]);
  IntegRawYieldKPlus->SetMarkerColor(Color[1]);
  IntegRawYieldKPlus->SetLineColor(Color[1]);
  IntegRawYieldKPlus->DrawCopy("E1X0Psame");
  lIntegRawYield->AddEntry(IntegRawYieldKPlus,"K^{+}","p");
  // Proton
  IntegRawYieldProton->SetMarkerStyle(Marker[2]);
  IntegRawYieldProton->SetMarkerColor(Color[2]);
  IntegRawYieldProton->SetLineColor(Color[2]);
  IntegRawYieldProton->DrawCopy("E1X0Psame");
  lIntegRawYield->AddEntry(IntegRawYieldProton,"p","p");
  // pi-
  IntegRawYieldPiMinus->SetMarkerStyle(Marker[3]);
  IntegRawYieldPiMinus->SetMarkerColor(Color[0]);
  IntegRawYieldPiMinus->SetLineColor(Color[0]);
  IntegRawYieldPiMinus->DrawCopy("E1X0Psame");
  lIntegRawYield->AddEntry(IntegRawYieldPiMinus,"#pi^{-}","p");
  // K-
  IntegRawYieldKMinus->SetMarkerStyle(Marker[4]);
  IntegRawYieldKMinus->SetMarkerColor(Color[1]);
  IntegRawYieldKMinus->SetLineColor(Color[1]);
  IntegRawYieldKMinus->DrawCopy("E1X0Psame");
  lIntegRawYield->AddEntry(IntegRawYieldKMinus,"K^{-}","p");
  // Antiproton
  IntegRawYieldAntiproton->SetMarkerStyle(Marker[5]);
  IntegRawYieldAntiproton->SetMarkerColor(Color[2]);
  IntegRawYieldAntiproton->SetLineColor(Color[2]);
  IntegRawYieldAntiproton->SetTitle("Raw Yield, Integrated over all p_{T}");
  IntegRawYieldAntiproton->DrawCopy("E1X0Psame");
  lIntegRawYield->AddEntry(IntegRawYieldAntiproton,"#bar{p}","p");
  // write the canvas to the file
  lIntegRawYield->DrawClone();
  fout->cd();
  cIntegRawYieldStability->Write();

  //------------------------------------------------------
  // Plot the efficiency as a function of the run number -
  //------------------------------------------------------
  if (useMC)
    {
      TCanvas * cEfficiencyRun = new TCanvas("cEfficiencyRun","cEfficiencyRun",200,100,700,500);
      TH2F * hAxesEfficiency = new TH2F("hAxesEfficiency","MC Correction Factor at p_{T} = 0.9 GeV/c",nRuns,0,nRuns,100,0,1);
      for (Int_t irun=0; irun<nRuns; irun++)
	hAxesEfficiency->GetXaxis()->SetBinLabel(irun+1, Form("%i",runs[irun]));
      hAxesEfficiency->SetStats(kFALSE);
      hAxesEfficiency->DrawCopy();
      TLegend * lEfficiency = new TLegend(.79,.69,.99,.99);
      // pi+
      EfficiencyPiPlus->SetMarkerStyle(Marker[0]);
      EfficiencyPiPlus->SetMarkerColor(Color[0]);
      EfficiencyPiPlus->SetLineColor(Color[0]);
      EfficiencyPiPlus->DrawCopy("E1X0Psame");
      lEfficiency->AddEntry(EfficiencyPiPlus,"#pi^{+}","lpe");
      // K+
      EfficiencyKPlus->SetMarkerStyle(Marker[1]);
      EfficiencyKPlus->SetMarkerColor(Color[1]);
      EfficiencyKPlus->SetLineColor(Color[1]);
      EfficiencyKPlus->DrawCopy("E1X0Psame");
      lEfficiency->AddEntry(EfficiencyKPlus,"K^{+}","p");
      // Proton
      EfficiencyProton->SetMarkerStyle(Marker[2]);
      EfficiencyProton->SetMarkerColor(Color[2]);
      EfficiencyProton->SetLineColor(Color[2]);
      EfficiencyProton->DrawCopy("E1X0Psame");
      lEfficiency->AddEntry(EfficiencyProton,"p","lpe");
      // pi-
      EfficiencyPiMinus->SetMarkerStyle(Marker[3]);
      EfficiencyPiMinus->SetMarkerColor(Color[0]);
      EfficiencyPiMinus->SetLineColor(Color[0]);
      EfficiencyPiMinus->DrawCopy("E1X0Psame");
      lEfficiency->AddEntry(EfficiencyPiMinus,"#pi^{-}","lpe");
      // K-
      EfficiencyKMinus->SetMarkerStyle(Marker[4]);
      EfficiencyKMinus->SetMarkerColor(Color[1]);
      EfficiencyKMinus->SetLineColor(Color[1]);
      EfficiencyKMinus->DrawCopy("E1X0Psame");
      lEfficiency->AddEntry(EfficiencyKMinus,"K^{-}","lpe");
      // Antiproton
      EfficiencyAntiproton->SetMarkerStyle(Marker[5]);
      EfficiencyAntiproton->SetMarkerColor(Color[2]);
      EfficiencyAntiproton->SetLineColor(Color[2]);
      EfficiencyAntiproton->DrawCopy("E1X0Psame");
      lEfficiency->AddEntry(EfficiencyAntiproton,"#bar{p}","lpe");
      // draw legend
      lEfficiency->SetFillColor(0);
      lEfficiency->DrawClone();
      // write the canvas to file
      fout->cd();
      cEfficiencyRun->Write();
      
      // write the last remaining canvas - cEfficienciesAllRuns - to file
      fout->cd();
      cEfficienciesAllRuns->Write();
    }

  //---------------------------------------------------------------
  // Plot the matching efficiency as a function of the run number -
  //---------------------------------------------------------------
  TCanvas * cMatchingEfficiency = new TCanvas("cMatchingEfficiency","cMatchingEfficiency",250,125,700,500);
  cMatchingEfficiency->cd();
  TH2F * hAxesMatchEff = new TH2F("hAxesMatchEff","",nRuns,0,nRuns,100,0,1);
  hAxesMatchEff->SetTitle("TOF Matching Efficiency at p_{T} = 0.9 GeV/c;;");
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesMatchEff->GetXaxis()->SetBinLabel(irun+1, Form("%i",runs[irun]));
  hAxesMatchEff->SetStats(kFALSE);
  hAxesMatchEff->DrawCopy();
  TLegend * lMatchEff = new TLegend(.7,.75,.85,.85);
  lMatchEff->SetFillColor(0);
  // positive particles
  MatchEffPos->SetMarkerStyle(Marker[0]);
  MatchEffPos->SetMarkerColor(kBlue);
  MatchEffPos->SetLineColor(kBlue);
  MatchEffPos->DrawCopy("E1X0Psame");
  lMatchEff->AddEntry(MatchEffPos,"Positive Particles","pe");
  // negative particles
  MatchEffNeg->SetMarkerStyle(Marker[1]);
  MatchEffNeg->SetMarkerColor(kRed);
  MatchEffNeg->SetLineColor(kRed);
  MatchEffNeg->DrawCopy("E1X0Psame");
  lMatchEff->AddEntry(MatchEffNeg,"Negative Particles","pe");
  // write to file
  lMatchEff->DrawClone();
  fout->cd();
  cMatchingEfficiency->Write();

  /*
  //-------------------------
  // nsigma Projection fits -
  //-------------------------
  TH1F * TPCnsigMeanPionProj = new TH1F("TPCnsigMeanPionProj","",100,-1,1);
  TH1F * TPCnsigMeanKaonProj = new TH1F("TPCnsigMeanKaonProj","",100,-1,1);
  TH1F * TPCnsigMeanProtonProj = new TH1F("TPCnsigMeanProtonProj","",100,-1,1);
  TH1F * TPCnsigSigmaPionProj = new TH1F("TPCnsigSigmaPionProj","",100,0,2);
  TH1F * TPCnsigSigmaKaonProj = new TH1F("TPCnsigSigmaKaonProj","",100,0,2);
  TH1F * TPCnsigSigmaProtonProj = new TH1F("TPCnsigSigmaProtonProj","",100,0,2);
  TH1F * TOFnsigMeanPionProj = new TH1F("TOFnsigMeanPionProj","",100,-1,1);
  TH1F * TOFnsigMeanKaonProj = new TH1F("TOFnsigMeanKaonProj","",100,-1,1);
  TH1F * TOFnsigMeanProtonProj = new TH1F("TOFnsigMeanProtonProj","",100,-1,1);
  TH1F * TOFnsigSigmaPionProj = new TH1F("TOFnsigSigmaPionProj","",100,0,2);
  TH1F * TOFnsigSigmaKaonProj = new TH1F("TOFnsigSigmaKaonProj","",100,0,2);
  TH1F * TOFnsigSigmaProtonProj = new TH1F("TOFnsigSigmaProtonProj","",100,0,2);
  for (Int_t irun=0; irun<nRuns; irun++)
    {
      TPCnsigMeanPionProj->Fill(TPCnsigMeanTrendPion->GetBinContent(irun+1));
      TPCnsigMeanKaonProj->Fill(TPCnsigMeanTrendKaon->GetBinContent(irun+1));
      TPCnsigMeanProtonProj->Fill(TPCnsigMeanTrendProton->GetBinContent(irun+1));
      TPCnsigSigmaPionProj->Fill(TPCnsigSigmaTrendPion->GetBinContent(irun+1));
      TPCnsigSigmaKaonProj->Fill(TPCnsigSigmaTrendKaon->GetBinContent(irun+1));
      TPCnsigSigmaProtonProj->Fill(TPCnsigSigmaTrendProton->GetBinContent(irun+1));
      TOFnsigMeanPionProj->Fill(TOFnsigMeanTrendPion->GetBinContent(irun+1));
      TOFnsigMeanKaonProj->Fill(TOFnsigMeanTrendKaon->GetBinContent(irun+1));
      TOFnsigMeanProtonProj->Fill(TOFnsigMeanTrendProton->GetBinContent(irun+1));
      TOFnsigSigmaPionProj->Fill(TOFnsigSigmaTrendPion->GetBinContent(irun+1));
      TOFnsigSigmaKaonProj->Fill(TOFnsigSigmaTrendKaon->GetBinContent(irun+1));
      TOFnsigSigmaProtonProj->Fill(TOFnsigSigmaTrendProton->GetBinContent(irun+1));
    }
  TCanvas * cNSigTrendProj = new TCanvas("cNSigTrendProj","cNSigTrendProj",300,150,700,500);
  cNSigTrendProj->Divide(2,2);
  // TPC Means
  cNSigTrendProj->cd(1);
  TH2F * hAxesTPCNSigMeanProj = new TH2F("hAxesTPCNSigMeanProj","",100,-1,1,100,0,80);
  hAxesTPCNSigMeanProj->SetStats(kFALSE);
  hAxesTPCNSigMeanProj->GetXaxis()->SetTitle("TPC nsigma Means");
  hAxesTPCNSigMeanProj->GetYaxis()->SetTitle("Number of Runs");
  hAxesTPCNSigMeanProj->DrawCopy();
  TPCnsigMeanPionProj->SetMarkerStyle(Marker[0]);
  TPCnsigMeanPionProj->SetMarkerColor(Color[0]);
  TPCnsigMeanPionProj->SetLineColor(Color[0]);
  TPCnsigMeanPionProj->SetStats(kFALSE);
  TPCnsigMeanPionProj->DrawCopy("LPsame");
  TPCnsigMeanKaonProj->SetMarkerStyle(Marker[1]);
  TPCnsigMeanKaonProj->SetMarkerColor(Color[1]);
  TPCnsigMeanKaonProj->SetLineColor(Color[1]);
  TPCnsigMeanKaonProj->SetStats(kFALSE);
  TPCnsigMeanKaonProj->DrawCopy("LPsame");
  TPCnsigMeanProtonProj->SetMarkerStyle(Marker[1]);
  TPCnsigMeanProtonProj->SetMarkerColor(Color[2]);
  TPCnsigMeanProtonProj->SetLineColor(Color[2]);
  TPCnsigMeanProtonProj->SetStats(kFALSE);
  TPCnsigMeanProtonProj->DrawCopy("LPsame");
  TLegend * lTPCNSigMeanProj = new TLegend(.65,.65,.85,.85);  
  lTPCNSigMeanProj->AddEntry(TPCnsigMeanPionProj,"#pi^{+}, #pi^{-}","lp");
  lTPCNSigMeanProj->AddEntry(TPCnsigMeanKaonProj,"K^{+}, K^{-}","lp");
  lTPCNSigMeanProj->AddEntry(TPCnsigMeanProtonProj,"p, #bar{p}","lp");
  lTPCNSigMeanProj->SetFillColor(0);
  lTPCNSigMeanProj->DrawClone();
  // TPC Sigmas
  cNSigTrendProj->cd(3);
  TH2F * hAxesTPCNSigSigmaProj = new TH2F("hAxesTPCNSigSigmaProj","",100,0,2,100,0,80);
  hAxesTPCNSigSigmaProj->SetStats(kFALSE);
  hAxesTPCNSigSigmaProj->GetXaxis()->SetTitle("TPC nsigma Sigmas");
  hAxesTPCNSigSigmaProj->GetYaxis()->SetTitle("Number of Runs");
  hAxesTPCNSigSigmaProj->DrawCopy();
  TPCnsigSigmaPionProj->SetMarkerStyle(Marker[0]);
  TPCnsigSigmaPionProj->SetMarkerColor(Color[0]);
  TPCnsigSigmaPionProj->SetLineColor(Color[0]);
  TPCnsigSigmaPionProj->SetStats(kFALSE);
  TPCnsigSigmaPionProj->DrawCopy("LPsame");
  TPCnsigSigmaKaonProj->SetMarkerStyle(Marker[1]);
  TPCnsigSigmaKaonProj->SetMarkerColor(Color[1]);
  TPCnsigSigmaKaonProj->SetLineColor(Color[1]);
  TPCnsigSigmaKaonProj->SetStats(kFALSE);
  TPCnsigSigmaKaonProj->DrawCopy("LPsame");
  TPCnsigSigmaProtonProj->SetMarkerStyle(Marker[2]);
  TPCnsigSigmaProtonProj->SetMarkerColor(Color[2]);
  TPCnsigSigmaProtonProj->SetLineColor(Color[2]);
  TPCnsigSigmaProtonProj->SetStats(kFALSE);
  TPCnsigSigmaProtonProj->DrawCopy("LPsame");
  TLegend * lTPCNSigSigmaProj = new TLegend(.65,.65,.85,.85);  
  lTPCNSigSigmaProj->AddEntry(TPCnsigSigmaPionProj,"#pi^{+}, #pi^{-}","lp");
  lTPCNSigSigmaProj->AddEntry(TPCnsigSigmaKaonProj,"K^{+}, K^{-}","lp");
  lTPCNSigSigmaProj->AddEntry(TPCnsigSigmaProtonProj,"p, #bar{p}","lp");
  lTPCNSigSigmaProj->SetFillColor(0);
  lTPCNSigSigmaProj->DrawClone();
  // TOF Means
  cNSigTrendProj->cd(2);
  TH2F * hAxesTOFNSigMeanProj = new TH2F("hAxesTOFNSigMeanProj","",100,-1,1,100,0,80);
  hAxesTOFNSigMeanProj->SetStats(kFALSE);
  hAxesTOFNSigMeanProj->GetXaxis()->SetTitle("TOF nsigma Means");
  hAxesTOFNSigMeanProj->GetYaxis()->SetTitle("Number of Runs");
  hAxesTOFNSigMeanProj->DrawCopy();
  TOFnsigMeanPionProj->SetMarkerStyle(Marker[0]);
  TOFnsigMeanPionProj->SetMarkerColor(Color[0]);
  TOFnsigMeanPionProj->SetLineColor(Color[0]);
  TOFnsigMeanPionProj->SetStats(kFALSE);
  TOFnsigMeanPionProj->DrawCopy("LPsame");
  TOFnsigMeanKaonProj->SetMarkerStyle(Marker[1]);
  TOFnsigMeanKaonProj->SetMarkerColor(Color[1]);
  TOFnsigMeanKaonProj->SetLineColor(Color[1]);
  TOFnsigMeanKaonProj->SetStats(kFALSE);
  TOFnsigMeanKaonProj->DrawCopy("LPsame");
  TOFnsigMeanProtonProj->SetMarkerStyle(Marker[2]);
  TOFnsigMeanProtonProj->SetMarkerColor(Color[2]);
  TOFnsigMeanProtonProj->SetLineColor(Color[2]);
  TOFnsigMeanProtonProj->SetStats(kFALSE);
  TOFnsigMeanProtonProj->DrawCopy("LPsame");
  TLegend * lTOFNSigMeanProj = new TLegend(.65,.65,.85,.85);  
  lTOFNSigMeanProj->AddEntry(TOFnsigMeanPionProj,"#pi^{+}, #pi^{-}","lp");
  lTOFNSigMeanProj->AddEntry(TOFnsigMeanKaonProj,"K^{+}, K^{-}","lp");
  lTOFNSigMeanProj->AddEntry(TOFnsigMeanProtonProj,"p, #bar{p}","lp");
  lTOFNSigMeanProj->SetFillColor(0);
  lTOFNSigMeanProj->DrawClone();
  // TOF Sigmas
  cNSigTrendProj->cd(4);
  TH2F * hAxesTOFNSigSigmaProj = new TH2F("hAxesTOFNSigSigmaProj","",100,0,2,100,0,80);
  hAxesTOFNSigSigmaProj->SetStats(kFALSE);
  hAxesTOFNSigSigmaProj->GetXaxis()->SetTitle("TOF nsigma Sigmas");
  hAxesTOFNSigSigmaProj->GetYaxis()->SetTitle("Number of Runs");
  hAxesTOFNSigSigmaProj->DrawCopy();
  TOFnsigSigmaPionProj->SetMarkerStyle(Marker[0]);
  TOFnsigSigmaPionProj->SetMarkerColor(Color[0]);
  TOFnsigSigmaPionProj->SetLineColor(Color[0]);
  TOFnsigSigmaPionProj->SetStats(kFALSE);
  TOFnsigSigmaPionProj->DrawCopy("LPsame");
  TOFnsigSigmaKaonProj->SetMarkerStyle(Marker[1]);
  TOFnsigSigmaKaonProj->SetMarkerColor(Color[1]);
  TOFnsigSigmaKaonProj->SetLineColor(Color[1]);
  TOFnsigSigmaKaonProj->SetStats(kFALSE);
  TOFnsigSigmaKaonProj->DrawCopy("LPsame");
  TOFnsigSigmaProtonProj->SetMarkerStyle(Marker[2]);
  TOFnsigSigmaProtonProj->SetMarkerColor(Color[2]);
  TOFnsigSigmaProtonProj->SetLineColor(Color[2]);
  TOFnsigSigmaProtonProj->SetStats(kFALSE);
  TOFnsigSigmaProtonProj->DrawCopy("LPsame");
  TLegend * lTOFNSigSigmaProj = new TLegend(.65,.65,.85,.85);  
  lTOFNSigSigmaProj->AddEntry(TOFnsigSigmaPionProj,"#pi^{+}, #pi^{-}","lp");
  lTOFNSigSigmaProj->AddEntry(TOFnsigSigmaKaonProj,"K^{+}, K^{-}","lp");
  lTOFNSigSigmaProj->AddEntry(TOFnsigSigmaProtonProj,"p, #bar{p}","lp");
  lTOFNSigSigmaProj->SetFillColor(0);
  lTOFNSigSigmaProj->DrawClone();
  // write to file
  fout->cd();
  cNSigTrendProj->Write();

  //-----------------------
  // Integrated Raw Yield -
  //-----------------------
  TH1F * IntegRawYieldProj = new TH1F("IntegRawYieldProj","",50,7500,9000);
  for (Int_t irun=0; irun<nRuns; irun++)
    {
      IntegRawYieldProj->Fill(IntegRawYieldAll->GetBinContent(irun+1));
    }
  TCanvas * cIntegRawYieldProj = new TCanvas("cIntegRawYieldProj","cIntegRawYieldProj",350,175,700,500);
  cIntegRawYieldProj->cd();
  IntegRawYieldProj->SetMarkerStyle(20);
  IntegRawYieldProj->SetMarkerColor(kBlack);
  IntegRawYieldProj->SetLineColor(kBlack);
  IntegRawYieldProj->SetStats(kFALSE);
  IntegRawYieldProj->GetXaxis()->SetTitle("Integrated Raw Yield (all particles)");
  IntegRawYieldProj->GetYaxis()->SetTitle("Number of Runs");
  IntegRawYieldProj->DrawCopy("hist");
  // write to file
  fout->cd();
  cIntegRawYieldProj->Write();

  //---------------
  // Efficiencies -
  //---------------
  if (useMC)
    {
      TH1F * EffPiPlusProj = new TH1F("EffPiPlusProj","",100,0,1);
      TH1F * EffKPlusProj = new TH1F("EffKPlusProj","",100,0,1);
      TH1F * EffProtonProj = new TH1F("EffProtonProj","",100,0,1);
      TH1F * EffPiMinusProj = new TH1F("EffPiMinusProj","",100,0,1);
      TH1F * EffKMinusProj = new TH1F("EffKMinusProj","",100,0,1);
      TH1F * EffAntiprotonProj = new TH1F("EffAntiprotonProj","",100,0,1);
      for (Int_t irun=0; irun<nRuns; irun++)
	{
	  EffPiPlusProj->Fill(EfficiencyPiPlus->GetBinContent(irun+1));
	  EffKPlusProj->Fill(EfficiencyKPlus->GetBinContent(irun+1));
	  EffProtonProj->Fill(EfficiencyProton->GetBinContent(irun+1));
	  EffPiMinusProj->Fill(EfficiencyPiMinus->GetBinContent(irun+1));
	  EffKMinusProj->Fill(EfficiencyKMinus->GetBinContent(irun+1));
	  EffAntiprotonProj->Fill(EfficiencyAntiproton->GetBinContent(irun+1));
	}
      TCanvas * cEffProj = new TCanvas("cEffProj","cEffProj",400,200,700,500);
      cEffProj->cd();
      TH2F * hAxesEffProj = new TH2F("hAxesEffProj","",100,0.35,1,100,0,60);
      hAxesEffProj->SetStats(kFALSE);
      hAxesEffProj->GetXaxis()->SetTitle("Efficiency");
      hAxesEffProj->GetYaxis()->SetTitle("Number of Runs");
      hAxesEffProj->DrawCopy();
      EffPiPlusProj->SetLineColor(Color[0]);
      EffPiPlusProj->SetMarkerStyle(Marker[0]);
      EffPiPlusProj->SetMarkerColor(Color[0]);
      EffPiPlusProj->SetStats(kFALSE);
      EffPiPlusProj->DrawCopy("LPSame");
      EffKPlusProj->SetLineColor(Color[1]);
      EffKPlusProj->SetMarkerStyle(Marker[1]);
      EffKPlusProj->SetMarkerColor(Color[1]);
      EffKPlusProj->SetStats(kFALSE);
      EffKPlusProj->DrawCopy("LPSame");
      EffProtonProj->SetLineColor(Color[2]);
      EffProtonProj->SetMarkerStyle(Marker[2]);
      EffProtonProj->SetMarkerColor(Color[2]);
      EffProtonProj->SetStats(kFALSE);
      EffProtonProj->DrawCopy("LPSame");
      EffPiMinusProj->SetLineColor(Color[0]);
      EffPiMinusProj->SetMarkerStyle(Marker[3]);
      EffPiMinusProj->SetMarkerColor(Color[0]);
      EffPiMinusProj->SetStats(kFALSE);
      EffPiMinusProj->DrawCopy("LPSame");
      EffKMinusProj->SetLineColor(Color[1]);
      EffKMinusProj->SetMarkerStyle(Marker[4]);
      EffKMinusProj->SetMarkerColor(Color[1]);
      EffKMinusProj->SetStats(kFALSE);
      EffKMinusProj->DrawCopy("LPSame");
      EffAntiprotonProj->SetLineColor(Color[2]);
      EffAntiprotonProj->SetMarkerStyle(Marker[5]);
      EffAntiprotonProj->SetMarkerColor(Color[2]);
      EffAntiprotonProj->SetStats(kFALSE);
      EffAntiprotonProj->DrawCopy("LPSame");
      TLegend * lEffProj = new TLegend(.65,.65,.85,.85);
      lEffProj->AddEntry(EffPiPlusProj,"#pi^{+}","lp");
      lEffProj->AddEntry(EffKPlusProj,"K^{+}","lp");
      lEffProj->AddEntry(EffProtonProj,"lp","lp");
      lEffProj->AddEntry(EffPiMinusProj,"#pi^{-}","lp");
      lEffProj->AddEntry(EffKMinusProj,"K^{-}","lp");
      lEffProj->AddEntry(EffAntiprotonProj,"#bar{p}","lp");
      lEffProj->SetFillColor(0);
      lEffProj->DrawClone();
      // save to file
      fout->cd();
      cEffProj->Write();
    }
  
  //------------------------
  // Matching Efficiencies -
  //------------------------
  // project the matching efficiency histograms onto the y-axis
  TH1F * MatchEffPosProj = new TH1F("MatchEffPosProj","",100,0.54,0.68);
  TH1F * MatchEffNegProj = new TH1F("MatchEffNegProj","",100,0.54,0.68);
  for (Int_t irun=0; irun<nRuns; irun++)
    {
      MatchEffPosProj->Fill(MatchEffPos->GetBinContent(irun+1));
      MatchEffNegProj->Fill(MatchEffNeg->GetBinContent(irun+1));
    }
  TCanvas * cMatchEffProj = new TCanvas("cMatchEffProj","cMatchEffProj",450,225,700,500);
  cMatchEffProj->cd();
  MatchEffPosProj->SetMarkerStyle(20);
  MatchEffPosProj->SetMarkerColor(kBlue);
  MatchEffPosProj->SetLineColor(kBlue);
  MatchEffPosProj->SetStats(kFALSE);
  MatchEffPosProj->GetXaxis()->SetTitle("Matching Efficiency");
  MatchEffPosProj->GetYaxis()->SetTitle("Number of Runs");
  MatchEffPosProj->GetYaxis()->SetRangeUser(0,50);
  MatchEffPosProj->DrawCopy("LP");
  MatchEffNegProj->SetMarkerStyle(21);
  MatchEffNegProj->SetMarkerColor(kRed);
  MatchEffNegProj->SetLineColor(kRed);
  MatchEffNegProj->SetStats(kFALSE);
  MatchEffNegProj->DrawCopy("LPsame");
  TLegend * lMatchEffProj = new TLegend(0.75,0.75,0.85,0.85);
  lMatchEffProj->AddEntry(MatchEffPosProj,"Pos","lp");
  lMatchEffProj->AddEntry(MatchEffNegProj,"Neg","lp");
  lMatchEffProj->SetFillColor(0);
  lMatchEffProj->DrawClone();
  // save to file
  fout->cd();
  cMatchEffProj->Write();
  */

  Printf("\n\n\n--- Leaving function PlotAndSave() ---\n\n\n");
}
