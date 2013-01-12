//////////////////////////////////////////////////////////////////////////////////
// CheckMatchingEfficiency.C (called by AODQAChecks.C)                          //
//                                                                              //
// written by John Groh                                                         //
//////////////////////////////////////////////////////////////////////////////////

void CheckMatchingEfficiency(AliSpectraAODTrackCuts * tcuts,
			     TH1F*& MatchEffPos,
			     TH1F*& MatchEffNeg,
			     Float_t FixedPtMatchEff,
			     Int_t runs[],
                             Int_t irun,
			     Int_t nRuns,
			     Bool_t useMC)
{

  // get matching efficiency histograms
  TH1F * hMatchEffPos = (TH1F*)tcuts->GetHistoNMatchedPos();
  hMatchEffPos->Sumw2();
  hMatchEffPos->Divide(hMatchEffPos,(TH1F*)tcuts->GetHistoNSelectedPos()->Clone(),1,1,"B"); // binomial error!  
  TH1F * hMatchEffNeg = (TH1F*)tcuts->GetHistoNMatchedNeg();
  hMatchEffNeg->Sumw2();
  hMatchEffNeg->Divide(hMatchEffNeg,(TH1F*)tcuts->GetHistoNSelectedNeg()->Clone(),1,1,"B"); // binomial error!  

  // get values for a fixed pt
  Int_t nMatchedPos = hMatchEffPos->GetBinContent(hMatchEffPos->FindBin(FixedPtMatchEff));
  Int_t nSelectedPos = tcuts->GetHistoNSelectedPos()->GetBinContent(tcuts->GetHistoNSelectedPos()->FindBin(FixedPtMatchEff));
  Int_t nMatchedNeg = tcuts->GetHistoNMatchedNeg()->GetBinContent(tcuts->GetHistoNMatchedNeg()->FindBin(FixedPtMatchEff));
  Int_t nSelectedNeg = tcuts->GetHistoNSelectedNeg()->GetBinContent(tcuts->GetHistoNSelectedNeg()->FindBin(FixedPtMatchEff));

  // fill matching eff vs. run number histograms
  MatchEffPos->SetBinContent(irun+1, hMatchEffPos->GetBinContent(hMatchEffPos->FindBin(FixedPtMatchEff)));
  MatchEffPos->SetBinError(irun+1, hMatchEffPos->GetBinError(hMatchEffPos->FindBin(FixedPtMatchEff)));
  MatchEffNeg->SetBinContent(irun+1, hMatchEffNeg->GetBinContent(hMatchEffNeg->FindBin(FixedPtMatchEff)));
  MatchEffNeg->SetBinError(irun+1, hMatchEffNeg->GetBinError(hMatchEffNeg->FindBin(FixedPtMatchEff)));

  // plot individual runs temporarily before saving to a pdf
  TCanvas * cMatchEffIndiv = new TCanvas("cMatchEffIndiv","cMatchEffIndiv");
  TH2F * hAxesMatchEffIndiv = new TH2F("hAxesMatchEffIndiv","",100,0,5,100,0,1);
  hAxesMatchEffIndiv->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hAxesMatchEffIndiv->GetYaxis()->SetTitle("TOF Matching Efficiency");
  hAxesMatchEffIndiv->SetStats(kFALSE);
  hAxesMatchEffIndiv->DrawCopy();
  // Pos
  hMatchEffPos->SetStats(kFALSE);
  hMatchEffPos->SetLineColor(kBlue);
  hMatchEffPos->SetMarkerColor(kBlue);
  hMatchEffPos->SetMarkerStyle(20);
  hMatchEffPos->SetMarkerSize(0.2);
  hMatchEffPos->DrawCopy("E1PSame");
  // Neg
  hMatchEffNeg->SetStats(kFALSE);
  hMatchEffNeg->SetLineColor(kRed);
  hMatchEffNeg->SetMarkerColor(kRed);
  hMatchEffNeg->SetMarkerStyle(21);
  hMatchEffNeg->SetMarkerSize(0.2);
  hMatchEffNeg->DrawCopy("E1PSame");
  TLegend * lMatchEffIndiv = new TLegend(.65,.15,.85,.35);
  lMatchEffIndiv->AddEntry(hMatchEffPos,"Pos","pe");
  lMatchEffIndiv->AddEntry(hMatchEffNeg,"Neg","pe");
  lMatchEffIndiv->AddEntry(hMatchEffPos,Form("Run %i",runs[irun]),"");
  if (useMC) lMatchEffIndiv->AddEntry(hMatchEffPos,"MC","");
  else lMatchEffIndiv->AddEntry(hMatchEffPos,"DATA","");
  lMatchEffIndiv->SetFillColor(0);
  lMatchEffIndiv->DrawClone();

  // save to a pdf and close the temporary canvas
  if (useMC)
    {
      if (irun == 0) cMatchEffIndiv->SaveAs("Plots/MC/AODMatchingEfficiencies.pdf(","pdf");
      else if (irun < nRuns-1) cMatchEffIndiv->SaveAs("Plots/MC/AODMatchingEfficiencies.pdf","pdf");
      else if (irun == nRuns-1) cMatchEffIndiv->SaveAs("Plots/MC/AODMatchingEfficiencies.pdf)","pdf");
      cMatchEffIndiv->Close();
    }
  else
    {
      if (irun == 0) cMatchEffIndiv->SaveAs("Plots/DATA/AODMatchingEfficiencies.pdf(","pdf");
      else if (irun < nRuns-1) cMatchEffIndiv->SaveAs("Plots/DATA/AODMatchingEfficiencies.pdf","pdf");
      else if (irun == nRuns-1) cMatchEffIndiv->SaveAs("Plots/DATA/AODMatchingEfficiencies.pdf)","pdf");
      cMatchEffIndiv->Close();
    }
}





