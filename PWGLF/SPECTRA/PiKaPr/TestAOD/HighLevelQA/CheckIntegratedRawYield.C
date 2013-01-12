//////////////////////////////////////////////////////////////////////////////
// CheckIntegratedRawYield.C (called by AODQAChecks.C)                      //
//                                                                          //
// Written by John Groh                                                     //
//////////////////////////////////////////////////////////////////////////////

void CheckIntegratedRawYield(AliSpectraAODEventCuts * ecuts,
			     AliSpectraAODHistoManager * hman,
			     TH1F*& IntegRawYieldAll,
			     TH1F*& IntegRawYieldPiPlus,
			     TH1F*& IntegRawYieldPiMinus,
			     TH1F*& IntegRawYieldKPlus,
			     TH1F*& IntegRawYieldKMinus,
			     TH1F*& IntegRawYieldProton,
			     TH1F*& IntegRawYieldAntiproton,
			     Int_t runs[],
			     Int_t irun,
			     Int_t nRuns,
			     TString names[],
			     Bool_t useMC)
{
  // get total raw yield histogram
  TH1F * hRaw_allCh = (TH1F*)((TH1F*) hman->GetPtHistogram1D("hHistPtRec",-1,-1))->Clone();
  // scale it by the number of events
  hRaw_allCh->Scale(1./ecuts->NumberOfEvents(), "width");
  
  // integrate and fill IntegRawYieldAll
  IntegRawYieldAll->SetBinContent(irun+1,hRaw_allCh->Integral());
  
  // get particle raw yield histograms
  TH1F * hRaw[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  TString hname = Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  hRaw[index] = (TH1F*)((TH1F*)hman->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
	  hRaw[index]->Scale(1./ecuts->NumberOfEvents(),"width");
	}
    }
  
  // integrate and fill IntegRawYield___ (for individual particles)
  IntegRawYieldPiPlus->Fill(irun,hRaw[0]->Integral());
  IntegRawYieldKPlus->Fill(irun,hRaw[1]->Integral());
  IntegRawYieldProton->Fill(irun,hRaw[2]->Integral());
  IntegRawYieldPiMinus->Fill(irun,hRaw[3]->Integral());
  IntegRawYieldKMinus->Fill(irun,hRaw[4]->Integral());
  IntegRawYieldAntiproton->Fill(irun,hRaw[5]->Integral());
  
  // temporarily draw individual runs to a canvas
  TCanvas * cRawYieldIndiv = new TCanvas("cRawYieldIndiv","cRawYieldIndiv");
  gPad->SetLogy();
  TLegend * lRawYieldIndiv = new TLegend(0.69,.59,.99,.99);
  lRawYieldIndiv->SetFillColor(0);
  hRaw_allCh->SetTitle("Raw Yield;p_{T} (GeV/c);");
  hRaw_allCh->GetXaxis()->CenterTitle();
  hRaw_allCh->SetStats(kFALSE);
  hRaw_allCh->SetMarkerStyle(34);
  hRaw_allCh->SetMarkerColor(kGreen);
  hRaw_allCh->SetLineColor(kGreen);
  hRaw_allCh->GetYaxis()->SetRangeUser(.001,1000);
  hRaw_allCh->DrawCopy("E1P");
  lRawYieldIndiv->AddEntry(hRaw_allCh,Form("Run %i",runs[irun]),"");
  if (useMC) lRawYieldIndiv->AddEntry(hRaw_allCh,"MC","");
  else lRawYieldIndiv->AddEntry(hRaw_allCh,"DATA","");
  lRawYieldIndiv->AddEntry(hRaw_allCh,"All Particles","pe");

  Int_t index = 0;
  Int_t colorIndex;
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      colorIndex = 0; // only 3 colors, so have to reset after inner loop runs
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  hRaw[index]->SetStats(kFALSE);
	  hRaw[index]->SetMarkerStyle(Marker[index]);
	  hRaw[index]->SetMarkerColor(Color[colorIndex]);
	  hRaw[index]->SetLineColor(Color[colorIndex]);
	  hRaw[index]->DrawCopy("E1Psame");
	  lRawYieldIndiv->AddEntry(hRaw[index],Names[index].Data(),"pe");
	  colorIndex++;
	  index++;
	}
    }
  lRawYieldIndiv->DrawClone();

  // close the canvas and save to a pdf file
  if (useMC)
    {
      if (irun == 0) cRawYieldIndiv->SaveAs("Plots/MC/RawYieldIndiv.pdf(","pdf");
      else if (irun < nRuns-1) cRawYieldIndiv->SaveAs("Plots/MC/RawYieldIndiv.pdf","pdf");
      else if (irun == nRuns-1) cRawYieldIndiv->SaveAs("Plots/MC/RawYieldIndiv.pdf)","pdf");
      cRawYieldIndiv->Close();
    }
  else
    {
      if (irun == 0) cRawYieldIndiv->SaveAs("Plots/DATA/RawYieldIndiv.pdf(","pdf");
      else if (irun < nRuns-1) cRawYieldIndiv->SaveAs("Plots/DATA/RawYieldIndiv.pdf","pdf");
      else if (irun == nRuns-1) cRawYieldIndiv->SaveAs("Plots/DATA/RawYieldIndiv.pdf)","pdf");
      cRawYieldIndiv->Close();
    }
}

