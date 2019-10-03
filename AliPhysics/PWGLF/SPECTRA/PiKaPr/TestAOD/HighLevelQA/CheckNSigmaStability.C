///////////////////////////////////////////////////////////////////////////////////////////
// CheckNSigmaStability.C (called by AODQAChecks.C)                                      //
//                                                                                       //
// Written by John Groh                                                                  //
///////////////////////////////////////////////////////////////////////////////////////////

void CheckNSigmaStability(AliSpectraAODHistoManager * hman,
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
			  Int_t runs[],
			  Int_t nRuns,
			  Int_t irun,
			  Bool_t useMC)
{
  // ranges for the projections
  const Double_t ProjTPC[2] = {0.4, 0.5};
  const Double_t ProjTOF[2] = {0.9, 1.0};

  // ranges for the gaussian fits
  const Double_t fitRangeTPC[nPart][2] = {{-3.5, 3.5},
					  {-1.4, 1.4},
					  {-2.5, 2.5}};
  const Double_t fitRangeTOF[nPart][2] = {{-2.3, 2.3},
					  {-2.0, 2.0},
					  {-2.0, 2.0}};

  //------------------------------------------------------
  //                        TPC                          -
  //------------------------------------------------------

  // canvas for printing the projections and fits to a pdf file (once per run)
  TCanvas * cNSigProjFits = new TCanvas("cNSigProjFits","cNSigProjFits");
  cNSigProjFits->Divide(3,2);

  // define the 2D nsigma histograms, the projections, and the fitting functions
  // (one each for pions, kaons, and protons)
  TH2F * nsig_TPC[nPart];
  TH1F * nsig_TPC_proj[nPart];
  TF1 * fitFuncTPC[nPart];
  for (Int_t ipart; ipart<nPart; ipart++)
    {
      nsig_TPC[ipart] = new TH2F;
      nsig_TPC_proj[ipart] = new TH1F;
    }

  // for plotting data/fit
  TGraph gDataOverFitTPC[nPart];
  TGraph gDataOverFitTOF[nPart];
  TCanvas * cDataOverFit = new TCanvas("cDataOverFit","cDataOverFit");
  cDataOverFit->Divide(2,1);
   
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      // create projection
      nsig_TPC[ipart] = (TH2F*)((TH2F*)hman->GetNSigHistogram(Form("hHistNSig%sTPC",Particle[ipart].Data())))->Clone();
      if (useMC)
	nsig_TPC_proj[ipart] = (TH1F*)nsig_TPC[ipart]->ProjectionY(Form("TPC NsigProjection %s, mc [%.1f,%.1f]",Particle[ipart].Data(),ProjTPC[0],ProjTPC[1]), nsig_TPC[ipart]->GetXaxis()->FindBin(ProjTPC[0]), nsig_TPC[ipart]->GetXaxis()->FindBin(ProjTPC[1]));
      else
	nsig_TPC_proj[ipart] = (TH1F*)nsig_TPC[ipart]->ProjectionY(Form("TPC NsigProjection %s, data [%.1f,%.1f]",Particle[ipart].Data(),ProjTPC[0],ProjTPC[1]), nsig_TPC[ipart]->GetXaxis()->FindBin(ProjTPC[0]), nsig_TPC[ipart]->GetXaxis()->FindBin(ProjTPC[1]));

      // fit the peak of interest with a gaussian
      fitFuncTPC[ipart] = new TF1("fitFuncTPC","gaus",fitRangeTPC[ipart][0],fitRangeTPC[ipart][1]);
      nsig_TPC_proj[ipart]->Fit(fitFuncTPC[ipart],"NRLQ");
    
      // draw the projections and fits
      cNSigProjFits->cd(ipart+1);
      gPad->SetLogy();
      nsig_TPC_proj[ipart]->GetXaxis()->SetTitle(Form("TPC nsigma for %.1f GeV/c < p < %.1f GeV/c",ProjTPC[0],ProjTPC[1]));
      nsig_TPC_proj[ipart]->GetXaxis()->SetRangeUser(-10,10);
      nsig_TPC_proj[ipart]->SetStats(kFALSE);
      nsig_TPC_proj[ipart]->DrawCopy();
      //fitFuncTPC[ipart]->SetLineWidth(1);
      fitFuncTPC[ipart]->DrawCopy("same");
      TLegend * lNSigTPC = new TLegend(0.59,0.59,0.99,0.99);
      lNSigTPC->SetFillColor(0);
      lNSigTPC->AddEntry(nsig_TPC_proj[ipart],Form("%ss, TPC",Particle[ipart].Data()),"");
      lNSigTPC->AddEntry(nsig_TPC_proj[ipart],Form("Run %i",runs[irun]),"");
      if (useMC) lNSigTPC->AddEntry(nsig_TPC_proj[ipart],"MC","");
      else lNSigTPC->AddEntry(nsig_TPC_proj[ipart],"DATA","");
      //      lNSigTPC->AddEntry(nsig_TPC_proj[ipart],Form("#chi^{2}/nDOF = %.2f",(Float_t)resultTPC->Chi2() / (Float_t)resultTPC->Ndf()),"");
      lNSigTPC->DrawClone();

      // fill gDataOverFitTPC[] graphs
      Int_t nGraphPoints = 0;
      for (Float_t iPt = fitRangeTPC[ipart][0]; iPt <= fitRangeTPC[ipart][1]; iPt += (fitRangeTPC[ipart][1] - fitRangeTPC[ipart][0])/15)
	{
	  gDataOverFitTPC[ipart].SetPoint(nGraphPoints,iPt, nsig_TPC_proj[ipart]->GetBinContent(nsig_TPC_proj[ipart]->FindBin(iPt)) / fitFuncTPC[ipart]->Eval(iPt));
	  nGraphPoints++;
	}

      // fill the histograms with the fit parameters and run numbers
      switch (ipart)
	{
	case 0: // pion
	  TPCnsigMeanTrendPion->SetBinContent(irun+1, fitFuncTPC[ipart]->GetParameter(1));
	  TPCnsigSigmaTrendPion->SetBinContent(irun+1, fitFuncTPC[ipart]->GetParameter(2));	 
	  TPCnsigMeanTrendPion->SetBinError(irun+1, fitFuncTPC[ipart]->GetParError(1));
	  TPCnsigSigmaTrendPion->SetBinError(irun+1, fitFuncTPC[ipart]->GetParError(2));
	  break;
	case 1: // kaon
	  TPCnsigMeanTrendKaon->SetBinContent(irun+1, fitFuncTPC[ipart]->GetParameter(1));
	  TPCnsigSigmaTrendKaon->SetBinContent(irun+1, fitFuncTPC[ipart]->GetParameter(2));
	  TPCnsigMeanTrendKaon->SetBinError(irun+1, fitFuncTPC[ipart]->GetParError(1));
	  TPCnsigSigmaTrendKaon->SetBinError(irun+1, fitFuncTPC[ipart]->GetParError(2));
	  break;
	case 2: // proton
	  TPCnsigMeanTrendProton->SetBinContent(irun+1, fitFuncTPC[ipart]->GetParameter(1));
	  TPCnsigSigmaTrendProton->SetBinContent(irun+1, fitFuncTPC[ipart]->GetParameter(2));
	  TPCnsigMeanTrendProton->SetBinError(irun+1, fitFuncTPC[ipart]->GetParError(1));
	  TPCnsigSigmaTrendProton->SetBinError(irun+1, fitFuncTPC[ipart]->GetParError(2));
	  break;
	default:
	  Printf("\n!!! ERROR in TPC switch control structure in CheckNSigmaStability.C !!!\n");
	  break;
	  }
    } // end loop over ipart

  // draw the data/fit
  cDataOverFit->cd(1);
  TH2F * hAxesDataOverFitTPC = new TH2F("hAxesDataOverFitTPC","",100,-5,5,100,-3,5);
  hAxesDataOverFitTPC->SetStats(kFALSE);
  hAxesDataOverFitTPC->SetTitle(";nsigma for 0.4 GeV/c < p < 0.5 GeV/c;Data/Fit (TPC)");
  hAxesDataOverFitTPC->DrawCopy();
  TLegend * lDataOverFitTPC = new TLegend(.59,.69,.99,.99);
  lDataOverFitTPC->SetFillColor(0);
  lDataOverFitTPC->AddEntry(&gDataOverFitTPC[0],Form("Run %i",runs[irun]),"");
  if (useMC) lDataOverFitTPC->AddEntry(&gDataOverFitTPC[0],"MC","");
  else lDataOverFitTPC->AddEntry(&gDataOverFitTPC[0],"DATA","");
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      gDataOverFitTPC[ipart].SetMarkerStyle(Marker[ipart]);
      gDataOverFitTPC[ipart].SetMarkerColor(Color[ipart]);
      gDataOverFitTPC[ipart].DrawClone("Psame");
      lDataOverFitTPC->AddEntry(&gDataOverFitTPC[ipart],Particle[ipart].Data(),"p");
    }
  lDataOverFitTPC->DrawClone();

  //------------------------------------------------------
  //                        TOF                          -
  //------------------------------------------------------

  // define the 2D nsigma histograms, the projections, and the fitting functions
  // (one each for pions, kaons, and protons)
  TH2F * nsig_TOF[nPart];
  TH1F * nsig_TOF_proj[nPart];
  TF1 * fitFuncTOF[nPart];
  for (Int_t ipart; ipart<nPart; ipart++)
    {
      nsig_TOF[ipart] = new TH2F;
      nsig_TOF_proj[ipart] = new TH1F;
    }
  
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      // create projection
      nsig_TOF[ipart] = (TH2F*)((TH2F*)hman->GetNSigHistogram(Form("hHistNSig%sTOF",Particle[ipart].Data())))->Clone();
      if (useMC)
	nsig_TOF_proj[ipart] = (TH1F*)nsig_TOF[ipart]->ProjectionY(Form("TOF NsigProjection %s, mc [%.1f,%.1f]",Particle[ipart].Data(),ProjTOF[0],ProjTOF[1]), nsig_TOF[ipart]->GetXaxis()->FindBin(ProjTOF[0]), nsig_TOF[ipart]->GetXaxis()->FindBin(ProjTOF[1]));
      else
	nsig_TOF_proj[ipart] = (TH1F*)nsig_TOF[ipart]->ProjectionY(Form("TOF NsigProjection %s, data [%.1f,%.1f]",Particle[ipart].Data(),ProjTOF[0],ProjTOF[1]), nsig_TOF[ipart]->GetXaxis()->FindBin(ProjTOF[0]), nsig_TOF[ipart]->GetXaxis()->FindBin(ProjTOF[1]));
	
      // fit the peak of interest with a gaussian
      fitFuncTOF[ipart] = new TF1("fitFuncTOF","gaus",fitRangeTOF[ipart][0],fitRangeTOF[ipart][1]);
      TFitResultPtr resultTOF = nsig_TOF_proj[ipart]->Fit(fitFuncTOF[ipart],"NRSLQ");

      // draw the projections and fits
      cNSigProjFits->cd(ipart+1+nPart);
      gPad->SetLogy();

      nsig_TOF_proj[ipart]->GetXaxis()->SetTitle(Form("TOF nsigma for %.1f GeV/c < p < %.1f GeV/c",ProjTOF[0],ProjTOF[1]));
      nsig_TOF_proj[ipart]->SetStats(kFALSE);
      nsig_TOF_proj[ipart]->DrawCopy();
      fitFuncTOF[ipart]->DrawCopy("same");
      TLegend * lNSigTOF = new TLegend(0.59,0.59,0.99,0.99);
      lNSigTOF->SetFillColor(0);
      lNSigTOF->AddEntry(nsig_TOF_proj[ipart],Form("%ss, TOF",Particle[ipart].Data()),"");
      lNSigTOF->AddEntry(nsig_TOF_proj[ipart],Form("Run %i",runs[irun]),"");
      if (useMC) lNSigTOF->AddEntry(nsig_TOF_proj[ipart],"MC","");
      else lNSigTOF->AddEntry(nsig_TOF_proj[ipart],"DATA","");
      //lNSigTOF->AddEntry(nsig_TOF_proj[ipart],Form("#chi^{2}/nDOF = %.2f",fitFuncTOF[ipart]->GetChisquare()/fitFuncTOF[ipart]->GetNDF()),"");
      lNSigTOF->DrawClone();

      // fill gDataOverFitTOF[] graphs
      Int_t nGraphPoints = 0;
      for (Float_t iPt = fitRangeTOF[ipart][0]; iPt <= fitRangeTOF[ipart][1]; iPt += (fitRangeTOF[ipart][1] - fitRangeTOF[ipart][0])/15)
	{
	  gDataOverFitTOF[ipart].SetPoint(nGraphPoints,iPt, nsig_TOF_proj[ipart]->GetBinContent(nsig_TOF_proj[ipart]->FindBin(iPt)) / fitFuncTOF[ipart]->Eval(iPt));
	  nGraphPoints++;
	}

    
      // fill the histograms with the fit parameters and run numbers
      switch (ipart)
	{
	case 0: // pion
	  TOFnsigMeanTrendPion->SetBinContent(irun+1, resultTOF->Parameter(1));
	  TOFnsigMeanTrendPion->SetBinError(irun+1, resultTOF->ParError(1));
	  TOFnsigSigmaTrendPion->SetBinContent(irun+1, resultTOF->Parameter(2));
	  TOFnsigSigmaTrendPion->SetBinError(irun+1, resultTOF->ParError(2));
	  break;
	case 1: // kaon
	  TOFnsigMeanTrendKaon->SetBinContent(irun+1, resultTOF->Parameter(1));
	  TOFnsigMeanTrendKaon->SetBinError(irun+1, resultTOF->ParError(1));
	  TOFnsigSigmaTrendKaon->SetBinContent(irun+1, resultTOF->Parameter(2));
	  TOFnsigSigmaTrendKaon->SetBinError(irun+1, resultTOF->ParError(2));
	  break;
	case 2: // proton
	  TOFnsigMeanTrendProton->SetBinContent(irun+1, resultTOF->Parameter(1));
	  TOFnsigMeanTrendProton->SetBinError(irun+1, resultTOF->ParError(1));
	  TOFnsigSigmaTrendProton->SetBinContent(irun+1, resultTOF->Parameter(2));
	  TOFnsigSigmaTrendProton->SetBinError(irun+1, resultTOF->ParError(2));
	  break;
	default:
	  Printf("\n!!! ERROR in TOF switch control structure in CheckNSigmaStability.C !!!\n");
	  break;
	}
    } // end loop over ipart

  // draw the data/fit
  cDataOverFit->cd(2);
  TH2F * hAxesDataOverFitTOF = new TH2F("hAxesDataOverFitTOF","",100,-5,5,100,-3,5);
  hAxesDataOverFitTOF->SetStats(kFALSE);
  hAxesDataOverFitTOF->SetTitle(";nsigma for 0.9 GeV/c < p < 1.0 GeV/c;Data/Fit (TOF)");
  hAxesDataOverFitTOF->DrawCopy();
  TLegend * lDataOverFitTOF = new TLegend(.59,.69,.99,.99);
  lDataOverFitTOF->SetFillColor(0);
  lDataOverFitTOF->AddEntry(&gDataOverFitTOF[0],Form("Run %i",runs[irun]),"");
  if (useMC) lDataOverFitTOF->AddEntry(&gDataOverFitTOF[0],"MC","");
  else lDataOverFitTOF->AddEntry(&gDataOverFitTOF[0],"DATA","");
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      gDataOverFitTOF[ipart].SetMarkerStyle(Marker[ipart]);
      gDataOverFitTOF[ipart].SetMarkerColor(Color[ipart]);
      gDataOverFitTOF[ipart].DrawClone("Psame");
      lDataOverFitTOF->AddEntry(&gDataOverFitTOF[ipart],Particle[ipart].Data(),"p");
    }
  lDataOverFitTOF->DrawClone();




  // save the projections and fits to a pdf file once per run
  if (useMC)
    {
      if (irun == 0) cNSigProjFits->Print("Plots/MC/AODnSigmaProjFits.pdf(","pdf");
      else if (irun < nRuns-1) cNSigProjFits->Print("Plots/MC/AODnSigmaProjFits.pdf","pdf");
      else if (irun == nRuns-1) cNSigProjFits->Print("Plots/MC/AODnSigmaProjFits.pdf)","pdf");
    } 
  else
    {
      if (irun == 0) cNSigProjFits->Print("Plots/DATA/AODnSigmaProjFits.pdf(","pdf");
      else if (irun < nRuns-1) cNSigProjFits->Print("Plots/DATA/AODnSigmaProjFits.pdf","pdf");
      else if (irun == nRuns-1) cNSigProjFits->Print("Plots/DATA/AODnSigmaProjFits.pdf)","pdf");
    }
  cNSigProjFits->Close();



  // also save the data/fit plots to a pdf file once per run
  if (useMC)
    {
      if (irun == 0) cDataOverFit->Print("Plots/MC/NSigDataOverFit.pdf(","pdf");
      else if (irun < nRuns-1) cDataOverFit->Print("Plots/MC/NSigDataOverFit.pdf","pdf");
      else if (irun == nRuns-1) cDataOverFit->Print("Plots/MC/NSigDataOverFit.pdf)","pdf");
    }
  else
    {
      if (irun == 0) cDataOverFit->Print("Plots/DATA/NSigDataOverFit.pdf(","pdf");
      else if (irun < nRuns-1) cDataOverFit->Print("Plots/DATA/NSigDataOverFit.pdf","pdf");
      else if (irun == nRuns-1) cDataOverFit->Print("Plots/DATA/NSigDataOverFit.pdf)","pdf");
    }
  cDataOverFit->Close();
}




