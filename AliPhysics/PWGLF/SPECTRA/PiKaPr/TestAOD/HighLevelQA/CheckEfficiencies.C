/////////////////////////////////////////////////////////////////////
// CheckEfficiencies.C (called by AODQAChecks.C)                   //
//                                                                 //
// Written by John Groh                                            //
/////////////////////////////////////////////////////////////////////

void CheckEfficiencies(AliSpectraAODHistoManager * hman,
		       TCanvas*& cEfficienciesAllRuns,
		       TH1F*& EfficiencyPiPlus,
		       TH1F*& EfficiencyKPlus,
		       TH1F*& EfficiencyProton,
		       TH1F*& EfficiencyPiMinus,
		       TH1F*& EfficiencyKMinus,
		       TH1F*& EfficiencyAntiproton,
		       Float_t FixedPtEff,
		       Int_t runs[],
		       Int_t nRuns,
		       Int_t irun)
{
  // canvas for printing individual runs seperately to Efficiencies.pdf
  TCanvas * cEfficienciesIndiv = new TCanvas("cEfficienciesIndiv","cEfficienciesIndiv");
  cEfficienciesIndiv->Divide(3,2);

  // calculate correction factors
  TH1F * CorrFact[nCharge*nPart];
  TLegend * lEfficiencyAllRuns[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  // get recreated histo and format it
	  Int_t index = ipart + nPart*icharge;
	  TString hname = Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  CorrFact[index] = (TH1F*)((TH1F*)hman->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
	  CorrFact[index]->SetName(Form("CorrFact_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
	  CorrFact[index]->SetTitle(Form("CorrFact_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
	  CorrFact[index]->SetMarkerStyle(Marker[index]);
	  CorrFact[index]->SetMarkerColor(Color[ipart]);
	  CorrFact[index]->SetLineColor(Color[ipart]);
	  CorrFact[index]->SetStats(kFALSE);
	  CorrFact[index]->GetYaxis()->SetTitle("MC Correction Factor");

	  // divide it by MC truth histo
	  hname = Form("hHistPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  CorrFact[index]->Divide(CorrFact[index],(TH1F*)((TH1F*)hman->GetPtHistogram1D(hname.Data(),-1,-1))->Clone(),1,1,"B"); // binomial error
	  
	  // draw efficiency histos for individual runs
	  cEfficienciesIndiv->cd(index+1);
	  TLegend * lEfficiencyIndiv = new TLegend(.55,.15,.85,.25);
	  lEfficiencyIndiv->AddEntry(CorrFact[index],Form("%s, Run %i",Names[index].Data(), runs[irun]),"lpe");
	  CorrFact[index]->DrawCopy();
	  lEfficiencyIndiv->DrawClone();

	  // superimpose efficiency histos on a different canvas for all runs
	  cEfficienciesAllRuns->cd(index+1);
	  if (irun == 0)
	    {
	      CorrFact[index]->DrawCopy(); 
	      lEfficiencyAllRuns[index] = new TLegend(.15,.75,.35,.85);
	      lEfficiencyAllRuns[index]->AddEntry(CorrFact[index],Names[index].Data(),"lpe");
	      lEfficiencyAllRuns[index]->SetFillColor(0);
	      lEfficiencyAllRuns[index]->DrawClone();
	    }
	  else CorrFact[index]->DrawCopy("same");

	  // using the values of FixedPtEff for a fixed pt, plot the correction factors vs the run #
	  switch (index)
	    {
	    case 0: // PiPlus
	      EfficiencyPiPlus->SetBinContent(irun+1,CorrFact[index]->GetBinContent(CorrFact[index]->FindBin(FixedPtEff)));
	      EfficiencyPiPlus->SetBinError(irun+1,CorrFact[index]->GetBinError(CorrFact[index]->FindBin(FixedPtEff)));
	      break;
	    case 1: // KPlus
	      EfficiencyKPlus->SetBinContent(irun+1,CorrFact[index]->GetBinContent(CorrFact[index]->FindBin(FixedPtEff)));
	      EfficiencyKPlus->SetBinError(irun+1,CorrFact[index]->GetBinError(CorrFact[index]->FindBin(FixedPtEff)));
	      break;
	    case 2: // Proton
	      EfficiencyProton->SetBinContent(irun+1,CorrFact[index]->GetBinContent(CorrFact[index]->FindBin(FixedPtEff)));
	      EfficiencyProton->SetBinError(irun+1,CorrFact[index]->GetBinError(CorrFact[index]->FindBin(FixedPtEff)));
	      break;
	    case 3: // PiMinus
	      EfficiencyPiMinus->SetBinContent(irun+1,CorrFact[index]->GetBinContent(CorrFact[index]->FindBin(FixedPtEff)));
	      EfficiencyPiMinus->SetBinError(irun+1,CorrFact[index]->GetBinError(CorrFact[index]->FindBin(FixedPtEff)));
	      break;
	    case 4: // KMinus
	      EfficiencyKMinus->SetBinContent(irun+1,CorrFact[index]->GetBinContent(CorrFact[index]->FindBin(FixedPtEff)));
	      EfficiencyKMinus->SetBinError(irun+1,CorrFact[index]->GetBinError(CorrFact[index]->FindBin(FixedPtEff)));
	      break;
	    case 5: // Antiproton
	      EfficiencyAntiproton->SetBinContent(irun+1,CorrFact[index]->GetBinContent(CorrFact[index]->FindBin(FixedPtEff)));
	      EfficiencyAntiproton->SetBinError(irun+1,CorrFact[index]->GetBinError(CorrFact[index]->FindBin(FixedPtEff)));
	      break;
	    default:
	      Printf("\n!!! ERROR in switch statement in CheckEfficiencies.C !!!\n");
	      break;
	    }
	} // end loop over ipart
    } // end loop over icharge
  
  // save the projections and fits to a pdf file once per run
  if (irun == 0) cEfficienciesIndiv->SaveAs("Plots/MC/AODEfficiencies.pdf(","pdf");
  else if (irun < nRuns-1) cEfficienciesIndiv->SaveAs("Plots/MC/AODEfficiencies.pdf","pdf");
  else if (irun == nRuns-1) cEfficienciesIndiv->SaveAs("Plots/MC/AODEfficiencies.pdf)","pdf");
   cEfficienciesIndiv->Close();

}










