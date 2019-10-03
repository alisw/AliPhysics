/*
 *  CalcET.C
 */

void CalcET(const char* fileName, const char* fileNameMC)
{
	// open the input file
	TFile *file = new TFile(fileName);
	
	// ********************************************
	// parameters to check before running the macro
	// ********************************************
	
	const Int_t NCUTS = 3; // Choose the number of cuts to check
	Float_t eCut[NCUTS] = {0.25,0.5,0.75}; // Choose the value of the cuts
	Float_t ptCut[NCUTS] = {0.25,0.5,0.75}; // Choose the value of the cuts
	Double_t matchCut = 0.02;
	
	const Int_t fgNumOfEBins = 78; // Check the number of eta bins in the histograms
	const Int_t fgNumOfEtaBins = 16; // Check the number of E bins in the histograms
	const Int_t fgNumOfPtBins = 111; // Check the number of pT bins in the histograms
	
	// declare histograms and graphs
	TGraphErrors *graph[NCUTS];
	
	// retrieve the input list of histogram. Check the TList name in the input file.
	TList *list = (TList*) file->Get("out1");
	
	// retrieve the histograms in the list. Check the name of the histograms
	THnSparseD* histAll = (THnSparseD*)list->FindObject("fHistAllRec_ETDep_EmcalRec");

	THnSparseD* hnSparse = (THnSparseD*)list->FindObject("fHistMuon_ETDep_EmcalRec");
	hnSparse->GetAxis(6)->SetRangeUser(0.,matchCut);
	TH3D* histMuon = hnSparse->Projection(2,0,1);
	
	hnSparse = (THnSparseD*)list->FindObject("fHistMuon_ETDep_EmcalRec");
	hnSparse->GetAxis(6)->SetRangeUser(0.,matchCut);
	
	hnSparse = (THnSparseD*)list->FindObject("fHistMuon_ETDep_EmcalRec");
	hnSparse->GetAxis(6)->SetRangeUser(0.,matchCut);
	TH3D* histPion = ((THnSparseD*)list->FindObject("fHistPion_ETDep_EmcalRec"))->Projection(2,0,1);
	
	hnSparse = (THnSparseD*)list->FindObject("fHistMuon_ETDep_EmcalRec");
	hnSparse->GetAxis(6)->SetRangeUser(0.,matchCut);
	TH3D* histKaon = ((THnSparseD*)list->FindObject("fHistKaon_ETDep_EmcalRec"))->Projection(2,0,1);
	
	hnSparse = (THnSparseD*)list->FindObject("fHistMuon_ETDep_EmcalRec");
	hnSparse->GetAxis(6)->SetRangeUser(0.,matchCut);
	TH3D* histProton = ((THnSparseD*)list->FindObject("fHistProton_ETDep_EmcalRec"))->Projection(2,0,1);
	
	TH3D* histCharged = new TH3D(*histMuon);
	histCharged->Add(histPion);
	histCharged->Add(histKaon);
	histCharged->Add(histProton);
	
	TEfficiency *efficClu = CalcEffic(fileNameMC);
	TEfficiency *efficChHad = CalcEfficHadrons(fileNameMC);
	Double_t f_acc = CalcCorrAcc(fileNameMC);
	Double_t f_Ecut = CalcCorrEcut(fileNameMC);
	Double_t f_neutral = CalcCorrNeutral(fileNameMC);
	
	// ********************************************
	
	Float_t x[fgNumOfEtaBins]={0}, ex[fgNumOfEtaBins]={0};
	Float_t y[fgNumOfEtaBins]={0}, ey[fgNumOfEtaBins]={0};
	
	// loop over different E cuts
	for (int iCut=0; iCut<NCUTS; iCut++)
	{
		// loop over eta bins
		for (int iy=0; iy<fgNumOfEtaBins; iy++)
		{
			// initialize ET variables for a new particle species
			Float_t E=0, pt=0, ET=0, ET_All=0, ET_Had=0, errorSqET_All=0, errorSqET_Had=0;
			x[iy]=0;
			y[iy]=0;
			ex[iy]=0;
			ey[iy]=0;
			
			// loop over E bins
			for (int ix=0; ix<fgNumOfEBins; ix++)
			{
				E = histAll->GetXaxis()->GetBinCenter(ix+1);
				if (E > eCut[iCut])
				{
					ET = histAll->GetBinContent(ix+1,iy+1); // sum over all E bins in order to get total ET
					errorSqET_All += pow(E_err(E)*(ET/E),2);
					ET_All += ET/efficClu->GetEfficiency(ix);
										
				}				
				
				// loop over pt bins
				for (int iz=0; iz<fgNumOfPtBins; iz++)
				{
					pt = histCharged->GetZaxis()->GetBinCenter(iz+1);
					if (pt > ptCut[iCut])
					{
						ET = histCharged->GetBinContent(ix+1,iy+1,iz+1); // sum over all pt bins in order to get total ET
						errorSqET_Had += pow(E_err(E)*(ET/E),2);
						ET_Had += ET/efficChHad->GetEfficiency(iz);
					}				
				} // end of loop over pt bins
				
			} // end of loop over E bins

			x[iy] = histChHad->GetYaxis()->GetBinCenter(iy+1); // x coordinate of eta bin
			
			y[iy] = (ET_All-(ET_Had/f_neutral))/(f_acc*f_Ecut); // y coordinate of eta bin (for a given particle species)
			
			ey[iy]=(sqrt(errorSqET_AllHad + errorSqET_Had/pow(f_neutral,2))/(f_acc*f_Ecut);
			
		} // end of loop over eta bins
		
		graph[iCut] = new TGraphErrors(fgNumOfEtaBins,x,y,ex,ey); // graphic of ET(>E_cut)/ET(total) for a given E cut (all particle species combined)
	} // end of loop over different E cuts
	
	
	// Draw the plot
	
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);	
	
	TCanvas *c = new TCanvas("c","c",500,400);
	//c->SetTopMargin(0.04);
	//c->SetRightMargin(0.04);
	//c->SetLeftMargin(0.181452);
	//c->SetBottomMargin(0.134409);
	c->SetBorderSize(0);
	c->SetFillColor(0);
	c->SetBorderMode(0);
	c->SetFrameFillColor(0);
	c->SetFrameBorderMode(0);
	
	Int_t style[NCUTS] = {20,21,22};
	Int_t color[NCUTS] = {2,3,4};
	
	for (int i=0; i<NCUTS; i++)
	{
		graph[i]->SetMarkerStyle(style[i]);
		graph[i]->SetMarkerColor(color[i]);
		graph[i]->SetLineColor(color[i]);
		graph[i]->SetFillColor(0);
		if (i == 0) 
		{
			graph[i]->GetXaxis()->SetTitle("#eta");
			graph[i]->GetYaxis()->SetTitle("E_{T}^{Neutral+Charged Hadrons}/E_{T}^{Charged hadrons}");
			//graph[i]->SetMaximum(1.0);
			graph[i]->SetMinimum(0.0);
			graph[i]->Draw("AP");
		}
		else
			graph[i]->Draw("P");
	}
	
	TLegend *leg = new TLegend(0.65,0.2,0.95,0.5);
	leg->AddEntry(graph[0],"E>250 MeV");
	leg->AddEntry(graph[1],"E>500 MeV");
	leg->AddEntry(graph[2],"E>750 MeV");
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->Draw();	
}

Float_t E_err(Float_t E)
{
	return (E*sqrt(pow(4.35/E,2)+pow(9.07,2)/E+pow(1.63,2)))/100;	
}
