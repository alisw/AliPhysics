/*
 *  CalcEffic.C
 */

TEfficiency* CalcEfficHadrons(const char* fileName, Bool_t makeDraw=kFALSE)
{
	// open the input file
	TFile *file = new TFile(fileName);
	
	// ********************************************
	// parameters to check before running the macro
	// ********************************************

	const Int_t NHISTS = 4; // Check the number of histograms for different particle species
	const Int_t NOUTPUTS = 3;
	const Int_t NHISTOUT[NOUTPUTS] = {1,1,1};
	const Int_t IHISTOUT[NOUTPUTS][NHISTS] = {{0,-1,-1,-1},{1,-1,-1,-1},{2,-1,-1,-1}};
	const Float_t CUT_RES = 0.02;
	
	Int_t style[NOUTPUTS] = {20,21,22};
	Int_t color[NOUTPUTS] = {1,2,4};
		
	const Int_t fgNumOfPtBins = 111; // Check the number of eta bins in the histograms
	const Int_t fgNumOfEtaBins = 16; // Check the number of E bins in the histograms
	const Int_t fgNumOfRBins = 45;
	
	Double_t fgPtAxis[117]= {0.0,0.01,0.02,0.03,0.04, 0.05, 0.06,0.07,0.08,0.09, 0.10,0.11, .12,0.13, .14,0.15, .16,0.17, .18,0.19,
		0.2, .22, .24, .26, .28, 0.30, 0.32, .34, .36, .38, 0.40, .42, .44, .46, .48,
		0.5, .52, .54, .56, .58, 0.60, 0.62, .64, .66, .68, 0.70, .72, .74, .76, .78,
		.80, .82, .84, .86, .88, 0.90, 0.92, .94, .96, .98, 1.00,1.05, 1.1,1.15, 1.2,
		1.25, 1.3,1.35,1.40,1.45, 1.50, 1.55, 1.6,1.65, 1.7, 1.75, 1.8,1.85, 1.9,1.95,
		2.0, 2.2, 2.4, 2.6, 2.8, 3.00, 3.20, 3.4, 3.6, 3.8, 4.00, 4.2, 4.4, 4.6, 4.8,
		5.0, 5.5, 6.0, 6.5, 7.0, 7.50, 8.00, 8.5, 9.0, 9.5, 10.0,12.0,14.0,16.0,18.0,
		20.0,25.0,30.0,35.0,40.0, 45.0, 50.0}; 
	
	// declare histograms and graphs
	TH2F *histNum[NHISTS];
	TH2F *histDen[NHISTS];
	TGraphErrors *graph[NOUTPUTS];
	TH1D* projYNum;
	TEfficiency *effic[NOUTPUTS];
	char efficName[50];
	
	// retrieve the input list of histogram. Check the TList name in the input file.
	TList *list = (TList*) file->Get("out1");
	
	// retrieve the histograms in the list. Check the name of the histograms
	histNum[0] = (TH2F*)list->FindObject("fHistPionRec_ResPt_EmcalMC");
	histNum[1] = (TH2F*)list->FindObject("fHistKaonRec_ResPt_EmcalMC");
	histNum[2] = (TH2F*)list->FindObject("fHistProtonRec_ResPt_EmcalMC");
	histNum[3] = (TH2F*)list->FindObject("fHistMuonRec_ResPt_EmcalMC");
	
	// retrieve the histograms in the list. Check the name of the histograms
	histDen[0] = (TH2F*)list->FindObject("fHistPionAcc_EtaPt_EmcalMC");
	histDen[1] = (TH2F*)list->FindObject("fHistKaonAcc_EtaPt_EmcalMC");
	histDen[2] = (TH2F*)list->FindObject("fHistProtonAcc_EtaPt_EmcalMC");
	histDen[3] = (TH2F*)list->FindObject("fHistMuonAcc_EtaPt_EmcalMC");
	
	// ********************************************

	Float_t x[fgNumOfPtBins]={0}, ex[fgNumOfPtBins]={0};
	Float_t y[fgNumOfPtBins]={0}, ey[fgNumOfPtBins]={0};
	Float_t num=0, den=0;
	//Int_t num=0, den=0;
	Float_t Res=0;
	
	// loop over different desired outputs
	for (int iOut=0; iOut<NOUTPUTS; iOut++)
	{
		sprintf(efficName,"effic_%d",iOut);
		effic[iOut] = new TEfficiency(efficName,efficName,fgNumOfPtBins,fgPtAxis);

		// loop over E bins
		for (int ix=0; ix<fgNumOfPtBins; ix++)
		{
			// initialize ET variables for a new particle species
			x[ix]=histNum[0]->GetXaxis()->GetBinCenter(ix+1);
			y[ix]=0;
			ex[ix]=0;
			ey[ix]=0;
			num = 0;
			den = 0;
			
			// loop over eta bins
			for (int iy=0; iy<fgNumOfEtaBins; iy++)
			{
				for (int iHist=0; iHist<NHISTOUT[iOut]; iHist++)
				{
					den += histDen[IHISTOUT[iOut][iHist]]->GetBinContent(ix+1,iy+1); // sum over all E bins in order to get total ET
				}
			}
			
			// loop over residual bins
			for (int iHist=0; iHist<NHISTOUT[iOut]; iHist++)
			{
				projYNum = histNum[IHISTOUT[iOut][iHist]]->ProjectionY();
				for (int iy=0; iy<fgNumOfRBins; iy++)
				{
					Res = projYNum->GetBinCenter(iy+1);
					if (Res<CUT_RES)
						num += histNum[IHISTOUT[iOut][iHist]]->GetBinContent(ix+1,iy+1); // sum over all E bins in order to get total ET
				}
			}
			
			if ((num>0) && (den>0))
			{
				effic[iOut]->SetTotalEvents(ix,den);
				effic[iOut]->SetPassedEvents(ix,num);
				y[ix] = num/den;
				ey[ix] = y[ix]*sqrt(1/num+1/den);
				//ey[ix] = ((num+1)*(num+2))/((den+2)*(den+3))-((num+1)*(num+1))/((den+2)*(den+2));
			}
			else
			{
				y[ix] = 0;	
				ey[ix] = 0;
			}
			
		} // end of loop over E bins

		graph[iOut] = new TGraphErrors(fgNumOfPtBins,x,y,ex,ey); // graphic of ET(>E_cut)/ET(total) for a given particle species and E cut

	} // end of loop over different outputs

	
	// Draw the plot
	
	if (makeDraw)
	{
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
		
		/*
		 for (int i=0; i<NOUTPUTS; i++)
		 {
		 graph[i]->SetMarkerStyle(style[i]);
		 graph[i]->SetMarkerColor(color[i]);
		 graph[i]->SetLineColor(color[i]);
		 graph[i]->SetFillColor(0);
		 if (i == 0) 
		 {
		 graph[i]->GetXaxis()->SetTitle("E (GeV)");
		 graph[i]->GetYaxis()->SetTitle("effic");
		 graph[i]->SetMaximum(1.0);
		 graph[i]->SetMinimum(0.0);
		 graph[i]->Draw("AP");
		 }
		 else 
		 graph[i]->Draw("P");
		 }
		 */
		for (int i=0; i<NOUTPUTS; i++)
		{
			effic[i]->SetMarkerStyle(style[i]);
			effic[i]->SetMarkerColor(color[i]);
			effic[i]->SetLineColor(color[i]);
			effic[i]->SetFillColor(0);
			effic[i]->SetTitle("efficiency; p_{T} (GeV/c); #epsilon");
			if (i == 0) 
			{
				effic[i]->Draw();
			}
			else 
				effic[i]->Draw("Psame");
		}
		
		TLegend *leg = new TLegend(0.65,0.2,0.95,0.5);
		leg->AddEntry(effic[0],"pions");
		leg->AddEntry(effic[1],"kaons");
		leg->AddEntry(effic[2],"protons");
		//leg->AddEntry(effic[3],"muons");
		leg->SetFillStyle(0);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->Draw();
	}
	
	return effic[0];
}
