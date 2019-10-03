/*
 *  CalcCorrAcc.C
 */

void CalcCorrAcc(const char* fileName, Bool_t makeDraw=kFALSE, Float_t* yy=0, Float_t* eyy=0)
{
	// open the input file
	TFile *file = new TFile(fileName);
	
	// ********************************************
	// parameters to check before running the macro
	// ********************************************
	
	const Int_t NCUTS = 3; // Choose the number of cuts to check
	Float_t eCut[NCUTS] = {0.1,0.25,0.5}; // Choose the value of the cuts
	//Float_t eCut[NCUTS] = {0.5}; // Choose the value of the cuts
	
	Int_t style[NCUTS] = {20,21,22};
	Int_t color[NCUTS] = {1,2,4};
	//Int_t style[NCUTS] = {20};
	//Int_t color[NCUTS] = {1};	
	
	const Int_t fgNumOfEBins = 78; // Check the number of eta bins in the histograms
	const Int_t fgNumOfEtaBins = 16; // Check the number of E bins in the histograms
	
	// declare histograms and graphs
	TH2F *hist;
	TH2F *histAllEM;
	TH2F *histAccEM;
	
	TGraphErrors *graph[NCUTS];
	
	// retrieve the input list of histogram. Check the TList name in the input file.
	TList *list = (TList*) file->Get("out1");
	
	// retrieve the histograms in the list. Check the name of the histograms
	hist = (TH2F*)list->FindObject("fHistElectron_EtaE_ET_EmcalMC");
	histAllEM = new TH2F(*hist);	
	hist = (TH2F*)list->FindObject("fHistConvElectron_EtaE_ET_EmcalMC");
	histAllEM->Add(hist);
	hist = (TH2F*)list->FindObject("fHistScatElectron_EtaE_ET_EmcalMC");
	histAllEM->Add(hist);
	hist = (TH2F*)list->FindObject("fHistGamma_EtaE_ET_EmcalMC");
	histAllEM->Add(hist);
	hist = (TH2F*)list->FindObject("fHistAnnihGamma_EtaE_ET_EmcalMC");
	histAllEM->Add(hist);
	hist = (TH2F*)list->FindObject("fHistScatGamma_EtaE_ET_EmcalMC");
	histAllEM->Add(hist);
	
	hist = (TH2F*)list->FindObject("fHistElectronAcc_EtaE_ET_EmcalMC");
	histAccEM = new TH2F(*hist);	
	hist = (TH2F*)list->FindObject("fHistConvElectronAcc_EtaE_ET_EmcalMC");
	histAccEM->Add(hist);
	hist = (TH2F*)list->FindObject("fHistScatElectronAcc_EtaE_ET_EmcalMC");
	histAccEM->Add(hist);
	hist = (TH2F*)list->FindObject("fHistGammaAcc_EtaE_ET_EmcalMC");
	histAccEM->Add(hist);
	hist = (TH2F*)list->FindObject("fHistScatGammaAcc_EtaE_ET_EmcalMC");
	histAccEM->Add(hist);
	hist = (TH2F*)list->FindObject("fHistAnnihGammaAcc_EtaE_ET_EmcalMC");
	histAccEM->Add(hist);
	
	
	// ********************************************
	
	Float_t x[fgNumOfEtaBins]={0}, ex[fgNumOfEtaBins]={0};
	Float_t y[fgNumOfEtaBins]={0}, ey[fgNumOfEtaBins]={0};
	
	// returning the output that is in memory
	yy=y;
	eyy = ey;
	
	// loop over different E cuts
	for (int iCut=0; iCut<NCUTS; iCut++)
	{
		// loop over eta bins
		for (int iy=0; iy<fgNumOfEtaBins; iy++)
		{
			// initialize ET variables for a new particle species
			Float_t E=0, ET=0, ET_AllEM=0, ET_AccEM=0, errorSqET_AllEM=0, errorSqET_AccEM=0;
			x[iy]=0;
			y[iy]=0;
			ex[iy]=0;
			ey[iy]=0;
			
			// loop over E bins
			for (int ix=0; ix<fgNumOfEBins; ix++)
			{
				E = histAccEM->GetXaxis()->GetBinLowEdge(ix+1);
				
				if (E >= eCut[iCut])
				{
					ET = histAccEM->GetBinContent(ix+1,iy+1); // sum over all E bins in order to get total ET
					errorSqET_AccEM += pow(E_err(E)*(ET/E),2);
					ET_AccEM += ET;
				}				
					
				ET = histAllEM->GetBinContent(ix+1,iy+1); // sum over all E bins in order to get total ET
				if (E > 0)
					errorSqET_AllEM += pow(E_err(E)*(ET/E),2);
				ET_AllEM += ET;
			} // end of loop over E bins
			
			x[iy] = histAccEM->GetYaxis()->GetBinCenter(iy+1); // x coordinate of eta bin
			
			if ( (ET_AllEM > 0) && (ET_AccEM > 0) )
			{
				y[iy] = ET_AccEM/ET_AllEM; // y coordinate of eta bin (for a given particle species)
				ey[iy]=y[iy]*sqrt(errorSqET_AllEM/pow(ET_AllEM,2) + errorSqET_AccEM/pow(ET_AccEM,2));
			}
			else
			{
				y[iy] = 0;
				ey[iy] = 0;
			}
		} // end of loop over eta bins
		
		graph[iCut] = new TGraphErrors(fgNumOfEtaBins,x,y,ex,ey); // graphic of ET(>E_cut)/ET(total) for a given E cut (all particle species combined)
	} // end of loop over different E cuts
	
	
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
		
		for (int i=0; i<NCUTS; i++)
		{
			graph[i]->SetMarkerStyle(style[i]);
			graph[i]->SetMarkerColor(color[i]);
			graph[i]->SetLineColor(color[i]);
			graph[i]->SetFillColor(0);
			if (i == 0) 
			{
				graph[i]->GetXaxis()->SetTitle("#eta");
				graph[i]->GetYaxis()->SetTitle("E_{T}^{EM Acceptance}/E_{T}^{EM All}");
				graph[i]->SetMaximum(0.1);
				graph[i]->SetMinimum(0.0);
				graph[i]->Draw("AP");
			}
			else
				graph[i]->Draw("P");
		}
		
		TLegend *leg = new TLegend(0.65,0.2,0.95,0.5);
		leg->AddEntry(graph[0],"E>100 MeV");
		leg->AddEntry(graph[1],"E>250 MeV");
		leg->AddEntry(graph[2],"E>500 MeV");
		leg->SetFillStyle(0);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->Draw();		
	}
}

Float_t E_err(Float_t E)
{
	return (E*sqrt(pow(4.35/E,2)+pow(9.07,2)/E+pow(1.63,2)))/100;	
}
