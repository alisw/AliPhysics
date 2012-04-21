/*
 *  CalcEffic.C
 */

TEfficiency* CalcEffic(const char* fileName, Bool_t makeDraw=kFALSE)
{
	// open the input file
	TFile *file = new TFile(fileName);
	
	// ********************************************
	// parameters to check before running the macro
	// ********************************************

	const Float_t E_MIN = 0.1;
	const Int_t NHISTS = 6; // Check the number of histograms for different particle species
	const Int_t NOUTPUTS = 3;
	const Int_t NHISTOUT[NOUTPUTS] = {6,3,3};
	Int_t IHISTOUT[NOUTPUTS][NHISTS] = {{0,1,2,3,4,5},{0,1,2,-1,-1,-1},{3,4,5,-1,-1,-1}};
	
	Int_t style[NOUTPUTS] = {20,21,22};
	Int_t color[NOUTPUTS] = {1,2,4};
		
	const Int_t fgNumOfEBins = 78; // Check the number of eta bins in the histograms
	const Int_t fgNumOfEtaBins = 16; // Check the number of E bins in the histograms
	Double_t fgEAxis[79]={0., 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
		1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,11.,
		12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,
		32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.};
	
	// declare histograms and graphs
	TH2F *histNum[NHISTS];
	TH2F *histDen[NHISTS];
	TGraphErrors *graph[NOUTPUTS];
	TGraphErrors *graphNum[NOUTPUTS];
	TGraphErrors *graphDen[NOUTPUTS];
	TEfficiency *effic[NOUTPUTS];
	char efficName[50];
	
	//define canvas
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);	
	
	TCanvas *c1;
	TCanvas *c2;
	
	if (makeDraw)
	{
		c1 = new TCanvas("c1","c1",500,400);
		//c1->SetTopMargin(0.04);
		//c1->SetRightMargin(0.04);
		//c1->SetLeftMargin(0.181452);
		//c1->SetBottomMargin(0.134409);
		c1->SetBorderSize(0);
		c1->SetFillColor(0);
		c1->SetBorderMode(0);
		c1->SetFrameFillColor(0);
		c1->SetFrameBorderMode(0);
		
		c2 = new TCanvas("c2","c2",500,400);
		c2->Divide(1,3);
	}
	
	// retrieve the input list of histogram. Check the TList name in the input file.
	TList *list = (TList*) file->Get("out1");
	
	// retrieve the histograms in the list. Check the name of the histograms
	histNum[0] = (TH2F*)list->FindObject("fHistElectronRec_EtaE_EmcalMC");
	histNum[1] = (TH2F*)list->FindObject("fHistConvElectronRec_EtaE_EmcalMC");
	histNum[2] = (TH2F*)list->FindObject("fHistScatElectronRec_EtaE_EmcalMC");
	histNum[3] = (TH2F*)list->FindObject("fHistGammaRec_EtaE_EmcalMC");
	histNum[4] = (TH2F*)list->FindObject("fHistAnnihGammaRec_EtaE_EmcalMC");
	histNum[5] = (TH2F*)list->FindObject("fHistScatGammaRec_EtaE_EmcalMC");
	
	// retrieve the histograms in the list. Check the name of the histograms
	histDen[0] = (TH2F*)list->FindObject("fHistElectronAcc_EtaE_EmcalMC");
	histDen[1] = (TH2F*)list->FindObject("fHistConvElectronAcc_EtaE_EmcalMC");
	histDen[2] = (TH2F*)list->FindObject("fHistScatElectronAcc_EtaE_EmcalMC");
	histDen[3] = (TH2F*)list->FindObject("fHistGammaAcc_EtaE_EmcalMC");
	histDen[4] = (TH2F*)list->FindObject("fHistAnnihGammaAcc_EtaE_EmcalMC");
	histDen[5] = (TH2F*)list->FindObject("fHistScatGammaAcc_EtaE_EmcalMC");
	
	// ********************************************

	Float_t x[fgNumOfEBins]={0}, ex[fgNumOfEBins]={0};
	Float_t y[fgNumOfEBins]={0}, ey[fgNumOfEBins]={0};
	Float_t num[fgNumOfEBins]={0}, den[fgNumOfEBins]={0};
	
	// loop over different desired outputs
	for (int iOut=0; iOut<NOUTPUTS; iOut++)
	{
		sprintf(efficName,"effic_%d",iOut);
		effic[iOut] = new TEfficiency(efficName,efficName,fgNumOfEBins,fgEAxis);

		// loop over E bins
		for (int ix=0; ix<fgNumOfEBins; ix++)
		{
			//check minimum energy
			if (histNum[0]->GetXaxis()->GetBinLowEdge(ix+1) < E_MIN)
				continue;
			
			// initialize ET variables for a new particle species
			x[ix]=histNum[0]->GetXaxis()->GetBinCenter(ix+1);
			y[ix]=0;
			ex[ix]=0;
			ey[ix]=0;
			num[ix] = 0;
			den[ix] = 0;
			
			// loop over eta bins
			for (int iy=0; iy<fgNumOfEtaBins; iy++)
			{
				for (int iHist=0; iHist<NHISTOUT[iOut]; iHist++)
				{
					num[ix] += histNum[IHISTOUT[iOut][iHist]]->GetBinContent(ix+1,iy+1); // sum over all E bins in order to get total ET
					den[ix] += histDen[IHISTOUT[iOut][iHist]]->GetBinContent(ix+1,iy+1); // sum over all E bins in order to get total ET
				}
			}
			
			if ((num[ix]>0) && (den[ix]>0))
			{
				effic[iOut]->SetTotalEvents(ix,den[ix]);
				effic[iOut]->SetPassedEvents(ix,num[ix]);

				y[ix] = num[ix]/den[ix];
				ey[ix] = y[ix]*sqrt(1/num[ix]+1/den[ix]);				
			}
			else
			{
				y[ix] = 0;	
				ey[ix] = 0;
			}
			
		} // end of loop over E bins

		graph[iOut] = new TGraphErrors(fgNumOfEBins,x,y,ex,ey); // graphic of ET(>E_cut)/ET(total) for a given particle species and E cut
		graphNum[iOut] = new TGraphErrors(fgNumOfEBins,x,num,ex,ey); // graphic of ET(>E_cut)/ET(total) for a given particle species and E cut	
		graphDen[iOut] = new TGraphErrors(fgNumOfEBins,x,den,ex,ey); // graphic of ET(>E_cut)/ET(total) for a given particle species and E cut

	} // end of loop over different outputs

	
	// Draw the plot
	
	if (makeDraw)
	{
		for (int i=0; i<NOUTPUTS; i++)
		{		
			c2->cd(i);
			
			graphDen[i]->SetMarkerStyle(style[i]);
			graphDen[i]->SetMarkerColor(color[i]);
			graphDen[i]->SetLineColor(color[i]);
			graphDen[i]->SetFillColor(0);
			if (i == 0) 
			{
				graphDen[i]->GetXaxis()->SetTitle("E (GeV)");
				graphDen[i]->GetYaxis()->SetTitle("effic");
				//graphDen[i]->SetMaximum(1.0);
				graphDen[i]->SetMinimum(0.0);
				graphDen[i]->Draw("AP");
			}
			else 
				graphDen[i]->Draw("P");
			
			graphNum[i]->SetMarkerStyle(style[i]+4);
			graphNum[i]->SetMarkerColor(color[i]);
			graphNum[i]->SetLineColor(color[i]);
			graphNum[i]->SetFillColor(0);
			graphNum[i]->Draw("P");		
		}
		
		c1->cd();
		/*
		for (int i=0; i<NOUTPUTS; i++)
		{
			effic[i]->SetMarkerStyle(style[i]);
			effic[i]->SetMarkerColor(color[i]);
			effic[i]->SetLineColor(color[i]);
			effic[i]->SetFillColor(0);
			effic[i]->SetTitle("efficiency;E (GeV); #epsilon");
			if (i == 0) 
			{
				effic[i]->Draw("AP");
			}
			else 
				effic[i]->Draw("Psame");
		}
		*/
		for (int i=0; i<NOUTPUTS; i++)
		{
			graph[i]->SetMarkerStyle(style[i]);
			graph[i]->SetMarkerColor(color[i]);
			graph[i]->SetLineColor(color[i]);
			graph[i]->SetFillColor(0);
			if (i == 0) 
			{
				graph[i]->Draw("AP");
			}
			else 
				graph[i]->Draw("P");
		}
		
		TLegend *leg = new TLegend(0.65,0.2,0.95,0.5);
		leg->AddEntry(effic[0],"electrons+gammas");
		leg->AddEntry(effic[1],"electrons");
		leg->AddEntry(effic[2],"gammas");
		leg->SetFillStyle(0);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->Draw();
	}
	
	return effic[0];
}
