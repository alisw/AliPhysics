void norm(TH1F *h)
{
	Int_t Nbins = h->GetNbinsX();	
	for (Int_t i = 1 ; i < Nbins+1 ; i++ )
	{
		h->SetBinContent(i,h->GetBinContent(i) / ( h->GetBinWidth(i) * h->GetBinCenter(i) *1.6 * 2 * TMath::Pi()  ) );
		h->SetBinError(  i,h->GetBinError(i)   / ( h->GetBinWidth(i) * h->GetBinCenter(i) *1.6 * 2 * TMath::Pi()  ) );

	}
}

void setError(TH1F *h , TH1F *n , TH1F *d )
{
	Int_t Nbins = h->GetNbinsX();	
	for (Int_t i = 1 ; i < Nbins+1 ; i++ )
	{
		h->SetBinError(  i, ( n->GetBinError(i)/n->GetBinContent(i) + d->GetBinError(i)/d->GetBinContent(i) ) /2 );

	}

}


void format( TH1F *h, int color , int marker)
{

	h->SetMarkerColor(color);
	h->SetLineColor(color);
	h->SetMarkerStyle(marker);
	h->SetXTitle("#it{p}_{T} (GeV/#it{c}) ");
	h->SetYTitle("#frac{1}{N_{ev}2#pi#it{p}_{T}} #frac{dN^{2}}{d#it{p}_{T}d#eta}");
	h->GetXaxis()->SetRangeUser(0.15,10);

	h->GetYaxis()->SetLabelSize(.045);
	h->GetXaxis()->SetLabelSize(.045);


	h->GetYaxis()->SetTitleSize(.065);
	h->GetXaxis()->SetTitleSize(.075);

	h->GetYaxis()->SetTitleOffset(1);
	h->GetXaxis()->SetTitleOffset(1);
	h->SetMarkerSize(1);
	

}

void formatCorrelations(TH2F *h)
{
	h->GetZaxis()->SetRangeUser(0.0,0.0012);
	h->SetXTitle("tracklets");
	h->SetYTitle("SPD clusters");
	h->SetZTitle("#frac{Entries}{N_{ev}}");

	h->GetYaxis()->SetLabelSize(.045);
	h->GetXaxis()->SetLabelSize(.045);
	h->GetZaxis()->SetLabelSize(.045);

	h->GetYaxis()->SetTitleSize(.065);
	h->GetXaxis()->SetTitleSize(.075);
	h->GetZaxis()->SetTitleSize(.045);

	h->GetYaxis()->SetTitleOffset(1);
	h->GetXaxis()->SetTitleOffset(1);
	h->GetZaxis()->SetTitleOffset(1.5);
}

const float small = 0.00001;
const double margin = 0.195;
const double left = margin, right = margin;



void processLeadingPtQA(const char *filePath = "AnalysisResults.root",
				TString suffix="eps",
				const char *outfile="LeadingPt_output.root")
{
	const char *dirName = "outputLeadingPt";
	TString pwd(gSystem->pwd());
	gSystem->cd(pwd.Data()); 
	if(gSystem->cd(dirName)) {
		gSystem->cd(pwd.Data());
	} else {
		gSystem->mkdir(dirName, kTRUE); // kTRUE means recursive    }
	}

// obtencion de los histos

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TFile *f = TFile::Open(filePath,"read");
	TList* list = new TList();
	list= (TList*) f->Get("outputLeadingPt");

	TH1F *nchH[7] = {0},*nchL[7] = {0};
	TH1F *ptLH[7] = {0},*ptLL[7] = {0};
	TH1F *ptIH[7] = {0},*ptIL[7] = {0};

	TH1F *nchH_SPD[7] = {0},*nchL_SPD[7] = {0};
	TH1F *ptLH_SPD[7] = {0},*ptLL_SPD[7] = {0};
	TH1F *ptIH_SPD[7] = {0},*ptIL_SPD[7] = {0};

	TH1F *nchH_MV[7] = {0},*nchL_MV[7] = {0};
	TH1F *ptLH_MV[7] = {0},*ptLL_MV[7] = {0};
	TH1F *ptIH_MV[7] = {0},*ptIL_MV[7] = {0};

	TH1F *nchH_SPDMV[7] = {0},*nchL_SPDMV[7] = {0};
	TH1F *ptLH_SPDMV[7] = {0},*ptLL_SPDMV[7] = {0};
	TH1F *ptIH_SPDMV[7] = {0},*ptIL_SPDMV[7] = {0};

	for ( int i = 0 ; i < 7 ; i++ )
	{
		// no pileup rej
		nchH[i] = (TH1F*) list->FindObject( Form("fnchH%d",i) );
		nchH[i]->Sumw2();
		nchL[i] = (TH1F*) list->FindObject( Form("fnchL%d",i) );
		nchL[i]->Sumw2();
		ptLH[i] = (TH1F*) list->FindObject( Form("fptLH%d",i) );
		ptLH[i]->Sumw2();
		ptLL[i] = (TH1F*) list->FindObject( Form("fptLL%d",i) );
		ptLL[i]->Sumw2();
		ptIH[i] = (TH1F*) list->FindObject( Form("fptH%d",i) );
		ptIH[i]->Sumw2();
		ptIL[i] = (TH1F*) list->FindObject( Form("fptL%d",i) );
		ptIL[i]->Sumw2();
		// SPD pileup rej
		nchH_SPD[i] = (TH1F*) list->FindObject( Form("fnchH_SPD%d",i) );
		nchH_SPD[i]->Sumw2();
		nchL_SPD[i] = (TH1F*) list->FindObject( Form("fnchL_SPD%d",i) );
		nchL_SPD[i]->Sumw2();
		ptLH_SPD[i] = (TH1F*) list->FindObject( Form("fptLH_SPD%d",i) );
		ptLH_SPD[i]->Sumw2();
		ptLL_SPD[i] = (TH1F*) list->FindObject( Form("fptLL_SPD%d",i) );
		ptLL_SPD[i]->Sumw2();
		ptIH_SPD[i] = (TH1F*) list->FindObject( Form("fptH_SPD%d",i) );
		ptIH_SPD[i]->Sumw2();
		ptIL_SPD[i] = (TH1F*) list->FindObject( Form("fptL_SPD%d",i) );
		ptIL_SPD[i]->Sumw2();
		// MV pile up rej
		nchH_MV[i] = (TH1F*) list->FindObject( Form("fnchH_MV%d",i) );
		nchH_MV[i]->Sumw2();
		nchL_MV[i] = (TH1F*) list->FindObject( Form("fnchL_MV%d",i) );
		nchL_MV[i]->Sumw2();
		ptLH_MV[i] = (TH1F*) list->FindObject( Form("fptLH_MV%d",i) );
		ptLH_MV[i]->Sumw2();
		ptLL_MV[i] = (TH1F*) list->FindObject( Form("fptLL_MV%d",i) );
		ptLL_MV[i]->Sumw2();
		ptIH_MV[i] = (TH1F*) list->FindObject( Form("fptH_MV%d",i) );
		ptIH_MV[i]->Sumw2();
		ptIL_MV[i] = (TH1F*) list->FindObject( Form("fptL_MV%d",i) );
		ptIL_MV[i]->Sumw2();
		// MV + SPD pileup rej
		nchH_SPDMV[i] = (TH1F*) list->FindObject( Form("fnchH_SPDMV%d",i) );
		nchH_SPDMV[i]->Sumw2();
		nchL_SPDMV[i] = (TH1F*) list->FindObject( Form("fnchL_SPDMV%d",i) );
		nchL_SPDMV[i]->Sumw2();
		ptLH_SPDMV[i] = (TH1F*) list->FindObject( Form("fptLH_SPDMV%d",i) );
		ptLH_SPDMV[i]->Sumw2();
		ptLL_SPDMV[i] = (TH1F*) list->FindObject( Form("fptLL_SPDMV%d",i) );
		ptLL_SPDMV[i]->Sumw2();
		ptIH_SPDMV[i] = (TH1F*) list->FindObject( Form("fptH_SPDMV%d",i) );
		ptIH_SPDMV[i]->Sumw2();
		ptIL_SPDMV[i] = (TH1F*) list->FindObject( Form("fptL_SPDMV%d",i) );
		ptIL_SPDMV[i]->Sumw2();
	}

	TH1F *nchAll[7] = {0}      , *ptLAll[7] = {0}       , *ptIAll[7] = {0};
	TH1F *nchAll_SPD[7] = {0}  , *ptLAll_SPD[7] = {0}   , *ptIAll_SPD[7] = {0};
	TH1F *nchAll_MV[7] = {0}   , *ptLAll_MV[7] = {0}    , *ptIAll_MV[7] = {0};
	TH1F *nchAll_SPDMV[7] = {0}, *ptLAll_SPDMV[7] = {0} , *ptIAll_SPDMV[7] = {0};

	for ( int i = 0 ;  i < 7 ; i++ )
	{
		// No pileup rej
		nchAll[i] = (TH1F*) nchL[i]->Clone(Form("nchAll%d",i));
		nchAll[i]->Add(nchH[i]);

		ptLAll[i] = (TH1F*) ptLL[i]->Clone(Form("ptLAll%d",i));
		ptLAll[i]->Add(ptLH[i]);

		ptIAll[i] = (TH1F*) ptIL[i]->Clone(Form("ptIAll%d",i));
		ptIAll[i]->Add(ptIH[i]);
		//SPD pileup rej
		nchAll_SPD[i] = (TH1F*) nchL_SPD[i]->Clone(Form("nchAll_SPD%d",i));
		nchAll_SPD[i]->Add(nchH_SPD[i]);

		ptLAll_SPD[i] = (TH1F*) ptLL_SPD[i]->Clone(Form("ptLAll_SPD%d",i));
		ptLAll_SPD[i]->Add(ptLH_SPD[i]);

		ptIAll_SPD[i] = (TH1F*) ptIL_SPD[i]->Clone(Form("ptIAll_SPD%d",i));
		ptIAll_SPD[i]->Add(ptIH_SPD[i]);
		//MV pileup rej
		nchAll_MV[i] = (TH1F*) nchL_MV[i]->Clone(Form("nchAll_MV%d",i));
		nchAll_MV[i]->Add(nchH_MV[i]);

		ptLAll_MV[i] = (TH1F*) ptLL_MV[i]->Clone(Form("ptLAll_MV%d",i));
		ptLAll_MV[i]->Add(ptLH_MV[i]);

		ptIAll_MV[i] = (TH1F*) ptIL_MV[i]->Clone(Form("ptIAll_MV%d",i));
		ptIAll_MV[i]->Add(ptIH_MV[i]);
		//SPD + MV pileup rej
		nchAll_SPDMV[i] = (TH1F*) nchL_SPDMV[i]->Clone(Form("nchAll_SPDMV%d",i));
		nchAll_SPDMV[i]->Add(nchH_SPDMV[i]);

		ptLAll_SPDMV[i] = (TH1F*) ptLL_SPDMV[i]->Clone(Form("ptLAll_SPDMV%d",i));
		ptLAll_SPDMV[i]->Add(ptLH_SPDMV[i]);

		ptIAll_SPDMV[i] = (TH1F*) ptIL_SPDMV[i]->Clone(Form("ptIAll_SPDMV%d",i));
		ptIAll_SPDMV[i]->Add(ptIH_SPDMV[i]);
	}


	TH2F *trackVsClusters = 0,*trackVsClusters_SPD = 0,*trackVsClusters_MV = 0,*trackVsClusters_SPDMV = 0;
	
	trackVsClusters = (TH2F*) list->FindObject("ftrackVsClusters");
	trackVsClusters_SPD = (TH2F*) list->FindObject("ftrackVsClusters_SPD");
	trackVsClusters_MV = (TH2F*) list->FindObject("ftrackVsClusters_MV");
	trackVsClusters_SPDMV = (TH2F*) list->FindObject("ftrackVsClusters_SPDMV");

// normalizaciones

	for ( int i = 0 ; i < 7 ; i++ )
	{
		// No pileup rej
		norm(ptIAll[i]);
		ptIAll[i]->Scale(1.0/nchAll[i]->GetEntries());

		norm(ptIH[i]);
		ptIH[i]->Scale(1.0/nchH[i]->GetEntries());

		norm(ptLAll[i]);
		ptLAll[i]->Scale(1.0/ptLAll[i]->GetEntries());

		// No pileup rej
		norm(ptIAll_SPD[i]);
		ptIAll_SPD[i]->Scale(1.0/nchAll_SPD[i]->GetEntries());

		norm(ptIH_SPD[i]);
		ptIH_SPD[i]->Scale(1.0/nchH_SPD[i]->GetEntries());

		norm(ptLAll_SPD[i]);
		ptLAll_SPD[i]->Scale(1.0/ptLAll_SPD[i]->GetEntries());

		// No pileup rej
		norm(ptIAll_MV[i]);
		ptIAll_MV[i]->Scale(1.0/nchAll_MV[i]->GetEntries());

		norm(ptIH_MV[i]);
		ptIH_MV[i]->Scale(1.0/nchH_MV[i]->GetEntries());

		norm(ptLAll_MV[i]);
		ptLAll_MV[i]->Scale(1.0/ptLAll_MV[i]->GetEntries());

		// No pileup rej
		norm(ptIAll_SPDMV[i]);
		ptIAll_SPDMV[i]->Scale(1.0/nchAll_SPDMV[i]->GetEntries());

		norm(ptIH_SPDMV[i]);
		ptIH_SPDMV[i]->Scale(1.0/nchH_SPDMV[i]->GetEntries());

		norm(ptLAll_SPDMV[i]);
		ptLAll_SPDMV[i]->Scale(1.0/ptLAll_SPDMV[i]->GetEntries());
		
	}

	trackVsClusters->Scale(1./trackVsClusters->GetEntries());
	trackVsClusters_SPD->Scale(1./trackVsClusters_SPD->GetEntries());
	trackVsClusters_MV->Scale(1./trackVsClusters_MV->GetEntries());
	trackVsClusters_SPDMV->Scale(1./trackVsClusters_SPDMV->GetEntries());

// ratios

	TH1F *ratios[7] = {0};
	TH1F *ratios_SPD[7] = {0};
	TH1F *ratios_MV[7] = {0};
	TH1F *ratios_SPDMV[7] = {0};

	for ( int i = 0 ; i < 7 ; i++ )
	{	
		// No pileup rej
		ratios[i] = (TH1F*) ptIAll[i]->Clone(Form("ratios%d",i));
		ratios[i]->Divide(ptIH[i]);

		// No pileup rej
		ratios_SPD[i] = (TH1F*) ptIAll_SPD[i]->Clone(Form("ratios_SPD%d",i));
		ratios_SPD[i]->Divide(ptIH_SPD[i]);

		// No pileup rej
		ratios_MV[i] = (TH1F*) ptIAll_MV[i]->Clone(Form("ratios_MV%d",i));
		ratios_MV[i]->Divide(ptIH_MV[i]);

		// No pileup rej
		ratios_SPDMV[i] = (TH1F*) ptIAll_SPDMV[i]->Clone(Form("ratios_SPDMV%d",i));
		ratios_SPDMV[i]->Divide(ptIH_SPDMV[i]);
	}

//ratios leading

	TH1F *ratiosleading[4] = {0};		
	
	ratiosleading[0] = (TH1F*) ptLAll[6]->Clone("ratiosleading0");
	ratiosleading[1] = (TH1F*) ptLAll_SPD[6]->Clone("ratiosleading1");
	ratiosleading[2] = (TH1F*) ptLAll_MV[6]->Clone("ratiosleading2");
	ratiosleading[3] = (TH1F*) ptLAll_SPDMV[6]->Clone("ratiosleading3");

	for ( int i = 0 ; i < 4 ; i++ )
	{
		ratiosleading[i]->Divide(ptLAll[6]);
		format(ratiosleading[i],i+1,21);
	}

	// formatos

	format(ptLAll[6]      ,1,21);
	format(ptLAll_SPD[6]  ,2,21);
	format(ptLAll_MV[6]   ,3,21);
	format(ptLAll_SPDMV[6],4,21);


	setError(ratiosleading[0],ptLAll[6]      ,ptLAll[6]);
	setError(ratiosleading[1],ptLAll_SPD[6]  ,ptLAll[6]);
	setError(ratiosleading[2],ptLAll_MV[6]   ,ptLAll[6]);
	setError(ratiosleading[3],ptLAll_SPDMV[6],ptLAll[6]);

	formatCorrelations(trackVsClusters);
	formatCorrelations(trackVsClusters_SPD);
	formatCorrelations(trackVsClusters_MV);
	formatCorrelations(trackVsClusters_SPDMV);

// Leyendas

	TLegend *lleading = new TLegend(0.62,0.6,0.805,0.9);
	lleading->AddEntry(ptLAll[6],"w/o pileup rej (Ref)");
	lleading->AddEntry(ptLAll_SPD[6],"SPD pileup rej");
	lleading->AddEntry(ptLAll_MV[6],"MV pileup rej");
	lleading->AddEntry(ptLAll_SPDMV[6],"SPD+MV pileup rej");
	lleading->SetFillColor(0);


// graficas

	TCanvas *cspectra = new TCanvas("cspectra","pt spectra",1300,700);
	cspectra->Divide(1,2,small,small);
	cspectra->cd(1);
	cspectra->cd(1)->SetLogy();
	cspectra->cd(1)->SetGrid();
	//cspectra->cd(1)->SetLogx();
	gPad->SetLeftMargin(margin);
	gPad->SetRightMargin(margin);
	gPad->SetBottomMargin(small);

	ptLAll[6]->Draw("esame");
	ptLAll_SPD[6]->Draw("esame");
	ptLAll_MV[6]->Draw("esame");
	ptLAll_SPDMV[6]->Draw("esame");

	lleading->Draw();

	cspectra->cd(2);
	cspectra->cd(2)->SetGrid();
	//cspectra->cd(2)->SetLogx();
	gPad->SetTopMargin(small);
	gPad->SetLeftMargin(margin);
	gPad->SetRightMargin(margin);
	gPad->SetBottomMargin(margin);

	for ( int i = 0 ; i < 4 ; i++ )
	{
		ratiosleading[i]->Draw("same e2");
		ratiosleading[i]->SetYTitle("Ratios to Ref.    ");
		ratiosleading[i]->GetYaxis()->SetRangeUser(0.1,1.2);
		ratiosleading[i]->SetFillStyle(1);
	}

	cspectra->SaveAs(Form("%s/Pt_leading.%s",dirName,suffix.Data()));

	TFile *fout = TFile::Open(Form("%s/%s",dirName,outfile),"UPDATE");
	
	ptLAll[6]->Write();
	ptLAll_SPD[6]->Write();
	ptLAll_MV[6]->Write();
	ptLAll_SPDMV[6]->Write();

	for ( int i = 0 ; i < 4 ; i++ )
		ratiosleading[i]->Write();

	fout->Close();

	TCanvas *ccorrelations = new TCanvas("ccorrelations","tracklets vs SPD closters",1300,700);
	ccorrelations->Divide(2,2,small,small);

	ccorrelations->cd(1);
	gPad->SetTopMargin(margin);
	gPad->SetLeftMargin(margin);
	gPad->SetRightMargin(margin);
	gPad->SetBottomMargin(margin);	
	ccorrelations->cd(1)->SetGrid();
	trackVsClusters->Draw("colz");

	TPaveText *label1 = new TPaveText(80,600,175,695);
	label1->AddText("w/o pileup rejection");
	label1->SetFillColor(0);
	label1->Draw();

	ccorrelations->cd(2);
	gPad->SetTopMargin(margin);
	gPad->SetLeftMargin(margin);
	gPad->SetRightMargin(margin);
	gPad->SetBottomMargin(margin);
	ccorrelations->cd(2)->SetGrid();	
	trackVsClusters_SPD->Draw("colz");

	TPaveText *label2 = new TPaveText(80,600,175,695);
	label2->AddText("SPD pileup rejection");
	label2->SetFillColor(0);
	label2->Draw();

	ccorrelations->cd(3);
	gPad->SetTopMargin(margin);
	gPad->SetLeftMargin(margin);
	gPad->SetRightMargin(margin);
	gPad->SetBottomMargin(margin);	
	ccorrelations->cd(3)->SetGrid();
	trackVsClusters_MV->Draw("colz");

	TPaveText *label3 = new TPaveText(80,600,175,695);
	label3->AddText("MV pileup rejection");
	label3->SetFillColor(0);
	label3->Draw();

	ccorrelations->cd(4);
	gPad->SetTopMargin(margin);
	gPad->SetLeftMargin(margin);
	gPad->SetRightMargin(margin);
	gPad->SetBottomMargin(margin);	
	ccorrelations->cd(4)->SetGrid();
	trackVsClusters_SPDMV->Draw("colz");

	TPaveText *label4 = new TPaveText(80,600,175,695);
	label4->AddText("SPD+MV pileup rejection");
	label4->SetFillColor(0);
	label4->Draw();

	ccorrelations->SaveAs(Form("%s/SPDClustersVsTracklets.%s",dirName,suffix.Data()));

	TFile *fout = TFile::Open(Form("%s/%s",dirName,outfile),"UPDATE");
	
	ptLAll[6]->Write();
	ptLAll_SPD[6]->Write();
	ptLAll_MV[6]->Write();
	ptLAll_SPDMV[6]->Write();
	trackVsClusters->Write();
	trackVsClusters_SPD->Write();
	trackVsClusters_MV->Write();
	trackVsClusters_SPDMV->Write();
	
	

}
