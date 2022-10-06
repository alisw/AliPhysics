#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the cross section
//-------------------------------------------------------

void QA_xs_N1N2(Int_t Fill,
		const char *rate_name, const char *rate_type,
		const char *sep_type, const char *intensity_type,
		Int_t fit_type, Int_t scan
		)
{
	// initialize
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// get number of bunch crossings
	const Int_t nIBC = GetNumberInteractingBunchCrossings();
	Int_t Bunches[nIBC]; //kimc
	GetBunchIndices(Bunches);

	//Get bad bunches list, kimc
	SetBCBlacklists();

	// --------------------
	// get intensity file and tree
	//---------------------
	// open the correct file
	char filename[120];
	sprintf(filename,"../Fill-%d/%s_Scan_%d.root",Fill,intensity_type,scan);
	TFile *IntensityFile = new TFile(filename);

	// get the info from the tree
	Double_t *bunch_intensity_1 = new Double_t[nIBC];
	Double_t *bunch_intensity_2 = new Double_t[nIBC];
	Double_t cf_dcct_1;
	Double_t cf_dcct_2;

	TTree *intensity_tree = (TTree *) IntensityFile->Get("Bunch_Intensity");
	intensity_tree->ResetBranchAddresses();
	intensity_tree->SetBranchAddress("cf_dcct_1",&cf_dcct_1);
	intensity_tree->SetBranchAddress("cf_dcct_2",&cf_dcct_2);    
	intensity_tree->SetBranchAddress("bunch_intensity_1",bunch_intensity_1);
	intensity_tree->SetBranchAddress("bunch_intensity_2",bunch_intensity_2);
	intensity_tree->GetEntry(0);

	// --------------------
	// get xs file and tree
	//---------------------
	char *xs_file_name = new char[kg_string_size];
	sprintf(xs_file_name,"../Fill-%d/xs_%sRate_%s_%sSep_%s_Scan_%d_Fit_%s.root",
			g_vdm_Fill,rate_type,rate_name,sep_type,intensity_type,scan,g_fit_model_name[fit_type]);
	TFile *xs_file = new TFile(xs_file_name);
	TTree *xs_tree = (TTree *) xs_file->Get("XS");

	// prepare tree branches
	Double_t xs=0;
	Double_t xs_error=0;  
	xs_tree->ResetBranchAddresses();
	xs_tree->SetBranchAddress("xs",&xs);
	xs_tree->SetBranchAddress("xs_error",&xs_error);    

	// reserve space
	Double_t *xs_all   = new Double_t [nIBC];
	Double_t *xse_all  = new Double_t [nIBC];
	Double_t *N1N2_all = new Double_t [nIBC];

	// define graph
	TGraphErrors* gr = new TGraphErrors();
	TGraphErrors* gr_o = new TGraphErrors();
	TGraphErrors* gr_e = new TGraphErrors(); 

	gr->SetMarkerColor(0);
	gr_o->SetLineColor(4);
	gr_o->SetMarkerColor(4);
	gr_o->SetMarkerStyle(20);
	gr_e->SetLineColor(210);
	gr_e->SetMarkerColor(210);
	gr_e->SetMarkerStyle(20);

	// fill info
	Int_t n = 0;
	float xMin = 0, xMax = 0;
	float yMin = 0, yMax = 0; //kimc
	for(Int_t i=0; i<nIBC; i++)
	{
		xs_tree->GetEntry(i);

		if (xs < 0) continue;
		if (OnBCBlacklist(Fill, Bunches[i])) //Rule out bad bunches, kimc
		{
			cout <<Form("Bad bunch detected, rule it out: %i (index %i)\n", Bunches[i], i);
			continue;
		}

		xs_all[i]   = xs;
		xse_all[i]  = xs_error;
		N1N2_all[i] = cf_dcct_1 * bunch_intensity_1[i] * cf_dcct_2 * bunch_intensity_2[i];

		//July 1, scale down x axis by 10^18 for PN not plot - kimc
		N1N2_all[i] = N1N2_all[i] / (1.E18);

		gr->SetPoint(gr->GetN(), N1N2_all[i], xs_all[i]);
		gr->SetPointError(gr->GetN()-1, 0, xse_all[i]);

		if (Bunches[i]%2 != 0)
		{
			gr_o->SetPoint(gr_o->GetN(), N1N2_all[i], xs_all[i]);
			gr_o->SetPointError(gr_o->GetN()-1, 0, xse_all[i]);
		}
		else
		{
			gr_e->SetPoint(gr_e->GetN(), N1N2_all[i], xs_all[i]);
			gr_e->SetPointError(gr_e->GetN()-1, 0, xse_all[i]);
		}

		if (yMin == 0.)
		{
			xMin = xMax = N1N2_all[i];
			yMin = yMax = xs_all[i];
		}
		else
		{
			if (N1N2_all[i] < xMin) xMin = N1N2_all[i]*0.97;
			if (N1N2_all[i] > xMax) xMax = N1N2_all[i]*1.03;
			if (xs_all[i] < yMin) yMin = xs_all[i]*0.95;
			if (xs_all[i] > yMax) yMax = xs_all[i]*1.05;
		}
	}

	#if 1
	//Plot for public notes (June 25, 2021)
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	float V0C = 0; 
	float T0C = 0;
	//if (Fill == 4937) { V0C = 58.5; T0C = 30.0; } //GlobalF
	if (Fill == 4937) { V0C = 57.0; T0C = 29.0; }
	if (Fill == 6012) { V0C = 56.5;	T0C = 27.7; }
	if (Fill == 6864) { V0C = 56.5;	T0C = 27.5; }
    if (!strcmp(rate_name, "VBAandVBC")) { yMin = V0C-1.0; yMax = V0C+1.5; }
    if (!strcmp(rate_name, "TVX"))       { yMin = T0C-1.0; yMax = T0C+1.5; }
	TCanvas *c1 = new TCanvas("c1","xs_histo", 1600, 800);
	c1->cd()->SetBottomMargin(0.125);

	TH1F* H1 = new TH1F("H1frame", ";N_{1} N_{2} #times 10^{ -18}; #sigma_{vis} [mb]", 100, xMin, xMax);
	H1->GetXaxis()->SetTitleOffset(1.3);
	H1->GetXaxis()->SetTitleSize(0.0375);
	H1->GetYaxis()->SetRangeUser(yMin, yMax);
	H1->GetYaxis()->SetTitleOffset(0.9);
	H1->GetYaxis()->SetTitleSize(0.05);
	H1->DrawCopy("");
	gPad->Update();

	gr->SetMarkerColor(1);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.5);
	gr->Draw("pe same");

	TF1* F1pol0 = new TF1("F1pol0", "pol0", xMin, xMax);
	F1pol0->SetParNames("Mean");
	F1pol0->SetLineColor(2);
	F1pol0->SetLineStyle(2);
	F1pol0->SetLineWidth(2);
	gr->Fit(F1pol0->GetName(), "ER", "", xMin, xMax);

	int YEAR = 0;
	if (Fill == 4937) YEAR = 2016;
	if (Fill == 6012) YEAR = 2017;
	if (Fill == 6864) YEAR = 2018;
	TLegend *L1 = new TLegend(0.145, 0.65, 0.35, 0.85);
	L1->SetMargin(0);
	L1->SetBorderSize(0);
	L1->SetFillStyle(3001);
	L1->AddEntry((TObject*)0, "ALICE", "");
	L1->AddEntry((TObject*)0, "pp #sqrt{s} = 13 TeV", "");
	L1->AddEntry((TObject*)0, Form("%i, %s, scan %i", YEAR, (!strcmp(rate_name,"VBAandVBC"))?"V0":"T0",scan), "");
	L1->Draw();

	TLegend *L2 = new TLegend(0.6, 0.735, 0.875, 0.85);
	L2->SetMargin(0.1);
	L2->SetNColumns(2);
	L2->AddEntry((TObject*)0, "#chi^{2}/NDF", "");
	L2->AddEntry((TObject*)0, Form("%5.2f/%i", F1pol0->GetChisquare(), F1pol0->GetNDF()), "");
	L2->AddEntry((TObject*)0, "Mean", "");
	L2->AddEntry((TObject*)0, Form("%5.3f #pm %4.3f", F1pol0->GetParameter(0), F1pol0->GetParError(0)), "");
	L2->Draw();

	c1->Print(Form("Fill%i_xsN1N2_%s_%i.eps", Fill, (!strcmp(rate_name,"VBAandVBC"))?"V0":"T0", scan));
	#endif

	#if 0
	// plot histo
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	//gStyle->SetOptTitle(0);

	const char* TYPE = Form("Fill%i_%s_%s_%s_%s_fit%i_scan%i",
			Fill, sep_type, intensity_type, rate_name, rate_type, fit_type, scan); //kimc

    if      (!strcmp(rate_name, "VBAandVBC")) { yMin = 54; yMax = 58; }
    else if (!strcmp(rate_name, "TVX"))       { yMin = 26; yMax = 30; }

	TCanvas *c1 = new TCanvas("xs_histo","xs_histo", 1200,800);
	c1->cd()->DrawFrame(xMin, yMin, xMax, yMax, Form("%s;N_{1}N_{2}; #sigma (mb)", TYPE));
	gr->Draw("pe same");
	gr->Fit("pol0");
	gr_o->Draw("pe same");
	gr_e->Draw("pe same");

	TLegend* L1 = new TLegend(0.15, 0.65, 0.35, 0.8);
    L1->AddEntry(gr_o, "Odd bunches", "lp");
    L1->AddEntry(gr_e, "Even bunches", "lp");
	L1->SetFillStyle(3001);
    L1->SetMargin(0.3);
    L1->Draw();	
	c1->Print(Form("c3b_XSvsN1N2_%s.%s", TYPE, FFormat)); //kimc
	#endif

	// clean
	delete [] xs_all;
	delete [] xse_all;  
	delete [] N1N2_all;
}
