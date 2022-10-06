#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the cross section
//-------------------------------------------------------

void QA_xs(Int_t Fill,
		const char *rate_name, const char *rate_type,
		const char *sep_type, const char *intensity_type,
		Int_t fit_type, Int_t scan
		)
{
	// initialize
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// get xs file and tree
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

	// get number of bunch crossings
	const Int_t nIBC = GetNumberInteractingBunchCrossings();
	Int_t Bunches[nIBC]; //kimc
	GetBunchIndices(Bunches);

	//Get bad bunches list, kimc
	SetBCBlacklists();//true);

	// reserve space
	Double_t *xs_all = new Double_t [nIBC];
	Double_t *xse_all = new Double_t [nIBC];
	TH1F* xs_H = new TH1F("xs_H", "", nIBC, -0.5, nIBC-0.5);
	TH1F* xs_Ho = (TH1F*)xs_H->Clone("xs_Hodd");
	TH1F* xs_He = (TH1F*)xs_H->Clone("xs_Heven");
	float yMin = 100;
	float yMax = 10;

	for (int i=0; i<nIBC; i++)
	{
		xs_tree->GetEntry(i);
		
		if (xs < 0) continue;
		if (OnBCBlacklist(Fill, Bunches[i])) //Rule out bad bunches, kimc
		{
			cout <<Form("Bad bunch detected, rule it out: %i (index %i)\n", Bunches[i], i);
			continue;
		}

		xs_H->SetBinContent(i+1, xs);
		xs_H->SetBinError(i+1, xs_error);
		if (Bunches[i]%2 != 0)
		{
			xs_Ho->SetBinContent(i+1, xs);
			xs_Ho->SetBinError(i+1, xs_error);
		}
		else
		{
			xs_He->SetBinContent(i+1, xs);
			xs_He->SetBinError(i+1, xs_error);
		}

		if (xs > yMax) yMax = xs * 1.05;
		if (xs < yMin) yMin = xs * 0.95;
	}

	// plot histo
	const char* TYPE = Form("Fill%i_%s_%s_%s_%s_fit%i_scan%i",
			Fill, sep_type, intensity_type, rate_name, rate_type, fit_type, scan); //kimc

	if      (!strcmp(rate_name, "VBAandVBC")) { yMin = 54; yMax = 58; }
	else if (!strcmp(rate_name, "TVX"))       { yMin = 26; yMax = 30; }

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetFitFormat("5.4f");
	TCanvas *bi_C = new TCanvas("xs_histo", "xs_histo", 1200, 800);

	xs_H->SetLineColor(0);
	xs_H->SetMarkerColor(0);
	xs_H->GetYaxis()->SetRangeUser(yMin, yMax);
	xs_H->SetTitle(Form("%s;BXing index;#sigma (mb)", TYPE));
	xs_H->DrawCopy("");
	xs_H->Fit("pol0");
	xs_Ho->SetLineColor(4);
	xs_Ho->SetMarkerColor(4);
	xs_Ho->SetMarkerStyle(20);
	xs_Ho->DrawCopy("pe same");
	xs_He->SetLineColor(210);
	xs_He->SetMarkerColor(210);
	xs_He->SetMarkerStyle(20);
	xs_He->DrawCopy("pe same");

	TLegend* L1 = new TLegend(0.15, 0.65, 0.35, 0.8);
	L1->AddEntry(xs_Ho, "Odd bunches", "lp");
	L1->AddEntry(xs_He, "Even bunches", "lp");
	L1->SetFillStyle(3001);
	L1->SetMargin(0.3);
	L1->Draw();
	bi_C->Print(Form("c3a_XSvsBC_%s.%s", TYPE, FFormat)); //kimc

	// clean
	delete [] xs_all;
	delete [] xse_all;  
}
