#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// accessing the rates

void GetRate(Double_t *rate, const char *rate_name, const char *rate_type, Int_t scan, Int_t scan_type, Int_t bc)
{
	char *file_name = new char[kg_string_size];
	if (scan_type == 1) sprintf(file_name,"../Fill-%d/%sRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_type,rate_name,scan);
	if (scan_type == 2) sprintf(file_name,"../Fill-%d/%sRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_type,rate_name,scan);
	TFile *rate_file = new TFile(file_name);
	TTree *rate_tree = (TTree *) rate_file->Get("Rate");
	rate_tree->ResetBranchAddresses();
	rate_tree->SetBranchAddress("rate",rate);
	rate_tree->GetEntry(bc);
	delete [] file_name;
}

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the corrections to the rate
// Note that it uses fBCT corrected rates.
//-------------------------------------------------------

void QA_corr_vs_sep(Int_t Fill, const char *rate_name, Int_t scan, Int_t scan_type, Int_t bc)
// scan_type: 1 => x-scan; 2 => y-scan
{
	// initialize
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// reserve space for rated
	Int_t n_sep = FindNumberSeparations(scan_type, scan);
	Double_t *raw = new Double_t[n_sep];
	Double_t *bkgd = new Double_t[n_sep];
	Double_t *pu = new Double_t[n_sep];
	Double_t *all = new Double_t[n_sep];      

	// get the rates
	GetRate(raw,rate_name,"Raw",scan,scan_type,bc);
	GetRate(bkgd,rate_name,"BkgdCorr",scan,scan_type,bc);
	GetRate(pu,rate_name,"PileupCorr",scan,scan_type,bc);
	GetRate(all,rate_name, Form("IntensityCorr%s", scan==0?"FBCT":"BPTX"), scan,scan_type,bc); //kimc
	//GetRate(all,rate_name,"IntensityCorrFBCT",scan,scan_type,bc);  

	// get the separations
	char *file_name = new char[kg_string_size];
	if (scan_type == 1) sprintf(file_name,"../Fill-%d/NomSep_x_Scan_%d.root",g_vdm_Fill,scan);
	if (scan_type == 2) sprintf(file_name,"../Fill-%d/NomSep_y_Scan_%d.root",g_vdm_Fill,scan);
	TFile *sep_file = new TFile(file_name);
	TTree *sep_tree = (TTree *) sep_file->Get("Separations");
	Double_t *sep = new Double_t[n_sep];
	sep_tree->ResetBranchAddresses();
	sep_tree->SetBranchAddress("separation",sep);
	sep_tree->GetEntry(bc);

	// print if needed
	// for(Int_t i=0;i<n_sep;i++) cout << i << " " << sep[i] << " " << raw[i] << " " << bkgd[i]<< " " << pu[i]<< " " << all [i] << endl;

	// make ratios
	Double_t *raw_bkgd = new Double_t[n_sep];
	Double_t *bkgd_pu = new Double_t[n_sep];
	Double_t *pu_all = new Double_t[n_sep];  
	for(Int_t i=0;i<n_sep;i++) {
		if (raw[i]>0 && bkgd[i]>0) raw_bkgd[i] = bkgd[i]/raw[i]; else raw_bkgd[i] = -1;
		if (bkgd[i]>0 && pu[i]>0) bkgd_pu[i] = pu[i]/bkgd[i]; else bkgd_pu[i] = -1;
		if (pu[i]>0 && all[i]>0) pu_all[i] = all[i]/pu[i]; else pu_all[i] = -1;    
	}

	// define the limits for the plot
	// --> separation
	Double_t sep_min = 0;
	Double_t sep_max = 0;
	for(Int_t i=0;i<n_sep;i++) {
		if(sep[i]<sep_min) sep_min=sep[i];
		if(sep[i]>sep_max) sep_max=sep[i];
	}
	sep_min*=1.2;
	sep_max*=1.2;  
	// --> ratio
	Double_t ratio_max = 0;
	Double_t ratio_min = 2;  
	for(Int_t i=0;i<n_sep;i++) {
		if(raw_bkgd[i] > -1 && raw_bkgd[i] > ratio_max) ratio_max=raw_bkgd[i];
		if(bkgd_pu[i] > -1 && bkgd_pu[i] > ratio_max) ratio_max=bkgd_pu[i];
		if(pu_all[i] > -1 && pu_all[i] > ratio_max) ratio_max=pu_all[i];    
		if(raw_bkgd[i] > -1 && raw_bkgd[i] < ratio_min) ratio_min=raw_bkgd[i];
		if(bkgd_pu[i] > -1 && bkgd_pu[i] < ratio_min) ratio_min=bkgd_pu[i];
		if(pu_all[i] > -1 && pu_all[i] < ratio_min) ratio_min=pu_all[i];    
	}
	ratio_max *= 1.5;
	ratio_min *= 0.8; 

	// make graphs
	TGraph *gr_raw_bkgd = new TGraph(n_sep,sep,raw_bkgd);
	gr_raw_bkgd->SetMarkerStyle(20);gr_raw_bkgd->SetMarkerColor(1);
	TGraph *gr_bkgd_pu = new TGraph(n_sep,sep,bkgd_pu);
	gr_bkgd_pu->SetMarkerStyle(24);gr_bkgd_pu->SetMarkerColor(2);  
	TGraph *gr_pu_all = new TGraph(n_sep,sep,pu_all);
	gr_pu_all->SetMarkerStyle(24);gr_pu_all->SetMarkerColor(4);

	/*
	// plot graphs
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TCanvas *corr_C = new TCanvas("corr_C","rate versus separation",600,400);
	TH1F* frame = gPad->DrawFrame(sep_min,ratio_min,sep_max,ratio_max);
	frame->SetTitle(";separation (mm); correction factor");
	gr_raw_bkgd->Draw("p,e1,same");
	gr_bkgd_pu->Draw("p,e1,same");
	gr_pu_all->Draw("p,e1,same");  
	TLegend *legend = new TLegend(0.3,0.7,0.7,0.9);
	legend->AddEntry(gr_raw_bkgd,"Raw/Bkgd","p");
	legend->AddEntry(gr_bkgd_pu,"Bkgd/Pileup","p");
	legend->AddEntry(gr_pu_all,"Pileup/Intensity","p");  
	legend->Draw();
	corr_C->Print(Form("c2b_QAcorrSep_Fill%i_%s_scanT%i_scan%i_bc%i.%s",
				Fill, rate_name, scan_type, scan, bc, FFormat));
				*/

    // Plot for public note, June 25
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1_corr", "", 800, 900); c1->cd();
    TH1F* H1 = new TH1F("H1frame", ";Separation [mm];Ratio", (sep_max - sep_min)*2, sep_min, sep_max);
	/*
	H1->SetTitle(Form("Fill %i, %s, scan %i, bunch pair index %i",
				Fill, !strcmp(rate_name, "TVX")?"T0":"V0", scan_type, bc));
				*/
    H1->GetXaxis()->SetRangeUser(-0.65, 0.65);
	if (Fill==4937) H1->GetYaxis()->SetRangeUser(0.7, 1.25);
	if (Fill==6012) H1->GetYaxis()->SetRangeUser(0.7, 1.30);
	if (Fill==6864) H1->GetYaxis()->SetRangeUser(0.7, 1.35);
    H1->DrawCopy();

    gr_raw_bkgd->SetMarkerSize(1.6);
    gr_raw_bkgd->SetMarkerStyle(28);
    gr_raw_bkgd->SetMarkerColor(1);
    gr_raw_bkgd->Draw("p,e1,same");
    gr_bkgd_pu->SetMarkerSize(1.6);
    gr_bkgd_pu->SetMarkerStyle(25);
    gr_bkgd_pu->SetMarkerColor(2);
    gr_bkgd_pu->Draw("p,e1,same");
    gr_pu_all->SetMarkerSize(1.6);
    gr_pu_all->SetMarkerStyle(24);
    gr_pu_all->SetMarkerColor(4);
    gr_pu_all->Draw("p,e1,same");

    TLegend *L1 = new TLegend(0.165, 0.75, 0.40, 0.85);
    L1->SetMargin(0);
    L1->SetBorderSize(0);
    L1->AddEntry((TObject*)0, "ALICE", "");
    L1->AddEntry((TObject*)0, "pp #sqrt{s} = 13 TeV", "");
    L1->Draw();

	TLegend* L2 = new TLegend(0.35, 0.2, 0.65, 0.4);
    L2->SetMargin(0.3);
    L2->SetBorderSize(0);
	L2->AddEntry(gr_raw_bkgd, "R_{BB} / R_{Raw}", "p");
	L2->AddEntry(gr_bkgd_pu, "R_{PU} / R_{BB}", "p");
	L2->AddEntry(gr_pu_all, "R_{DC} / R_{PU}", "p");
	L2->Draw();

	c1->Print(Form("Fill%i_Ratio_%s_%i.eps", Fill, !strcmp(rate_name, "TVX")?"T0":"V0", scan_type));

	// clean up
	delete [] raw;
	delete [] bkgd;
	delete [] pu;
	delete [] all;
	delete [] raw_bkgd;
	delete [] bkgd_pu;
	delete [] pu_all;  
}
