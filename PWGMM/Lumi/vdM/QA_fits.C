#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the fits
//-------------------------------------------------------

void QA_fits(Int_t Fill, const char *rate_name, const char *rate_type,
	     const char *sep_type, Int_t fit_type,
	     Int_t scan)
{
	// initialize
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// get number of bunch crossings
	const Int_t nIBC = GetNumberInteractingBunchCrossings();

	// create hx/hy file names
	char *hx_file_name = new char[kg_string_size];
	sprintf(hx_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_x_Scan_%d_Fit_%s.root",
			g_vdm_Fill,rate_type,rate_name,sep_type,scan,g_fit_model_name[fit_type]);
	char *hy_file_name = new char[kg_string_size];
	sprintf(hy_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_y_Scan_%d_Fit_%s.root",
			g_vdm_Fill,rate_type,rate_name,sep_type,scan,g_fit_model_name[fit_type]);

	// open files
	TFile *hx_file = new TFile(hx_file_name);
	TFile *hy_file = new TFile(hy_file_name);

	// get the trees
	TTree *hx_tree = (TTree *) hx_file->Get("AreaRate");
	TTree *hy_tree = (TTree *) hy_file->Get("AreaRate");  

	// set branches
	const Int_t npar = Get_number_par(fit_type);
	Double_t chi2_dof_x; 
	Double_t chi2_dof_y; 
	Double_t Rmax_x[2]; 
	Double_t Rmax_y[2]; 
	Double_t par_x[npar];
	Double_t par_y[npar];  
	Double_t par_err_x[npar];
	Double_t par_err_y[npar];  
	hx_tree->ResetBranchAddresses();
	hx_tree->SetBranchAddress("chi2_dof",&chi2_dof_x);  
	hx_tree->SetBranchAddress("rate_zero",Rmax_x);   
	hx_tree->SetBranchAddress("par",par_x);
	hx_tree->SetBranchAddress("par_err",par_err_x);    
	hy_tree->ResetBranchAddresses();
	hy_tree->SetBranchAddress("chi2_dof",&chi2_dof_y);
	hy_tree->SetBranchAddress("rate_zero",Rmax_y); 
	hy_tree->SetBranchAddress("par",par_y);
	hy_tree->SetBranchAddress("par_err",par_err_y);

	// histograms
	TH1D *chi2_dof_x_H = new TH1D("chi2_dof_x_H","chi2 horizontal scan",nIBC,-0.5,nIBC-0.5);
	TH1D *chi2_dof_y_H = new TH1D("chi2_dof_y_H","chi2 vertical scan",nIBC,-0.5,nIBC-0.5);
	TH1D *par_x_H[npar];
	TH1D *par_y_H[npar];
	for (Int_t i=0;i<npar;i++) {
		par_x_H[i] = new TH1D(Form("par_x %d",i),Form("par_x %d",i),nIBC,-0.5,nIBC-0.5);
		par_y_H[i] = new TH1D(Form("par_y %d",i),Form("par_y %d",i),nIBC,-0.5,nIBC-0.5);    
	}
	// graphs
	TGraphErrors *gr = new TGraphErrors();

	// fill info
	for(Int_t i=0;i<nIBC;i++) {
		hx_tree->GetEntry(i);
		hy_tree->GetEntry(i);    
		chi2_dof_x_H->SetBinContent(i+1,chi2_dof_x); 
		chi2_dof_y_H->SetBinContent(i+1,chi2_dof_y);
		if (Rmax_y[0]>0 && Rmax_x[0]>0) {
			Int_t n = gr->GetN();
			gr->SetPoint(n,i,(Rmax_x[0]/Rmax_y[0]));
			gr->SetPointError(n,0,(Rmax_x[1]/Rmax_x[0]));
		}
		for(Int_t j=0;j<npar;j++) {
			par_x_H[j]->SetBinContent(i+1,par_x[j]);
			par_x_H[j]->SetBinError(i+1,par_err_x[j]);      
			par_y_H[j]->SetBinContent(i+1,par_y[j]);
			par_y_H[j]->SetBinError(i+1,par_err_y[j]);      
		}
	}

	// plot histos
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(1);
	gStyle->SetOptFit(1);

	TCanvas *chi2_C = new TCanvas("chi2_C","chi2_C",800,800);
	chi2_C->Divide(1,2);
	chi2_C->cd(1);
	chi2_dof_x_H->SetTitle("Horizontal scan;bunch crossing;#chi^{2}");
	chi2_dof_x_H->SetMarkerStyle(20);
	chi2_dof_x_H->Draw("p");
	chi2_C->cd(2);
	chi2_dof_y_H->SetTitle("Vertical scan;bunch crossing;#chi^{2}");
	chi2_dof_y_H->SetMarkerStyle(20);
	chi2_dof_y_H->Draw("p");

	for(Int_t j=0;j<npar;j++) {
		TCanvas *par_C = new TCanvas(Form("par_%d_C",j),Form("par_%d_C",j),800,800);
		par_C->Divide(1,2);
		par_C->cd(1);
		par_x_H[j]->SetTitle(Form("Vertical scan;bunch crossing;par[%d]",j));
		par_x_H[j]->SetMarkerStyle(20);
		par_x_H[j]->Draw("p,e");
		par_C->cd(2);
		par_y_H[j]->SetTitle(Form("Vertical scan;bunch crossing;par[%d]",j));
		par_y_H[j]->SetMarkerStyle(20);
		par_y_H[j]->Draw("p,e");
	}

	TCanvas *r_C = new TCanvas("r_C","r_C",800,400);
	gr->SetMarkerStyle(20);
	TH1* h = (TH1*) gr->GetHistogram();
	h->GetXaxis()->SetRangeUser(-0.5,nIBC-0.5);
	h->SetTitle(";bunch crossing;Rmax ratio ver/hor");
	gr->Draw("ap");
	gr->Fit("pol0");
}
