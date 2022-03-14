#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots for the
// beam center correction; that is, the ratio of the rate
// at mu to the rate at zero, where mu is the mean of the
// gaussian and the rate is evaluated using the model
//-------------------------------------------------------

void QA_bcc(Int_t Fill, const char *rate_name, const char *rate_type,
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
	Double_t par_x[npar];
	Double_t par_y[npar];  
	hx_tree->ResetBranchAddresses();
	hx_tree->SetBranchAddress("chi2_dof",&chi2_dof_x);  
	hx_tree->SetBranchAddress("par",par_x);
	hy_tree->ResetBranchAddresses();
	hy_tree->SetBranchAddress("chi2_dof",&chi2_dof_y);
	hy_tree->SetBranchAddress("par",par_y);

	// histograms
	TH1D *bcc_x_H = new TH1D("bcc_x_H","bcc horizontal scan",nIBC,-0.5,nIBC-0.5);
	TH1D *bcc_y_H = new TH1D("bcc_y_H","bcc horizontal scan",nIBC,-0.5,nIBC-0.5);

	// set up the model fit
	TF1 *fit_model_x = NULL;
	if      (fit_type == 0) { fit_model_x = new TF1("fit_model", fit_GP2, -1,1,Get_number_par(fit_type)); }
	else if (fit_type == 1) { fit_model_x = new TF1("fit_model", fit_GP6, -1,1,Get_number_par(fit_type)); }
	else if (fit_type == 2) { fit_model_x = new TF1("fit_model", fit_G,   -1,1,Get_number_par(fit_type)); }
	else if (fit_type == 4) { fit_model_x = new TF1("fit_model", fit_DG,  -1,1,Get_number_par(fit_type)); }
	else { 	cout << " Fit model " << fit_type << " not known " << endl; exit(-105); }
	TF1 *fit_model_y = NULL;
	if      (fit_type == 0) { fit_model_y = new TF1("fit_model", fit_GP2, -1,1,Get_number_par(fit_type)); }
	else if (fit_type == 1) { fit_model_y = new TF1("fit_model", fit_GP6, -1,1,Get_number_par(fit_type)); }
	else if (fit_type == 2) { fit_model_y = new TF1("fit_model", fit_G,   -1,1,Get_number_par(fit_type)); }
	else if (fit_type == 4) { fit_model_y = new TF1("fit_model", fit_DG,  -1,1,Get_number_par(fit_type)); }
	else { 	cout << " Fit model " << fit_type << " not known " << endl; exit(-105); }
	
	// fill info
	for(Int_t i=0;i<nIBC;i++) {
		hx_tree->GetEntry(i);
		hy_tree->GetEntry(i);
		fit_model_x->SetParameters(par_x);
		fit_model_y->SetParameters(par_y);
		if (chi2_dof_x>0) bcc_x_H->SetBinContent(i+1,(fit_model_x->Eval(par_x[1]))/(fit_model_x->Eval(0)));
		if (chi2_dof_y>0) bcc_y_H->SetBinContent(i+1,(fit_model_y->Eval(par_y[1]))/(fit_model_y->Eval(0)));
	}

	// plot histos
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(1);
	gStyle->SetOptFit(1);

	TCanvas *bcc_C = new TCanvas("bcc_C","bcc_C",800,800);
	bcc_C->Divide(1,2);
	bcc_C->cd(1);
	bcc_x_H->SetTitle("Horizontal scan;bunch crossing;f(#mu)/f(0)");
	bcc_x_H->SetMarkerStyle(20);
	bcc_x_H->Draw("p");
	bcc_x_H->Fit("pol0");
	bcc_C->cd(2);
	bcc_y_H->SetTitle("Vertical scan;bunch crossing;f(#mu)/f(0)");
	bcc_y_H->SetMarkerStyle(20);
	bcc_y_H->Draw("p");
	bcc_y_H->Fit("pol0");

}
