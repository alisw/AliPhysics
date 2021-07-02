#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"
#include "Plotting.h"


//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the cross section
//-------------------------------------------------------

TH1D* GetHistogram(Int_t Fill,  const char *rate_name, const char *rate_type,
                  const char *sep_type, const char *intensity_type, Int_t fit_type, Int_t scan )
{

    // initialize
    Set_input_file_names(Fill);
    Set_pointers_to_input_files_and_trees();
    
    // get hxhy file and tree
    char *hx_file_name = new char[kg_string_size];
    sprintf(hx_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_x_Scan_%d_Fit_%s.root",
            g_vdm_Fill,rate_type,rate_name,sep_type,scan,g_fit_model_name[fit_type]);
    char *hy_file_name = new char[kg_string_size];
    sprintf(hy_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_y_Scan_%d_Fit_%s.root",
            g_vdm_Fill,rate_type,rate_name,sep_type,scan,g_fit_model_name[fit_type]);
    
    // --> open files
    TFile *hx_file = new TFile(hx_file_name);
    TFile *hy_file = new TFile(hy_file_name);
    
    // --> get the trees
    TTree *hx_tree = (TTree *) hx_file->Get("AreaRate");
    TTree *hy_tree = (TTree *) hy_file->Get("AreaRate");
    
    // Next step: prepare variables to store the info for each bunch crossing
    // --> variables
    Double_t chi2_dof_x;
    Double_t chi2_dof_y;
    Double_t *rate_zero_x = new Double_t[2]; // rate at zero and its error
    Double_t *rate_zero_y = new Double_t[2]; // rate at zero and its error
    // --> set branches for hx, hy
    hx_tree->ResetBranchAddresses();
    hy_tree->ResetBranchAddresses();
    hx_tree->SetBranchAddress("chi2_dof",&chi2_dof_x);
    hy_tree->SetBranchAddress("chi2_dof",&chi2_dof_y);
    hx_tree->SetBranchAddress("rate_zero",rate_zero_x);
    hy_tree->SetBranchAddress("rate_zero",rate_zero_y);
    
    
    // get number of bunch crossings
    Int_t nIBC = GetNumberInteractingBunchCrossings();
    
    // reserve space
    TH1D *ratex_H = new TH1D("ratex_H","Rate zero x ",nIBC,-0.5,nIBC-0.5);
    TH1D *ratey_H = new TH1D("ratey_H","Rate zero y ",nIBC,-0.5,nIBC-0.5);
    
    // fill info
    for(Int_t i=0;i<nIBC;i++) {
        //Get info
        hx_tree->GetEntry(i);
        hy_tree->GetEntry(i);
        
        if (!(UseBunchCrossing( i))) continue;
         if (chi2_dof_x>0 && chi2_dof_y>0) { //check that fit is valid
             //fill histograms
             ratex_H->SetBinContent(i+1,rate_zero_x[0]);
             ratex_H->SetBinError(i+1,rate_zero_x[1]);
             ratey_H->SetBinContent(i+1,rate_zero_y[0]);
             ratey_H->SetBinError(i+1,rate_zero_y[1]);
         }
    }
    
    TH1D* ratio_H = (TH1D*)ratex_H->Clone("ratio");
    ratio_H->Divide(ratey_H);
    
    BeutifyTH1(ratio_H,Form("Scan %d;BC;Ratio: #frac{R^{x}(0,0)}{R^{y}(0,0)}",scan), 2, 1,20,1,1.2);
    
    return ratio_H;
    
}

void QA_headonrate(Int_t Fill, const char *rate_name,const char *rate_type, const char *sep_type,
              const char *intensity_type, Int_t fit_type,
              Bool_t save = kTRUE)
{
    
    TString name = "Ratio V0: #frac{R^{x}(0,0)}{R^{y}(0,0)} ";
    if (strncmp(rate_name,"TVX",3) == 0)
        name = "Ratio T0: #frac{R^{x}(0,0)}{R^{y}(0,0)} ";
    
    //Set histograms

    TH1D* ratio_0 = GetHistogram(Fill,rate_name,rate_type,sep_type,intensity_type,fit_type,0);
    TH1D* ratio_1 = GetHistogram(Fill,rate_name,rate_type,sep_type,intensity_type,fit_type,1);
  
    
    // plot histo
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetOptFit(1);
    gStyle->SetPadLeftMargin(0.15);
    TCanvas *bi_C = new TCanvas("ratio_histo","ratio_histo",1200,800);
    TPaveLabel* title = new TPaveLabel(0.06,0.92,0.9,0.99,name.Data());
    title->SetTextSize(0.5);
    title->Draw();
    TPad* graphPad = new TPad("Graphs","Graphs",0.01,0.05,0.95,0.90);
    graphPad->Draw();
    graphPad->cd();
    graphPad->Divide(2);
    

    graphPad->cd(1);
    
    ratio_0->SetMaximum(1.04); ratio_0->SetMinimum(0.96);
    ratio_0->Fit("pol0");
    ratio_0->Draw("p");
    
    graphPad->cd(2);
    ratio_1->SetMaximum(1.04); ratio_1->SetMinimum(0.96);
    ratio_1->Fit("pol0");
    ratio_1->Draw("p");
    
    // save plot
    if (save){
		const char* TYPE = Form("Fill%i_%s_%s_%s_%s_fit%i",
				Fill, sep_type, intensity_type, rate_name, rate_type, fit_type); //kimc
        TString plotName = Form("c3c_headonrate_%s.%s", TYPE, FFormat);
        //TString plotName = Form("../Fill-%d/Plots/%s_headonrate.pdf",Fill,rate_name);
        bi_C->SaveAs(plotName.Data());
    }
    
	return;
}





