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


TH1D* GetHistogram(const char *rate_name, const char *rate_type,
                  const char *sep_type, const char *intensity_type, Int_t fit_type, Int_t scan )
{
    
    // first get the files and trees
    // --> create hx/hy file names
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
    Double_t *area_x = new Double_t[2]; // area and its error
    Double_t *rate_zero_x = new Double_t[2]; // rate at zero and its error
    Double_t *area_y = new Double_t[2]; // area and its error
    Double_t *rate_zero_y = new Double_t[2]; // rate at zero and its error
    // --> set branches for hx, hy
    hx_tree->ResetBranchAddresses();
    hy_tree->ResetBranchAddresses();
    hx_tree->SetBranchAddress("chi2_dof",&chi2_dof_x);
    hy_tree->SetBranchAddress("chi2_dof",&chi2_dof_y);
    hx_tree->SetBranchAddress("area",area_x);
    hy_tree->SetBranchAddress("area",area_y);
    hx_tree->SetBranchAddress("rate_zero",rate_zero_x);
    hy_tree->SetBranchAddress("rate_zero",rate_zero_y);
    
    // get number of bunch crossings
    Int_t nIBC = GetNumberInteractingBunchCrossings();
    
    // reserve space
    TH1D *hxhy_H = new TH1D("hxhy_H","hxhy ",nIBC,-0.5,nIBC-0.5);
    
    // fill info
    for(Int_t i=0;i<nIBC;i++) {
        //Get info
        hx_tree->GetEntry(i);
        hy_tree->GetEntry(i);
        
        if (!(UseBunchCrossing( i))) continue;
        if (chi2_dof_x>0 && chi2_dof_y>0) {
            //fill histograms
            Double_t hxhy = GetHxHy(area_x[0],area_y[0],rate_zero_x[0],rate_zero_y[0]);
            Double_t hxhy_err = GetHxHyerr(area_x[0],area_x[1],area_y[0],area_y[1],rate_zero_x[0],rate_zero_x[1],rate_zero_y[0],rate_zero_y[1]);
            hxhy_H->SetBinContent(i+1,hxhy); cout <<"bin: "<< i<< " " << hxhy<< " +/-" <<  hxhy_err << endl;
            hxhy_H->SetBinError(i+1,hxhy_err);
        }
    }
    
    return hxhy_H;
    
}

void QA_hxhy_V0T0(Int_t Fill, const char *rate_type, const char *sep_type,
              const char *intensity_type, Int_t fit_type,
              Bool_t save = kTRUE)
{
    // initialize
    Set_input_file_names(Fill);
    Set_pointers_to_input_files_and_trees();
    TString name = "Ratio of h_{x}h_{y} #frac{T0}{V0}";
    
    //Set histograms
    TH1D* hxhy_V0_0 = GetHistogram("VBAandVBC",rate_type,sep_type,intensity_type,fit_type,0);
    TH1D* hxhy_V0_1 = GetHistogram("VBAandVBC",rate_type,sep_type,intensity_type,fit_type,1);
    TH1D* hxhy_T0_0 = GetHistogram("TVX",rate_type,sep_type,intensity_type,fit_type,0);
    TH1D* hxhy_T0_1 = GetHistogram("TVX",rate_type,sep_type,intensity_type,fit_type,1);

    TH1D* ratio[2];
	ratio[0] = (TH1D*)hxhy_T0_0->Clone("ratio_0");
    ratio[1] = (TH1D*)hxhy_T0_1->Clone("ratio_1");
    ratio[0]->Divide(hxhy_V0_0);
    ratio[1]->Divide(hxhy_V0_1);
    BeutifyTH1(ratio[0], Form(";Bunch pair ID; %s", name.Data()), 2, 1, 20, 1, 1.2);
    BeutifyTH1(ratio[1], Form(";Bunch pair ID; %s", name.Data()), 2, 1, 20, 1, 1.2);

#if 0    
    TH1D* ratio_0 = (TH1D*)hxhy_T0_0->Clone("ratio_0");
    TH1D* ratio_1 = (TH1D*)hxhy_T0_1->Clone("ratio_1");
    ratio_0->Divide(hxhy_V0_0);
    ratio_1->Divide(hxhy_V0_1);
    BeutifyTH1(ratio_0, Form("Scan 0;Bunch pair index; %s", name.Data()), 2, 1, 20, 1, 1.2);
    BeutifyTH1(ratio_1, Form("Scan 1;Bunch pair index; %s", name.Data()), 2, 1, 20, 1, 1.2);

    // plot histo
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetOptTitle(1);
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
    if (save)
	{
        const char* TYPE = Form("Fill%i_%s_%s_%s_fit%i",
                Fill, sep_type, intensity_type, rate_type, fit_type); //kimc
        TString plotName = Form("c3d_hxhyT0V0_%s.%s", TYPE, FFormat);

        //TString plotName = Form("../Fill-%d/hxhy_T0V0.pdf",Fill);
        bi_C->SaveAs(plotName.Data());
    }
#endif

	//Plot for public note, June 25
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	TCanvas* c1 = new TCanvas("c1_RatioT0V0", "", 1600, 1000); c1->Divide(1, 2);

	for (int a=0; a<2; a++)
	{
		c1->cd(a+1)->SetBottomMargin(0.125);
		ratio[a]->GetXaxis()->SetLabelSize(0.050);
		ratio[a]->GetXaxis()->SetTitleSize(0.050);
		ratio[a]->GetXaxis()->SetTitleOffset(1.2);
		ratio[a]->GetYaxis()->SetLabelSize(0.050);
		ratio[a]->GetYaxis()->SetTitleSize(0.055);
		ratio[a]->GetYaxis()->SetTitleOffset(0.8);
		ratio[a]->GetYaxis()->SetRangeUser(0.96, 1.04);
		ratio[a]->SetMarkerSize(1.4);
		ratio[a]->SetMarkerStyle(20);
		ratio[a]->DrawCopy("pe");

		const float x0 = ratio[a]->GetXaxis()->GetBinCenter(1);
		const float x1 = ratio[a]->GetXaxis()->GetBinCenter(ratio[a]->GetNbinsX());
		TF1* F1 = new TF1(Form("F1_ratio%i", a), "pol0", x0, x1);
		F1->SetParNames("Mean");
		F1->SetLineStyle(2);
		ratio[a]->Fit(F1->GetName(), "ER", "", x0, x1);

		TLegend *L1 = new TLegend(0.15, 0.7, 0.35, 0.825);
		L1->SetMargin(0);
		L1->SetBorderSize(0);
		L1->SetFillStyle(3001);
		L1->AddEntry((TObject*)0, "ALICE", "");
		L1->AddEntry((TObject*)0, "pp #sqrt{s} = 13 TeV", "");
		L1->Draw();

		TLegend *L2 = new TLegend(0.15, 0.2, 0.35, 0.275);
		L2->SetMargin(0);
		L2->SetBorderSize(0);
		L2->SetFillStyle(3001);
		L2->AddEntry((TObject*)0, Form("Scan %i", a), ""); 
		L2->Draw();

		TLegend *L3 = new TLegend(0.675, 0.72, 0.875, 0.85);
		L3->SetMargin(0.1);
		L3->SetNColumns(2);
		L3->AddEntry((TObject*)0, "#chi^{2}/NDF", "");
		L3->AddEntry((TObject*)0, Form("%5.2f/%i", F1->GetChisquare(), F1->GetNDF()), "");
		L3->AddEntry((TObject*)0, "Mean", "");
		L3->AddEntry((TObject*)0, Form("%5.3f #pm %4.3f", F1->GetParameter(0), F1->GetParError(0)), "");
		L3->Draw();
	}

	c1->Print(Form("Fill%i_RatioV0T0.eps", Fill));

	return;
}





