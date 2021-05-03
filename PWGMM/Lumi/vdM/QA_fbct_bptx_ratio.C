#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"
#include "Plotting.h"


TH1D* ratio_H[4] = {NULL};

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the cross section
//-------------------------------------------------------


void GetHistogram(Int_t Fill,Int_t scan)
{
	// initialize
    Set_input_file_names(Fill);
    Set_pointers_to_input_files_and_trees();

	// --> create intensity file name
    char *fbct_file_name = new char[kg_string_size];
    sprintf(fbct_file_name,"../Fill-%d/FBCT_Scan_%d.root",g_vdm_Fill,scan);
	TFile *fbct_file = new TFile(fbct_file_name);
	TTree *fbct_tree = (TTree *) fbct_file->Get("Bunch_Intensity");

	char *bptx_file_name = new char[kg_string_size];
	sprintf(bptx_file_name,"../Fill-%d/BPTX_Scan_%d.root",g_vdm_Fill,scan);
	TFile *bptx_file = new TFile(bptx_file_name);
	TTree *bptx_tree = (TTree *) bptx_file->Get("Bunch_Intensity");

	// get number of bunch crossings
	Int_t nIBC = GetNumberInteractingBunchCrossings();
	Int_t nIBCUsed = GetNumberOfUsedInteractingBunchCrossings();
	
  // --> get intensity information
	Double_t *fbct_bunch_intensity_1 = new Double_t[nIBC];
	Double_t *fbct_bunch_intensity_2 = new Double_t[nIBC];
	Double_t fbct_cf_dcct_1;
	Double_t fbct_cf_dcct_2;
	fbct_tree->ResetBranchAddresses();
    fbct_tree->SetBranchAddress("cf_dcct_1",&fbct_cf_dcct_1);
    fbct_tree->SetBranchAddress("cf_dcct_2",&fbct_cf_dcct_2);
	fbct_tree->SetBranchAddress("bunch_intensity_1",fbct_bunch_intensity_1);
	fbct_tree->SetBranchAddress("bunch_intensity_2",fbct_bunch_intensity_2);
	fbct_tree->GetEntry(0);

	Double_t *bptx_bunch_intensity_1 = new Double_t[nIBC];
	Double_t *bptx_bunch_intensity_2 = new Double_t[nIBC];
    Double_t bptx_cf_dcct_1;
    Double_t bptx_cf_dcct_2;
	bptx_tree->ResetBranchAddresses();
    bptx_tree->SetBranchAddress("cf_dcct_1",&bptx_cf_dcct_1);
    bptx_tree->SetBranchAddress("cf_dcct_2",&bptx_cf_dcct_2);
	bptx_tree->SetBranchAddress("bunch_intensity_1",bptx_bunch_intensity_1);
	bptx_tree->SetBranchAddress("bunch_intensity_2",bptx_bunch_intensity_2);
	bptx_tree->GetEntry(0);

	TH1D* fbct_beam1_H = new TH1D("fbct_beam1_H","fbct_beam1",nIBCUsed,-0.5,nIBCUsed-0.5);
	TH1D* fbct_beam2_H = new TH1D("fbct_beam2_H","fbct_beam2",nIBCUsed,-0.5,nIBCUsed-0.5);
	TH1D* bptx_beam1_H = new TH1D("bptx_beam1_H","fbct_beam1",nIBCUsed,-0.5,nIBCUsed-0.5);
	TH1D* bptx_beam2_H = new TH1D("bptx_beam2_H","bptx_beam2",nIBCUsed,-0.5,nIBCUsed-0.5);


	Int_t bin = 1;
	for (Int_t i=0;i<nIBC;i++) {
	  if (! UseBunchCrossing(i)) continue;
        fbct_beam1_H->SetBinContent(bin,fbct_cf_dcct_1*fbct_bunch_intensity_1[i]);
        fbct_beam2_H->SetBinContent(bin,fbct_cf_dcct_2*fbct_bunch_intensity_2[i]);
        bptx_beam1_H->SetBinContent(bin,bptx_cf_dcct_1*bptx_bunch_intensity_1[i]);
        bptx_beam2_H->SetBinContent(bin,bptx_cf_dcct_2*bptx_bunch_intensity_2[i]);
        bin++;
    }


    TH1D* ratio_beam1_H = (TH1D*)bptx_beam1_H->Clone("ratio_beam1_H");
	ratio_beam1_H->Divide(fbct_beam1_H);
	TH1D* ratio_beam2_H = (TH1D*)bptx_beam2_H->Clone("ratio_beam2_H");
	ratio_beam2_H->Divide(fbct_beam2_H);

    BeutifyTH1(ratio_beam1_H,Form("Scan %d Beam 1;BC;Ratio:BPTX/fBCT",scan), 2, 1,20,1,1.2);
	BeutifyTH1(ratio_beam2_H,Form("Scan %d Beam 2;BC;Ratio:BPTX/fBCT",scan), 2, 1,20,1,1.2);
    

	switch (scan) {
		case 0: ratio_H[0] = ratio_beam1_H;
                ratio_H[1] = ratio_beam2_H;
                break;
		case 1: ratio_H[2] = ratio_beam1_H;
                ratio_H[3] = ratio_beam2_H;
                break;
	}


}

void QA_fbct_bptx(Int_t Fill, Bool_t save = kFALSE)
{
	//Set histograms
	GetHistogram(Fill,0);
	GetHistogram(Fill,1);

    // plot histo
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
	gStyle->SetOptFit(1);
    TCanvas *bi_C = new TCanvas("ratio_histo","ratio_histo",1200,800);
	TString name = "Ratio: BPTX/FBCT";
	TPaveLabel* title = new TPaveLabel(0.06,0.96,0.9,0.99,name.Data());
	title->Draw();
	TPad* graphPad = new TPad("Graphs","Graphs",0.01,0.05,0.95,0.95);
	graphPad->Draw();
	graphPad->cd();
	graphPad->Divide(2,2);

	graphPad->cd(1);
	ratio_H[0]->SetMaximum(1.05); ratio_H[0]->SetMinimum(0.95);
	ratio_H[0]->Fit("pol0");
	ratio_H[0]->Draw("p");
	graphPad->cd(3);
	ratio_H[1]->SetMaximum(1.05); ratio_H[1]->SetMinimum(0.95);
	ratio_H[1]->Fit("pol0");
	ratio_H[1]->Draw("p");

	graphPad->cd(2);
	ratio_H[2]->SetMaximum(1.05); ratio_H[2]->SetMinimum(0.95);
	ratio_H[2]->Fit("pol0");
	ratio_H[2]->Draw("p");

	graphPad->cd(4);
	ratio_H[3]->SetMaximum(1.05); ratio_H[3]->SetMinimum(0.95);
	ratio_H[3]->Draw("p");
	ratio_H[3]->Fit("pol0");


  //chi2x_H->SetTitle("#chi^{2}/dof ;bunch crossing; #chi^{2}/dof (x)");

	// save plot
    if (save) {
        TString plotName = Form("../Fill-%d/Plots/FBCT_BPTX_ratio.pdf",Fill);
        bi_C->SaveAs(plotName.Data());
    }


}
