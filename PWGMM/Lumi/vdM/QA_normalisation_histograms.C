#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Create the files with the information on the beam
//-------------------------------------------------------

void Create_single_normalisation_histogram(Int_t scan, Int_t opt)
// opt = 0 => fBCT; opt = 1 => BPTX
{
	// get the branches with intensity information (also, set the name of the output file)
	char *file_name = new char[kg_string_size];
	TFile *NormFile = NULL;
	g_vdm_Tree->ResetBranchAddresses();
	Double_t timerel;
	g_vdm_Tree->SetBranchAddress("timerel",&timerel);
	if (opt == 0) {	sprintf(file_name,"../Fill-%d/FBCT_norm.root",g_vdm_Fill); }
	if (opt == 1) { sprintf(file_name,"../Fill-%d/BPTX_norm.root",g_vdm_Fill); }

	// --> set up tree with info on normalisation
	NormFile = new TFile(file_name);      
	TTree *norm_tree = (TTree *) NormFile->Get("Beam_Normalisation");
	Double_t cf_dcct_1 = 0;
	Double_t cf_dcct_2 = 0;  
	norm_tree->ResetBranchAddresses();  
	norm_tree->SetBranchAddress("cf_dcct_1",&cf_dcct_1);
	norm_tree->SetBranchAddress("cf_dcct_2",&cf_dcct_2);

	// get indices for start and end of scan
	Int_t *indices = new Int_t [2];
	indices[0]=indices[1]=-1;
	FindIdicesOfScanStartEnd(scan, indices);

	// create histograms
	g_vdm_Tree->GetEntry(indices[0]);
	Double_t t0 = timerel;
	g_vdm_Tree->GetEntry(indices[1]);
	Double_t t1 = timerel;
	char *histo_name = new char[kg_string_size];
	sprintf(histo_name,"Norm1_%d_%d_%d",g_vdm_Fill,scan,opt);
	TH1D *h1 = new TH1D(histo_name,histo_name,((Int_t) (t1-t0)),t0,t1);
	sprintf(histo_name,"Norm2_%d_%d_%d",g_vdm_Fill,scan,opt);
	TH1D *h2 = new TH1D(histo_name,histo_name,((Int_t) (t1-t0)),t0,t1);

	// loop over tree to normalise each timerel
	Double_t min = 1e6;
	Double_t max = 1e-6;  
	for (Int_t i=indices[0]; i<indices[1]; i++)
	{
		// get relative time
		g_vdm_Tree->GetEntry(i);
		norm_tree->GetEntry(i);

		// fill histos
		Int_t ibin=h1->FindBin(timerel);
		h1->SetBinContent(ibin,cf_dcct_1);
		h2->SetBinContent(ibin,cf_dcct_2);

		if (min>cf_dcct_1) min = cf_dcct_1;
		if (min>cf_dcct_2) min = cf_dcct_2;    
		if (max<cf_dcct_1) max = cf_dcct_1;
		if (max<cf_dcct_2) max = cf_dcct_2;    

		/*
		// print info
		if ( fabs(cf_dcct_1 - 1.0) > 1.e-3 || fabs(cf_dcct_2 - 1.0) > 1.e-3 )
		{
			cout <<Form("timerel %4.0f (i = %i): cf_dcct_1 %5.4f cf_dcct_2 %5.4f\n",
					timerel,i,cf_dcct_1,cf_dcct_2);
		}
		*/
	}
	cout <<endl;

	// plot canvas
	sprintf(histo_name,"Fill%d_opt%d_scan%d", g_vdm_Fill, opt, scan);
	TCanvas *c = new TCanvas(histo_name,histo_name,900,600);
	c->Divide(1,2);
	c->cd(1);
	h1->SetTitle(";timerel;correction for beam 1");
	h1->SetMinimum(min*0.99);
	h1->SetMaximum(max*1.01);
	h1->Draw("p");
	c->cd(2);
	h2->SetTitle(";timerel;correction for beam 2");
	h2->SetMinimum(min*0.99);
	h2->SetMaximum(max*1.01);
	h2->Draw("p");
	c->Print(Form("c1b_sngNorm_%s.%s", c->GetName(), FFormat));

	// clean
	delete [] file_name;  
	delete [] histo_name;
	delete [] indices;
	return;
}

//-------------------------------------------------------
// Create normalisation histograms 
//-------------------------------------------------------

void QA_normalisation_histograms(Int_t Fill, Int_t opt)
// opt = 0 => fBCT; opt = 1 => BPTX
{
	// get name of files and set pointers to trees
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// create histograms for all scans
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	for(Int_t i=0;i<g_n_Scans_in_Fill;i++)
	{
		Create_single_normalisation_histogram(i,opt);
	}

	return;
}
