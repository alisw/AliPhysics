#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Create the file with the normalisation for each  beam
//-------------------------------------------------------

void Create_beam_normalisation_tree(Int_t Fill, Int_t opt)
// opt = 0 => fBCT; opt = 1 => BPTX
{
	// get name of files and set pointers to trees
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// create tree for the output
	char *file_name = new char[kg_string_size];
	if      (opt == 0) sprintf(file_name, "../Fill-%d/FBCT_norm.root", g_vdm_Fill);
	else if (opt == 1) sprintf(file_name, "../Fill-%d/BPTX_norm.root", g_vdm_Fill);
	else { cout << " Option " << opt << " not recognised. Bye! " << endl; return; }

	TFile *NormFile = new TFile(file_name,"recreate");
	TTree *norm_tree = new TTree("Beam_Normalisation", "Beam_Normalisation");
	Double_t cf_dcct_1a = 0;
	Double_t cf_dcct_2a = 0;
	Double_t cf_dcct_1b = 0;
	Double_t cf_dcct_2b = 0;
	Double_t cf_dcct_1 = 0;
	Double_t cf_dcct_2 = 0;
	norm_tree->Branch("cf_dcct_1a", &cf_dcct_1a, "cf_dcct_1a/D");
	norm_tree->Branch("cf_dcct_1b", &cf_dcct_1b, "cf_dcct_1b/D");
	norm_tree->Branch("cf_dcct_1",  &cf_dcct_1,  "cf_dcct_1/D");
	norm_tree->Branch("cf_dcct_2a", &cf_dcct_2a, "cf_dcct_2a/D");
	norm_tree->Branch("cf_dcct_2b", &cf_dcct_2b, "cf_dcct_2b/D");
	norm_tree->Branch("cf_dcct_2",  &cf_dcct_2,  "cf_dcct_2/D");

	// get the branches with intensity information
	g_vdm_Tree->ResetBranchAddresses();
	Double_t timerel;
	Double_t *bunch1 = new Double_t[3564];  
	Double_t *bunch2 = new Double_t[3564];
	g_vdm_Tree->SetBranchAddress("timerel", &timerel);
	if (opt == 0)
	{
		g_vdm_Tree->SetBranchAddress("bunch1", bunch1); //kimc: this must be the FBCT intensity
		g_vdm_Tree->SetBranchAddress("bunch2", bunch2);
	}
	else if (opt == 1)
	{
		g_vdm_Tree->SetBranchAddress("bptx1", bunch1);
		g_vdm_Tree->SetBranchAddress("bptx2", bunch2);
	}

	// get the histograms with DCCT information, two devices: A and B
	TH1D* DCCT1AH = (TH1D*) g_vdm_File->Get("Beam1A");
	TH1D* DCCT2AH = (TH1D*) g_vdm_File->Get("Beam2A");
	TH1D* DCCT1BH = (TH1D*) g_vdm_File->Get("Beam1B");
	TH1D* DCCT2BH = (TH1D*) g_vdm_File->Get("Beam2B");  

	// loop over tree to normalise each timerel
	Int_t nEntries = g_vdm_Tree->GetEntries();
	for (Int_t i=0; i<nEntries; i++)
	{
		if (i%1000 == 0) cout << " Working on entry " << i << " of " << nEntries << endl;

		// get relative time
		g_vdm_Tree->GetEntry(i);

		// get the index of each histogram corresponding to the relative time
		Int_t idx_DCCT1AH = GetHistogramIndex(DCCT1AH, timerel);
		Int_t idx_DCCT1BH = GetHistogramIndex(DCCT1BH, timerel);  
		Int_t idx_DCCT2AH = GetHistogramIndex(DCCT2AH, timerel);
		Int_t idx_DCCT2BH = GetHistogramIndex(DCCT2BH, timerel);

		// get total intensity
		Double_t total_beam1 = 0;
		Double_t total_beam2 = 0;    
		for (Int_t j=0; j<3564; j++) // get total intensity - kimc: this must be related to a bunch crossing
		{
			total_beam1 += bunch1[j];
			total_beam2 += bunch2[j];      
		}

		// get correction factor 
		cf_dcct_1a = ((DCCT1AH->GetBinContent(idx_DCCT1AH))/total_beam1);
		cf_dcct_1b = ((DCCT1BH->GetBinContent(idx_DCCT1BH))/total_beam1);
		cf_dcct_2a = ((DCCT2AH->GetBinContent(idx_DCCT2AH))/total_beam2);
		cf_dcct_2b = ((DCCT2BH->GetBinContent(idx_DCCT2BH))/total_beam2);

		// correction factor is average over A and B factors
		cf_dcct_1 = 0.5*(cf_dcct_1a+cf_dcct_1b);
		cf_dcct_2 = 0.5*(cf_dcct_2a+cf_dcct_2b);

		// fill tree
		norm_tree->Fill();
	}

	// save  tree
	NormFile->cd();
	norm_tree->SetDirectory(NormFile);
	norm_tree->Write();
	NormFile->Close();

	// clean
	delete [] file_name;
	delete [] bunch1;
	delete [] bunch2;

	return;
}

