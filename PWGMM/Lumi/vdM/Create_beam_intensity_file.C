#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// get indices between x and y scans
//-------------------------------------------------------

void GetIndices(Int_t *idx)
{
	for (Int_t i=0; i<g_n_Scans_in_Fill; i++) idx[i] = FindIdxBetweenScans(i);
	return;
}

//-------------------------------------------------------
// Create the files with the information on the beam
//-------------------------------------------------------

void Create_single_beam_intensity_file(Int_t Fill, Int_t scan, Int_t opt)
// opt = 0 => fBCT; opt = 1 => BPTX
// first compute the normalization factors
// then compute the intensities at the reference point
{
	// get the branches with intensity and normalisation information
	// (also, set the name of the output file)
	TFile *NormFile = NULL;
	char *file_name = new char[kg_string_size];
	g_vdm_Tree->ResetBranchAddresses();
	Double_t *bunch1 = new Double_t[3564];  
	Double_t *bunch2 = new Double_t[3564];
	Double_t timerel;
	g_vdm_Tree->SetBranchAddress("timerel",&timerel);

	if (opt == 0)
	{
		g_vdm_Tree->SetBranchAddress("bunch1", bunch1);
		g_vdm_Tree->SetBranchAddress("bunch2", bunch2);
		sprintf(file_name, "../Fill-%d/FBCT_norm.root", g_vdm_Fill);

		NormFile = new TFile(file_name);
		sprintf(file_name, "../Fill-%d/FBCT_Scan_%d.root", g_vdm_Fill,scan);
	}
	else if (opt == 1)
	{
		g_vdm_Tree->SetBranchAddress("bptx1", bunch1);
		g_vdm_Tree->SetBranchAddress("bptx2", bunch2);
		sprintf(file_name, "../Fill-%d/BPTX_norm.root", g_vdm_Fill);

		NormFile = new TFile(file_name);
		sprintf(file_name, "../Fill-%d/BPTX_Scan_%d.root", g_vdm_Fill,scan);
	}
	else { cout <<"Unknown opt! Stop." <<endl; return; }

	// --> set up tree with info on normalisation
	TTree *norm_tree = (TTree*)NormFile->Get("Beam_Normalisation");
	Double_t cf_dcct_1 = 0;
	Double_t cf_dcct_2 = 0;  
	Double_t cf_dcct_1a = 0;
	Double_t cf_dcct_2a = 0;  
	Double_t cf_dcct_1b = 0;
	Double_t cf_dcct_2b = 0;
	norm_tree->ResetBranchAddresses();
	norm_tree->SetBranchAddress("cf_dcct_1", &cf_dcct_1);
	norm_tree->SetBranchAddress("cf_dcct_2", &cf_dcct_2);
	norm_tree->SetBranchAddress("cf_dcct_1a", &cf_dcct_1a);
	norm_tree->SetBranchAddress("cf_dcct_2a", &cf_dcct_2a);
	norm_tree->SetBranchAddress("cf_dcct_1b", &cf_dcct_1b);
	norm_tree->SetBranchAddress("cf_dcct_2b", &cf_dcct_2b);

	// get the index of the tree and get the entry
	Int_t idx_tree = FindIdxBetweenScans(scan);
	g_vdm_Tree->GetEntry(idx_tree);
	norm_tree->GetEntry(idx_tree);

	// print info
	cout << " For scan " << scan << " using idx_tree " << idx_tree << " with relative time " << timerel << endl;
	cout << " DCCT 1A = " << cf_dcct_1a << " DCCT 1B = " << cf_dcct_1b << endl;
	cout << " DCCT 2A = " << cf_dcct_2a << " DCCT 2B = " << cf_dcct_2b << endl;
	cout << " DCCT 1  = " << cf_dcct_1  << " DCCT 2  = " << cf_dcct_2  << endl;

	// get bunch crossing information
	// -- number of bc
	const Int_t nIBC = GetNumberInteractingBunchCrossings();
	// -- index <-> bunch id, kimc
	// Int_t Bunches[nIBC];
	// GetBunchIndices(Bunches);
	// -- bucket info
	Int_t *BucketA = new Int_t[nIBC];
	Int_t *BucketC = new Int_t[nIBC];  
	GetBucketInfo(BucketA, BucketC);

	// set up output tree for rates
	TFile *IntensityFile = new TFile(file_name,"recreate");
	TTree *intensity_tree = new TTree("Bunch_Intensity","Bunch_Intensity");
	Double_t *bunch_intensity_1 = new Double_t[nIBC];
	Double_t *bunch_intensity_2 = new Double_t[nIBC];
	for (Int_t i=0; i<nIBC; i++) bunch_intensity_1[i] = bunch_intensity_2[i] = 0.0;

	char *txt_tmp = new char[kg_string_size];
	sprintf(txt_tmp,"bunch_intensity_1[%d]/D", nIBC);
	intensity_tree->Branch("bunch_intensity_1", bunch_intensity_1,txt_tmp);
	sprintf(txt_tmp,"bunch_intensity_2[%d]/D" ,nIBC);
	intensity_tree->Branch("bunch_intensity_2", bunch_intensity_2,txt_tmp);
	intensity_tree->Branch("cf_dcct_1a",&cf_dcct_1a,"cf_dcct_1a/D");
	intensity_tree->Branch("cf_dcct_1b",&cf_dcct_1b,"cf_dcct_1b/D");
	intensity_tree->Branch("cf_dcct_1",&cf_dcct_1,"cf_dcct_1/D");      
	intensity_tree->Branch("cf_dcct_2a",&cf_dcct_2a,"cf_dcct_2a/D");
	intensity_tree->Branch("cf_dcct_2b",&cf_dcct_2b,"cf_dcct_2b/D");
	intensity_tree->Branch("cf_dcct_2",&cf_dcct_2,"cf_dcct_2/D");      

	// get the intensities
	g_vdm_Tree->GetEntry(idx_tree);
	for (Int_t j=0; j<nIBC; j++)
	{
		Int_t idx1 = ((Int_t)(BucketA[j]/10.0));
		Int_t idx2 = ((Int_t)(BucketC[j]/10.0));
		bunch_intensity_1[j] = (bunch1[idx1]/gBeamA);
		bunch_intensity_2[j] = (bunch2[idx2]/gBeamB);
	}
	intensity_tree->Fill();

	// save  tree
	IntensityFile->cd();
	intensity_tree->SetDirectory(IntensityFile);
	intensity_tree->Write();
	IntensityFile->Close();

	// free memory
	delete [] bunch1;
	delete [] bunch2;     
	delete [] bunch_intensity_1;
	delete [] bunch_intensity_2;     
	delete [] file_name;
	delete [] txt_tmp;  
	delete [] BucketA;
	delete [] BucketC;  

	return;
}

//-------------------------------------------------------
// Create the files with the information on the beam
// intensities for all bunch crossings 
//-------------------------------------------------------

void Create_beam_intensity_file(Int_t Fill, Int_t opt)
// opt = 0 => fBCT; opt = 1 => BPTX
{
	// get name of files and set pointers to trees
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// create files for all scans
	for (Int_t i=0; i<g_n_Scans_in_Fill; i++) Create_single_beam_intensity_file(Fill, i, opt);

	return;
}
