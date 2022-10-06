#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Compute the intensity correction for each scan
//-------------------------------------------------------

void Compute_intensity_correction(Int_t Fill, Int_t scan_type, Int_t scan, Int_t intensity_type)
// scan_type: 1 => x-scan; 2 => y-scan
// intensity_type: 1 => BPTX; 2 => FBCT

// There are three steps
// 1.- Find the start and end of each step in the scan
// 2.- Upload the intensities to be used as normalization
// 3.- loop over the range of each step and get the correction and its error
{
	// ** 1.- Find the start and end of each step in the scan ** //
	// -- find the number of separations (ie of steps in the scan)
	Int_t n_separations = FindNumberSeparations(scan_type, scan);
	// -- reserve memory to store the start and end of each step
	Int_t *idx_start = new Int_t[n_separations];
	Int_t *idx_end   = new Int_t[n_separations];
	// -- find indices of steps
	FindStepStartAndEnd(scan_type, scan, n_separations, idx_start, idx_end);

	// ** 2.- Upload the intensities to be used as normalization ** //
	// (first get bunch crossing information)
	// -- number of bc
	const Int_t nIBC = GetNumberInteractingBunchCrossings();
	// -- index <-> bunch id, kimc
	// Int_t Bunches[nIBC];
	// GetBunchIndices(Bunches);
	// -- bucket info
	Int_t *BucketA = new Int_t [nIBC];
	Int_t *BucketC = new Int_t [nIBC];  
	GetBucketInfo(BucketA,BucketC);

	// --> create intensity file name
	char *intensity_file_name = new char[kg_string_size];
	if (intensity_type == 1) sprintf(intensity_file_name,"../Fill-%d/BPTX_Scan_%d.root",g_vdm_Fill,scan);
	if (intensity_type == 2) sprintf(intensity_file_name,"../Fill-%d/FBCT_Scan_%d.root",g_vdm_Fill,scan);
	// --> open file
	TFile *intensity_file = new TFile(intensity_file_name);
	// --> get intensity information
	TTree *intensity_tree = (TTree *) intensity_file->Get("Bunch_Intensity");
	Double_t *bunch_intensity_1 = new Double_t[nIBC];
	Double_t *bunch_intensity_2 = new Double_t[nIBC];  
	Double_t cf_dcct_1_ref = 0;
	Double_t cf_dcct_2_ref = 0;  
	intensity_tree->ResetBranchAddresses();
	intensity_tree->SetBranchAddress("bunch_intensity_1",bunch_intensity_1);
	intensity_tree->SetBranchAddress("bunch_intensity_2",bunch_intensity_2);
	intensity_tree->SetBranchAddress("cf_dcct_1",&cf_dcct_1_ref);
	intensity_tree->SetBranchAddress("cf_dcct_2",&cf_dcct_2_ref);
	intensity_tree->GetEntry(0);

	// --> set up tree with info on normalisation
	char *norm_file_name = new char[kg_string_size];
	if (intensity_type == 1) sprintf(norm_file_name,"../Fill-%d/BPTX_norm.root",g_vdm_Fill);
	if (intensity_type == 2) sprintf(norm_file_name,"../Fill-%d/FBCT_norm.root",g_vdm_Fill);
	// --> open file
	TFile *NormFile = new TFile(norm_file_name);      
	TTree *norm_tree = (TTree *) NormFile->Get("Beam_Normalisation");
	Double_t cf_dcct_1 = 0;
	Double_t cf_dcct_2 = 0;  
	norm_tree->ResetBranchAddresses();  
	norm_tree->SetBranchAddress("cf_dcct_1", &cf_dcct_1);
	norm_tree->SetBranchAddress("cf_dcct_2", &cf_dcct_2);

	// ** 3.- loop over the range of each step and get the correction and its error ** //
	// set up branch addresses for incoming data
	Int_t aqflag;
	Double_t *bunch1 = new Double_t[3564];  
	Double_t *bunch2 = new Double_t[3564];
	g_vdm_Tree->ResetBranchAddresses();
	g_vdm_Tree->SetBranchAddress("aqflag",&aqflag);
	if (intensity_type == 1)
	{
		g_vdm_Tree->SetBranchAddress("bptx1",bunch1);      
		g_vdm_Tree->SetBranchAddress("bptx2",bunch2);
	}
	else if (intensity_type == 2)
	{
		g_vdm_Tree->SetBranchAddress("bunch1",bunch1);      
		g_vdm_Tree->SetBranchAddress("bunch2",bunch2);
	}

	// set up output tree for rates
	char *file_name = new char[kg_string_size];
	if (scan_type == 1)
	{
		if      (intensity_type == 1) { sprintf(file_name,"../Fill-%d/BPTXCorr_x_Scan_%d.root", g_vdm_Fill,scan); }
		else if (intensity_type == 2) { sprintf(file_name,"../Fill-%d/FBCTCorr_x_Scan_%d.root",	g_vdm_Fill,scan); }
	}
	else if (scan_type == 2)
	{
		if      (intensity_type == 1) { sprintf(file_name,"../Fill-%d/BPTXCorr_y_Scan_%d.root", g_vdm_Fill,scan); }
		else if (intensity_type == 2) { sprintf(file_name,"../Fill-%d/FBCTCorr_y_Scan_%d.root",	g_vdm_Fill,scan); }
	}

	TFile *CorrFile = new TFile(file_name,"recreate");
	Double_t *b1_correction = new Double_t[n_separations];
	Double_t *b1_correction_error = new Double_t[n_separations];  
	Double_t *b2_correction = new Double_t[n_separations];
	Double_t *b2_correction_error = new Double_t[n_separations];  
	TTree *correction_tree = new TTree("IntensityCorrection","IntensityCorrection");
	char txt_tmp[kg_string_size];
	sprintf(txt_tmp,"b1_correction[%d]/D",n_separations);
	correction_tree->Branch("b1_correction",b1_correction,txt_tmp);
	sprintf(txt_tmp,"b1_correction_error[%d]/D",n_separations);
	correction_tree->Branch("b1_correction_error",b1_correction_error,txt_tmp);
	sprintf(txt_tmp,"b2_correction[%d]/D",n_separations);
	correction_tree->Branch("b2_correction",b2_correction,txt_tmp);
	sprintf(txt_tmp,"b2_correction_error[%d]/D",n_separations);
	correction_tree->Branch("b2_correction_error",b2_correction_error,txt_tmp);

	// loop over input data to fill the info
	for (Int_t k=0; k<nIBC; k++) // loop over bunches
	{
		if (k%10 == 0) cout <<Form("Scan %i, ScanType %i, Bunch %2i/%i\n", scan, scan_type, k, nIBC) <<flush;
		//if (OnBCBlacklist(Fill, Bunches[k])) continue;

		for (Int_t i=0; i<n_separations; i++) // loop over steps
		{
			// initalize to no correction and small error
			b1_correction[i]=b2_correction[i]=1.0;
			b1_correction_error[i]=b2_correction_error[i]=0.0001;

			// loop within a step
			Double_t bi1 = 0;
			Double_t bi2 = 0;
			Double_t bi1_2 = 0;
			Double_t bi2_2 = 0;
			Double_t ni = 0;

			for (Int_t j=idx_start[i]; j<=idx_end[i]; j++)
			{
				g_vdm_Tree->GetEntry(j);
				norm_tree->GetEntry(j);
				if (aqflag==0) continue;

				Int_t idx1 = ((Int_t)(BucketA[k]/10.0));
				Int_t idx2 = ((Int_t)(BucketC[k]/10.0)); 

				bi1 += (cf_dcct_1*bunch1[idx1]/gBeamA);
				bi2 += (cf_dcct_2*bunch2[idx2]/gBeamB);

				bi1_2 += ((cf_dcct_1*bunch1[idx1]/gBeamA)*(cf_dcct_1*bunch1[idx1]/gBeamA));
				bi2_2 += ((cf_dcct_2*bunch2[idx2]/gBeamB)*(cf_dcct_2*bunch2[idx2]/gBeamB));
				ni += 1.0;
			}//end loop within a step

			// avoid dividing by zero
			if (ni<3) continue;

			// compute correction and error
			b1_correction[i] = bi1/ni;
			b2_correction[i] = bi2/ni;
			b1_correction_error[i] = TMath::Sqrt((1./(ni-1.0))*(bi1_2-bi1*bi1/ni));
			b2_correction_error[i] = TMath::Sqrt((1./(ni-1.0))*(bi2_2-bi2*bi2/ni));

			b1_correction[i] = (b1_correction[i]/(cf_dcct_1_ref*bunch_intensity_1[k]));
			b1_correction_error[i] = (b1_correction_error[i]/(cf_dcct_1_ref*bunch_intensity_1[k]));

			b2_correction[i] = (b2_correction[i]/(cf_dcct_2_ref*bunch_intensity_2[k]));
			b2_correction_error[i] = (b2_correction_error[i]/(cf_dcct_2_ref*bunch_intensity_2[k]));

			// avoid extremely large correction factors
			if (b1_correction[i]>2.0 || b1_correction[i]<0.2)
			{
				b1_correction[i]=1.0; b1_correction_error[i]=0.0001;
				cout << " Extreme correction factor for b1 " << b1_correction[i]
					<< " separation " << i << " scan " << scan << endl;
			}
			if (b2_correction[i]>2.0 || b2_correction[i]<0.2)
			{
				b2_correction[i]=1.0; b2_correction_error[i]=0.0001;
				cout << " Extreme correction factor for b1 " << b1_correction[i]
					<< " separation " << i << " scan " << scan << endl;
			}
		} // end loop over steps

		// fill tree for each bunch
		correction_tree->Fill();
	} // end loop over bunches

	// save  tree
	CorrFile->cd();
	correction_tree->SetDirectory(CorrFile);
	correction_tree->Write();
	CorrFile->Close();

	// free memory
	delete [] idx_start;
	delete [] idx_end;   
	delete [] file_name;
	delete [] norm_file_name;  
	delete [] intensity_file_name;
	delete [] bunch_intensity_1; 
	delete [] bunch_intensity_2;
	delete [] bunch1;
	delete [] bunch2;
	delete [] b1_correction;
	delete [] b1_correction_error;
	delete [] b2_correction;
	delete [] b2_correction_error;
	delete [] BucketA;
	delete [] BucketC;

	return;
}

void Create_intensity_correction_file(Int_t Fill, Int_t intensity_type)
 // intensity_type: 1 => BPTX; 2 => FBCT
{  
	// get name of files and set pointers to trees
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// create files for all scans
	for (Int_t i=0; i<g_n_Scans_in_Fill; i++)
	{
		Compute_intensity_correction(Fill, 1, i, intensity_type); //x-scans
		Compute_intensity_correction(Fill, 2, i, intensity_type); //y-scans
	}
}
