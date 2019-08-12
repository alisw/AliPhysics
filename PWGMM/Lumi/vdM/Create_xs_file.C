
#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"

//-------------------------------------------------------
// Computes the cross section file for a scan
//-------------------------------------------------------
void Compute_xs(Int_t scan, const char *rate_name, const char *rate_type,
		const char *sep_type, const char *intensity_type, Int_t fit_type,
		Double_t lsc, Double_t ghost,
		Double_t satellite, Double_t non_fact)
{

  // first get the files and trees
  // --> create hx/hy file names
  char *hx_file_name = new char[kg_string_size];
  sprintf(hx_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_x_Scan_%d_Fit_%s.root",
	  g_vdm_Fill,rate_type,rate_name,sep_type,scan,g_fit_model_name[fit_type]);
  char *hy_file_name = new char[kg_string_size];
  sprintf(hy_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_y_Scan_%d_Fit_%s.root",
	  g_vdm_Fill,rate_type,rate_name,sep_type,scan,g_fit_model_name[fit_type]);
  // --> create intensity file name
  char *intensity_file_name = new char[kg_string_size];
  sprintf(intensity_file_name,"../Fill-%d/%s_Scan_%d.root",g_vdm_Fill,intensity_type,scan);

  // --> open files
  TFile *hx_file = new TFile(hx_file_name);
  TFile *hy_file = new TFile(hy_file_name);
  TFile *intensity_file = new TFile(intensity_file_name);
  // --> get the trees
  TTree *hx_tree = (TTree *) hx_file->Get("AreaRate");
  TTree *hy_tree = (TTree *) hy_file->Get("AreaRate");  
  TTree *intensity_tree = (TTree *) intensity_file->Get("Bunch_Intensity");

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

  // Next step: get intensity related information
  // --> get intensity information
  Int_t nIBC = GetNumberInteractingBunchCrossings();

  Double_t *bunch_intensity_1 = new Double_t[nIBC];
  Double_t *bunch_intensity_2 = new Double_t[nIBC];
  Double_t cf_dcct_1;
  Double_t cf_dcct_2;
  intensity_tree->ResetBranchAddresses();
  intensity_tree->SetBranchAddress("cf_dcct_1",&cf_dcct_1);
  intensity_tree->SetBranchAddress("cf_dcct_2",&cf_dcct_2);    
  intensity_tree->SetBranchAddress("bunch_intensity_1",bunch_intensity_1);
  intensity_tree->SetBranchAddress("bunch_intensity_2",bunch_intensity_2);
  intensity_tree->GetEntry(0);

  // Next step: compute global correction factor
  Double_t total_correction = lsc/(cf_dcct_1*cf_dcct_2*ghost*satellite*non_fact);
  
  // Next step: create file and tree for output
  // --> file
  char *xs_file_name = new char[kg_string_size];
  sprintf(xs_file_name,"../Fill-%d/xs_%sRate_%s_%sSep_%s_Scan_%d_Fit_%s.root",
	  g_vdm_Fill,rate_type,rate_name,sep_type,intensity_type,scan,g_fit_model_name[fit_type]);  
  TFile *xs_file = new TFile(xs_file_name,"recreate");
  // --> tree
  TTree *xs_tree = new TTree("XS","XS");
  xs_tree->Branch("chi2_dof_x",&chi2_dof_x,"chi2_dof_x/D");
  xs_tree->Branch("chi2_dof_y",&chi2_dof_y,"chi2_dof_y/D");
  Double_t xs=0;
  Double_t xs_error=0;  
  xs_tree->Branch("xs",&xs,"xs/D");
  xs_tree->Branch("xs_error",&xs_error,"xs_error/D");    
  xs_tree->Branch("lsc",&lsc,"lsc/D");
  xs_tree->Branch("ghost",&ghost,"ghost/D");  
  xs_tree->Branch("cf_dcct_1",&cf_dcct_1,"cf_dcct_1/D");
  xs_tree->Branch("cf_dcct_2",&cf_dcct_2,"cf_dcct_2/D");  
  xs_tree->Branch("satellite",&satellite,"satellite/D");
  xs_tree->Branch("non_fact",&non_fact,"non_fact/D");  

  // Next step: loop over the interacting bunch crossings (entries in the trees)
  // and compute the cross section
  for(Int_t k=0;k<nIBC;k++) { // loop over bunches
    xs = -1;
    xs_error = -1;
    // get info
    hx_tree->GetEntry(k);
    hy_tree->GetEntry(k);
    // check a valid fit
    if (chi2_dof_x>0 && chi2_dof_y>0) {
      // compute cross section
      xs = GetXS(area_x[0], area_y[0], rate_zero_x[0], rate_zero_y[0], bunch_intensity_1[k], bunch_intensity_2[k]);
      xs_error = GetXSerr(area_x[0], area_x[1],area_y[0], area_y[1],
			  rate_zero_x[0], rate_zero_x[1],rate_zero_y[0], rate_zero_y[1], 
			  bunch_intensity_1[k], bunch_intensity_2[k]);
      // correct cross section
      xs *= total_correction;
      xs_error *= total_correction;    
    //    cout << " scan " << scan << " bc " << k << " xs " << xs << " +/- " << xs_error << endl;
    }
    // save output
    xs_tree->Fill();
  } // end loop over bunches

  // save
  xs_file->cd();
  xs_tree->SetDirectory(xs_file);
  xs_tree->Write();
  xs_file->Close();

  // clean up
  delete [] xs_file_name;
  delete [] hx_file_name;
  delete [] hy_file_name;      
  delete [] intensity_file_name;
  delete [] area_x;
  delete [] rate_zero_x;
  delete [] area_y;
  delete [] rate_zero_y;
  delete [] bunch_intensity_1; 
  delete [] bunch_intensity_2; 

}


//-------------------------------------------------------
// Prepare the cross section file for a fill
//-------------------------------------------------------

void Create_xs_file(Int_t Fill, const char *rate_name, const char *rate_type,
		    const char *sep_type, const char *intensity_type, Int_t fit_type,
		    Double_t lsc, Double_t ghost,
		    Double_t satellite, Double_t non_fact)
{
  // initialize
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // create files for all scans
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    Compute_xs(i, rate_name, rate_type, sep_type,intensity_type, fit_type,
	       lsc, ghost, satellite, non_fact);
  }

}
