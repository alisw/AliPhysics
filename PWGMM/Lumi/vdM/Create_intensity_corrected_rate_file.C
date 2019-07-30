

#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Create a intensity corrected rate file for a given
// trigger, intensity correction and scan
//-------------------------------------------------------

void Create_one_intensity_corrected_rate_file(Int_t scan_type, Int_t scan, const char *rate_name, const char *intensity_corr_name)
// scan_type: 1 => x-scan; 2 => y-scan

{
  // create names for files
  char *file_name_rate = new char[kg_string_size];  
  char *file_name_corr = new char[kg_string_size];
  char *file_name_new_rate = new char[kg_string_size];  
  if (scan_type == 1) {
    sprintf(file_name_rate,"../Fill-%d/PileupCorrRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_name,scan);
    sprintf(file_name_corr,"../Fill-%d/%sCorr_x_Scan_%d.root",g_vdm_Fill,intensity_corr_name,scan);
    sprintf(file_name_new_rate,"../Fill-%d/IntensityCorr%sRate_%s_x_Scan_%d.root",g_vdm_Fill,intensity_corr_name,rate_name,scan);
  }
  if (scan_type == 2) {
    sprintf(file_name_rate,"../Fill-%d/PileupCorrRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_name,scan);
    sprintf(file_name_corr,"../Fill-%d/%sCorr_y_Scan_%d.root",g_vdm_Fill,intensity_corr_name,scan);
    sprintf(file_name_new_rate,"../Fill-%d/IntensityCorr%sRate_%s_y_Scan_%d.root",g_vdm_Fill,intensity_corr_name,rate_name,scan);
  }
  
  // open files and get trees
  TFile *RateFile = new TFile(file_name_rate);
  TFile *CorrFile = new TFile(file_name_corr);
  TTree *rate_tree = (TTree *) RateFile->Get("Rate");
  TTree *corr_tree = (TTree *) CorrFile->Get("IntensityCorrection");  

  // set branch adresses
  Int_t n_separations = FindNumberSeparations(scan_type, scan);
  Double_t *rate = new Double_t[n_separations];
  Double_t *rate_error = new Double_t[n_separations];  
  rate_tree->SetBranchAddress("rate",rate);
  rate_tree->SetBranchAddress("rate_error",rate_error);   
  Double_t *b1_correction = new Double_t[n_separations];
  Double_t *b1_correction_error = new Double_t[n_separations];  
  Double_t *b2_correction = new Double_t[n_separations];
  Double_t *b2_correction_error = new Double_t[n_separations];  
  corr_tree->SetBranchAddress("b1_correction",b1_correction);
  corr_tree->SetBranchAddress("b1_correction_error",b1_correction_error);   
  corr_tree->SetBranchAddress("b2_correction",b2_correction);
  corr_tree->SetBranchAddress("b2_correction_error",b2_correction_error);   

  // prepare new tree
  TFile *NewRateFile = new TFile(file_name_new_rate,"recreate");
  TTree *new_rate_tree = new TTree("Rate","Rate");
  Double_t *new_rate = new Double_t[n_separations];
  Double_t *new_rate_error = new Double_t[n_separations];
  char *txt_tmp = new char[kg_string_size];
  sprintf(txt_tmp,"rate[%d]/D",n_separations);
  new_rate_tree->Branch("rate",new_rate,txt_tmp);
  sprintf(txt_tmp,"rate_error[%d]/D",n_separations);
  new_rate_tree->Branch("rate_error",new_rate_error,txt_tmp);

  // loop over input data to compute new rate
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  for(Int_t k=0;k<nIBC;k++) { // loop over bunches
    // get entries
    rate_tree->GetEntry(k);
    corr_tree->GetEntry(k);
    // compute new rate
    for(Int_t i=0;i<n_separations;i++) { // loop over steps
      new_rate[i] = rate[i]/(b1_correction[i]*b2_correction[i]);
      new_rate_error[i] = rate_error[i]/(b1_correction[i]*b2_correction[i]);
    }
    // store the info
    new_rate_tree->Fill();
  }

  // save  tree
  NewRateFile->cd();
  new_rate_tree->SetDirectory(NewRateFile);
  new_rate_tree->Write();
  NewRateFile->Close();

  // free memory
  delete [] file_name_rate;
  delete [] file_name_corr;  
  delete [] file_name_new_rate;
  delete [] txt_tmp;  
  delete [] rate;
  delete [] rate_error;
  delete [] new_rate;
  delete [] new_rate_error;
  delete [] b1_correction;
  delete [] b1_correction_error;
  delete [] b2_correction;
  delete [] b2_correction_error;


}

//-------------------------------------------------------
// Create intensity corrected  rate files
//-------------------------------------------------------

void Create_intensity_corrected_rate_file(Int_t Fill, const char *rate_name, const char *intensity_corr_name)
{  
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // create files for all scans
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    //x-scans
    Create_one_intensity_corrected_rate_file(1, i, rate_name, intensity_corr_name);
    //y-scans
    Create_one_intensity_corrected_rate_file(2, i, rate_name, intensity_corr_name);    
  }
}
