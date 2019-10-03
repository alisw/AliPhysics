#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"

//-------------------------------------------------------
// Compute hx or hy for a scan
//-------------------------------------------------------
void Compute_RateIntegral(Int_t scan_type, Int_t scan, const char *rate_name,
			  const char *rate_type, const char *sep_type, Int_t fit_type, Int_t bc)
// scan_type: 1 => x-scan; 2 => y-scan
//
// First get the right files and the trees inside them
// then loop over all bunch crossings and for each of them make a fit
{
  // first get the files and trees
  // --> create rate file name
  char *rate_file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(rate_file_name,"../Fill-%d/%sRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_type,rate_name,scan);
  if (scan_type == 2) sprintf(rate_file_name,"../Fill-%d/%sRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_type,rate_name,scan);

  // --> create separation file name
  char *sep_file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(sep_file_name,"../Fill-%d/%sSep_x_Scan_%d.root",g_vdm_Fill,sep_type,scan);
  if (scan_type == 2) sprintf(sep_file_name,"../Fill-%d/%sSep_y_Scan_%d.root",g_vdm_Fill,sep_type,scan);

  // --> open files
  TFile *rate_file = new TFile(rate_file_name);
  TFile *sep_file = new TFile(sep_file_name);
  // --> get the trees
  TTree *rate_tree = (TTree *) rate_file->Get("Rate");
  TTree *sep_tree = (TTree *) sep_file->Get("Separations");

   // Next step: prepare variables to be fitted
  // --> get number of separations
  Int_t n_sep = FindNumberSeparations(scan_type, scan);
  // --> reserve space for rates and separations
  Double_t *rate = new Double_t[n_sep];
  Double_t *rate_error = new Double_t[n_sep];  
  Double_t *sep = new Double_t[n_sep];
  // --> set branches
  rate_tree->ResetBranchAddresses();
  rate_tree->SetBranchAddress("rate",rate);
  rate_tree->SetBranchAddress("rate_error",rate_error);
  sep_tree->ResetBranchAddresses();
  sep_tree->SetBranchAddress("separation",sep);

  // Now reserve space for fit output
  Double_t *area = new Double_t[2]; // area and its error
  Double_t *rate_zero = new Double_t[2]; // rate at zero and its error
  Double_t chi2_dof = -1;
  Double_t *par = new Double_t[Get_number_par(fit_type)];
  Double_t *par_err = new Double_t[Get_number_par(fit_type)];

  // create file and tree for output
  // --> file
  char *h_file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(h_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_x_Scan_%d_Fit_%s.root",
			      g_vdm_Fill,rate_type,rate_name,sep_type,scan,g_fit_model_name[fit_type]);
  if (scan_type == 2) sprintf(h_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_y_Scan_%d_Fit_%s.root",
			      g_vdm_Fill,rate_type,rate_name,sep_type,scan,g_fit_model_name[fit_type]);  
  TFile *hFile = new TFile(h_file_name,"recreate");
  // --> tree
  TTree *h_tree = new TTree("AreaRate","AreaRate");
  char txt_tmp[kg_string_size];
  h_tree->Branch("chi2_dof",&chi2_dof,"chi2_dof/D");
  h_tree->Branch("area",area,"area[2]/D");
  h_tree->Branch("rate_zero",rate_zero,"rate_zero[2]/D");  
  sprintf(txt_tmp,"par[%d]/D",Get_number_par(fit_type));
  h_tree->Branch("par",par,txt_tmp);
  sprintf(txt_tmp,"par_err[%d]/D",Get_number_par(fit_type));
  h_tree->Branch("par_err",par_err,txt_tmp);

  // Next step: loop over the interacting bunch crossings (entries in the trees)
  // and make a fit to the separation vs rate curve
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  for(Int_t k=0;k<nIBC;k++) { // loop over bunches
    if (bc > -1 && bc != k) continue; // this is to study one specific bc ...
    // get info
    rate_tree->GetEntry(k);
    sep_tree->GetEntry(k);
    // do fit
    chi2_dof = Fit_rate_separation(n_sep,sep,rate,rate_error,fit_type,
				   area, rate_zero, par, par_err, scan, scan_type, bc);
    // save output
    h_tree->Fill();
  } // end loop over bunches

  // save
  hFile->cd();
  h_tree->SetDirectory(hFile);
  h_tree->Write();
  hFile->Close();

  // clean up
  delete [] h_file_name;    
  delete [] rate_file_name;
  delete [] sep_file_name;
  delete [] rate;
  delete [] rate_error;  
  delete [] sep;
  delete [] area;
  delete [] rate_zero;
  delete [] par;
  delete [] par_err;  
}


//-------------------------------------------------------
// Compute hx and hy for a fill
//-------------------------------------------------------

void Create_hxhy_file(Int_t Fill, const char *rate_name,
		      const char *rate_type, const char *sep_type, Int_t fit_type, Int_t bc = -1)
// if bc = -1, all bunches are analysed. if it is > -1, the corresponding bc is plotted
{
  // initialize
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // create files for all scans
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    //x-scans
    Compute_RateIntegral(1, i, rate_name, rate_type, sep_type, fit_type, bc);
    //y-scans
    Compute_RateIntegral(2, i, rate_name, rate_type, sep_type, fit_type, bc);
  }

}
