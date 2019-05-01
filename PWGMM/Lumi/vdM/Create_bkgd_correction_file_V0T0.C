

#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Compute the background correction factor for each scan
//-------------------------------------------------------

void Create_bkgd_correction_file(Int_t scan_type, Int_t scan, Int_t rate_type)
// scan_type: 1 => x-scan; 2 => y-scan
// rate_type: 1 => V0; 2 => T0

// There are two steps"
// 1.- Find the start and end of each step in the scan
// 2.- loop over the range of each step and get the correction and its error
{

  // ** 1.- Find the start and end of each step in the scan ** //
  // -- find the number of separations (ie of steps in the scan)
  Int_t n_separations = FindNumberSeparations(scan_type, scan);
  // -- reserve memory to store the start and end of each step
  Int_t *idx_start = new Int_t[n_separations];
  Int_t *idx_end = new Int_t[n_separations];
  // -- find indices of steps
  FindStepStartAndEnd(scan_type, scan, n_separations, idx_start, idx_end);

  // ** 2.- loop over the range of each step and get the correction and its error ** //
   
  // set up branch addresses for incoming data
  Int_t aqflag;
  Double_t *acc_all = new Double_t[3564];
  Double_t *tot_all = new Double_t[3564];
  g_vdm_Tree->ResetBranchAddresses();
  g_vdm_Tree->SetBranchAddress("aqflag",&aqflag);
  if (rate_type == 1) {
    g_vdm_Tree->SetBranchAddress("bv0acc_setv0",acc_all);      
    g_vdm_Tree->SetBranchAddress("bv0tot_setv0",tot_all);
  } else if  (rate_type == 2) {
    g_vdm_Tree->SetBranchAddress("bt0acc_sett0",acc_all);      
    g_vdm_Tree->SetBranchAddress("bt0tot_sett0",tot_all);
  } 
  // set up output tree for rates
  char *file_name = new char[kg_string_size];
  if (scan_type == 1) {
    if (rate_type == 1) {
      sprintf(file_name,"../Fill-%d/BkgdCorr_V0_x_Scan_%d.root",
	      g_vdm_Fill,scan);
    } else if  (rate_type == 2) {
      sprintf(file_name,"../Fill-%d/BkgdCorr_T0_x_Scan_%d.root",
	      g_vdm_Fill,scan);
    }
  } else if (scan_type == 2) {
    if (rate_type == 1) {
      sprintf(file_name,"../Fill-%d/BkgdCorr_V0_y_Scan_%d.root",
	      g_vdm_Fill,scan);
    } else if  (rate_type == 2) {
      sprintf(file_name,"../Fill-%d/BkgdCorr_T0_y_Scan_%d.root",
	      g_vdm_Fill,scan);
    }
  }

  TFile *CorrFile = new TFile(file_name,"recreate");
  Double_t *correction = new Double_t[n_separations];
  Double_t *correction_error = new Double_t[n_separations];  
  TTree *correction_tree = new TTree("BkgdCorr","BkgdCorr");
  char txt_tmp[kg_string_size];
  sprintf(txt_tmp,"correction[%d]/D",n_separations);
  correction_tree->Branch("correction",correction,txt_tmp);
  sprintf(txt_tmp,"correction_error[%d]/D",n_separations);
  correction_tree->Branch("correction_error",correction_error,txt_tmp);

  // get bunch crossing information
  // -- number of bc
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  // -- bucket info
  Int_t *bunches = new Int_t [nIBC];
  GetBunchIndices(bunches);

  // loop over input data to fill the info
  for(Int_t k=0;k<nIBC;k++) { // loop over bunches
    if (k%10 == 0) cout << " Scan " << scan << "," << scan_type
		   << " working on bunch " << k << " of " <<  nIBC<< endl << flush;
    for(Int_t i=0;i<n_separations;i++) { // loop over steps
      // initalize
      correction_error[i]=correction[i]=0;
      Double_t acc = 0;
      Double_t tot = 0;
      // loop within a step
      for(Int_t j=idx_start[i];j<=idx_end[i];j++) {
	g_vdm_Tree->GetEntry(j);
	if (aqflag==0) continue;
	acc+=acc_all[bunches[k]];
	tot+=tot_all[bunches[k]];		
      }  // end loop within a step
      // compute correction and error
      correction[i] = BkgdCorrection(acc,tot);
      correction_error[i] = BkgdCorrectionError(acc,tot);
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
  delete [] acc_all;
  delete [] tot_all;  
  delete [] correction;
  delete [] correction_error;
  delete [] idx_start;
  delete [] idx_end;   
  delete [] file_name;
  delete [] bunches;  

}

void Create_bkgd_correction_file_V0T0(Int_t Fill, Int_t rate_type)
  // rate_type: 1 => V0; 2 => T0
{
  cout << " Starting, this will take some time " << endl;
  
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // create files for all scans
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    //x-scans
    Create_bkgd_correction_file(1, i, rate_type);
    //y-scans
    Create_bkgd_correction_file(2, i, rate_type);
  }
}
