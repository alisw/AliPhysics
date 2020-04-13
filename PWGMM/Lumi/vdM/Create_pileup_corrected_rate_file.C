

#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Create a pileup corrected rate file for a given
// trigger,  and scan
//-------------------------------------------------------

void Create_one_pileup_corrected_rate_file(Int_t scan_type, Int_t scan, const char *rate_name, Double_t ratioA, Double_t ratioC)
// scan_type: 1 => x-scan; 2 => y-scan
// ratioA/ratioC -> parameters for the pileup correction

{
  // create names for files
  char *file_name_rate = new char[kg_string_size];  
  char *file_name_new_rate = new char[kg_string_size];  
  if (scan_type == 1) {
    sprintf(file_name_rate,"../Fill-%d/BkgdCorrRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_name,scan);
    sprintf(file_name_new_rate,"../Fill-%d/PileupCorrRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_name,scan);    
  }
  if (scan_type == 2) {
    sprintf(file_name_rate,"../Fill-%d/BkgdCorrRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_name,scan);
    sprintf(file_name_new_rate,"../Fill-%d/PileupCorrRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_name,scan);    
  }

  // open file and get tree
  TFile *RateFile = new TFile(file_name_rate);
  TTree *rate_tree = (TTree *) RateFile->Get("Rate");

  // set branch adresses
  Int_t n_separations = FindNumberSeparations(scan_type, scan);
  Double_t *rate = new Double_t[n_separations];
  Double_t *rate_error = new Double_t[n_separations];  
  rate_tree->SetBranchAddress("rate",rate);
  rate_tree->SetBranchAddress("rate_error",rate_error);   

  // prepare new tree
  TFile *NewRateFile = new TFile(file_name_new_rate,"recreate");
  TTree *ratio_tree = new TTree("Ratio","Ratio");
  ratio_tree->Branch("ratioA",&ratioA,"ratioA/D");
  ratio_tree->Branch("ratioC",&ratioC,"ratioC/D");
  ratio_tree->Fill(); // save parameters used in the correction
  TTree *new_rate_tree = new TTree("Rate","Rate");
  Double_t *new_rate = new Double_t[n_separations];
  Double_t *new_rate_error = new Double_t[n_separations];
  char *txt_tmp = new char[kg_string_size];
  sprintf(txt_tmp,"rate[%d]/D",n_separations);
  new_rate_tree->Branch("rate",new_rate,txt_tmp);
  sprintf(txt_tmp,"rate_error[%d]/D",n_separations);
  new_rate_tree->Branch("rate_error",new_rate_error,txt_tmp);

  // create a function to compute pileup
  // procedure and code from Martino
  TF1 * fPileUp = new TF1("fPileUp",GetPileUp,0.,4.,2);
  fPileUp->SetParameter(0,ratioA);
  fPileUp->SetParameter(1,ratioC);	


  // loop over input data to compute new rate
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  for(Int_t k=0;k<nIBC;k++) { // loop over bunches
    // get entries
    rate_tree->GetEntry(k);
    // compute new rate
    for(Int_t i=0;i<n_separations;i++) { // loop over steps
      new_rate[i] = RatePileUp(rate[i],fPileUp);
      new_rate_error[i] = RatePileUpErr(rate[i],rate_error[i],fPileUp);
    }
    // store the info
    new_rate_tree->Fill();
  }

  // save  tree
  NewRateFile->cd();
  new_rate_tree->SetDirectory(NewRateFile);
  new_rate_tree->Write();
  ratio_tree->SetDirectory(NewRateFile);
  ratio_tree->Write();
  NewRateFile->Close();

  // free memory
  delete [] file_name_rate;
  delete [] file_name_new_rate;
  delete [] txt_tmp;  
  delete [] rate;
  delete [] rate_error;
  delete [] new_rate;
  delete [] new_rate_error;

}

//-------------------------------------------------------
// Create pileup corrected  rate files
//-------------------------------------------------------

void Create_pileup_corrected_rate_file(Int_t Fill, const char *rate_name, Double_t ratioA, Double_t ratioC)
{
  cout << " This will take a while, be patient " << endl;
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // create files for all scans
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    //x-scans
    Create_one_pileup_corrected_rate_file(1, i, rate_name, ratioA, ratioC);
    //y-scans
    Create_one_pileup_corrected_rate_file(2, i, rate_name, ratioA, ratioC);    
  }
}
