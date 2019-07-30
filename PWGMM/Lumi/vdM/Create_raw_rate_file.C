

#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Create the raw rate files for a given trigger and scan
//-------------------------------------------------------

void Create_DDL2_raw_rate_file(Int_t scan_type, Int_t scan, const char *rate_name)
// scan_type: 1 => x-scan; 2 => y-scan

// There are two steps"
// 1.- Find the start and end of each step in the scan
// 2.- loop over the range of each step and get the raw rate and its error
{
  // check if branch exists
  char *ddl2_counter = new char[kg_string_size];
  sprintf(ddl2_counter,"c2%s",rate_name);
  if (g_vdm_Tree->FindBranch(ddl2_counter) == NULL) {
    cout << " Requested branch ("<< ddl2_counter <<") does not exists in file " << g_Input_vdm_DDL2_File << endl;
    exit(-104);
  }

  // ** 1.- Find the start and end of each step in the scan ** //
  // -- find the number of separations (ie of steps in the scan)
  Int_t n_separations = FindNumberSeparations(scan_type, scan);
  // -- reserve memory to store the start and end of each step
  Int_t *idx_start = new Int_t[n_separations];
  Int_t *idx_end = new Int_t[n_separations];
  // -- find indices of steps
  FindStepStartAndEnd(scan_type, scan, n_separations, idx_start, idx_end);

  // ** 2.- loop over the range of each step and get the raw rate and its error ** //
   
  // set up branch addresses for incoming data
  Int_t orbit;
  Int_t aqflag;
  Double_t *raw_rate = new Double_t[3564];
  Double_t *counter = new Double_t[3564];
  Double_t *r2V0M = new Double_t[3564];//--jgc
  g_vdm_Tree->ResetBranchAddresses();
  g_vdm_Tree->SetBranchAddress("c2orbit",&orbit);
  g_vdm_Tree->SetBranchAddress("aqflag",&aqflag);
  g_vdm_Tree->SetBranchAddress("r2V0M",r2V0M);  //--jgc
  g_vdm_Tree->SetBranchAddress(ddl2_counter,counter);

  // set up output tree for rates
  char *file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/RawRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_name,scan);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/RawRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_name,scan);
  TFile *RateFile = new TFile(file_name,"recreate");
  Double_t *rate_orbits = new Double_t[n_separations];
  Double_t *rate = new Double_t[n_separations];
  Double_t *rate_error = new Double_t[n_separations];  
  Double_t *rate_counts = new Double_t[n_separations];
  TTree *rate_tree = new TTree("Rate","Rate");
  char *txt_tmp = new char[kg_string_size];
  sprintf(txt_tmp,"rate_orbits[%d]/D",n_separations);
  rate_tree->Branch("rate_orbits",rate_orbits,txt_tmp);
  sprintf(txt_tmp,"rate[%d]/D",n_separations);
  rate_tree->Branch("rate",rate,txt_tmp);
  sprintf(txt_tmp,"rate_error[%d]/D",n_separations);
  rate_tree->Branch("rate_error",rate_error,txt_tmp);
  sprintf(txt_tmp,"rate_counts[%d]/D",n_separations);
  rate_tree->Branch("rate_counts",rate_counts,txt_tmp);

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
      rate_orbits[i]=rate[i]=rate_error[i]=rate_counts[i]=0;
      // loop within a step
      for(Int_t j=idx_start[i];j<=idx_end[i];j++) {
	g_vdm_Tree->GetEntry(j);
	if (aqflag==0) continue;
	rate_orbits[i]+=orbit;
	rate_counts[i]+=counter[bunches[k]];
      }  // end loop within a step
      // compute rate and error
      rate[i]=RateRaw(rate_counts[i],rate_orbits[i]);
      rate_error[i]=RateRawErr(rate_counts[i],rate_orbits[i]);
    } // end loop over steps
    // fill tree for each bunch
    rate_tree->Fill();
  } // end loop over bunches

  // save  tree
  RateFile->cd();
  rate_tree->SetDirectory(RateFile);
  rate_tree->Write();
  RateFile->Close();

  // free memory
  delete [] rate_orbits;
  delete [] raw_rate;
  delete [] counter;
  delete [] rate;
  delete [] rate_error;
  delete [] rate_counts;
  delete [] idx_start;
  delete [] idx_end;   
  delete [] file_name;
  delete [] txt_tmp;  
  delete [] ddl2_counter;
  delete [] bunches;  

}

void Create_raw_rate_file(Int_t Fill, const char *rate_name)
{
  cout << " Starting, this will take some time " << endl;
  
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // create files for all scans
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    //x-scans
    Create_DDL2_raw_rate_file(1, i, rate_name);
    //y-scans
    Create_DDL2_raw_rate_file(2, i, rate_name);
  }
}
