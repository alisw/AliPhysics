
#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Check that the input trees have the same number of entries
//-------------------------------------------------------

void Check_size_of_input_trees(Int_t Fill)
{
  cout << " Comparing number of entries in the different trees " << endl;
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();
  // get number of entries
  Int_t n_vdm = (g_vdm_Tree ? g_vdm_Tree->GetEntries() : 0);
  Int_t n_vdm_DDL2 = (g_vdm_DDL2_Tree?g_vdm_DDL2_Tree->GetEntries() : 0);
  Int_t n_vdm_BPTX = (g_vdm_BPTX_Tree?g_vdm_BPTX_Tree->GetEntries() : 0);
  // compare number of entries
  if ((n_vdm != n_vdm_DDL2) || (n_vdm != n_vdm_BPTX) || (n_vdm_DDL2 != n_vdm_BPTX)) { 
    cout << " Entries in " << g_Input_vdm_File << ": " << n_vdm << endl;
    cout << " Entries in " << g_Input_vdm_DDL2_File << ": " << n_vdm_DDL2 << endl; 
    cout << " Entries in " << g_Input_vdm_BPTX_File << ": " << n_vdm_BPTX << endl;
  } else {
    cout << " All input trees have the same number of entries: " << n_vdm << endl;
  }
}

//-------------------------------------------------------
// Check that the input tree is ordered in time
//-------------------------------------------------------
void Check_time_order_of_input_tree(Int_t Fill)
{
  cout << " Checking time ordering in input tree " << endl;
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();
  // set time branches
  Double_t time;
  Double_t timerel;  
  g_vdm_Tree->SetBranchAddress("time",&time);  
  g_vdm_Tree->SetBranchAddress("timerel",&timerel);

  // loop over tree
  Double_t time_before = -1;
  Double_t timerel_before = -1;
  Bool_t not_ordered = kFALSE;
  Int_t n_vdm = g_vdm_Tree->GetEntries();
  for(Int_t i=0;i<n_vdm;i++) {
    g_vdm_Tree->GetEntry(i);
    if (time<time_before) {
      not_ordered = kTRUE;
      cout << " Entry " << i << " time_before " << time_before << "  time " << time << endl;
    }
    time_before = time;
    if (timerel<timerel_before) {
      cout << " Entry " << i << " timerel_before " << timerel_before << "  timerel " << timerel << endl;
      not_ordered = kTRUE;
    }
    timerel_before = timerel;
  }
  if (!not_ordered) cout << " Time order in tree corresponds to order of entries " << endl;
}

void Checks(Int_t Fill, Int_t opt=0)
{
  if (opt == 0) Check_size_of_input_trees(Fill);
  if (opt == 1) Check_time_order_of_input_tree(Fill);
}
