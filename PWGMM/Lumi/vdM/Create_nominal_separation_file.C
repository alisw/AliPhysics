
#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Find the nominal separations in the vdm tree
// within the given inidices 
//-------------------------------------------------------
void Find_separations(Int_t scan_type, Int_t scan_num, Int_t IdxStart, Int_t IdxEnd)
// scan_type: 1 => x-scan; 2 => y-scan
{

  // set  branches
  Double_t time, nsep;
  Int_t aqflag;
  g_vdm_Tree->ResetBranchAddresses();
  g_vdm_Tree->SetBranchAddress("aqflag",&aqflag);
  g_vdm_Tree->SetBranchAddress("time",&time);
  g_vdm_Tree->SetBranchAddress("nsep",&nsep);

  // prepare output tree with general info
  Int_t idx_separation_start;
  Int_t idx_separation_end;  
  Long_t time_separation_start;
  Long_t time_separation_end;
  TTree *sep_info_tree = new TTree("SepInfo","SepInfo");
  sep_info_tree->Branch("idx_separation_start",&idx_separation_start,"idx_separation_start/I");
  sep_info_tree->Branch("time_separation_start",&time_separation_start,"time_separation_start/L");
  sep_info_tree->Branch("idx_separation_end",&idx_separation_end,"idx_separation_end/I");
  sep_info_tree->Branch("time_separation_end",&time_separation_end,"time_separation_end/L");

  // get separations
  Double_t separation = 0;
  Double_t nsep_old = 0;
  Double_t small = 1e-12;
  Int_t counter = 0; // internal use only
  Double_t *sep_array = new Double_t[100]; // internal use only
  Double_t slope_change = 0;  // internal use only
  for(Int_t j=IdxStart;j<IdxEnd;j++) {
    // cout << " idx = " << IdxStart << " -- " << IdxEnd << endl;
    g_vdm_Tree->GetEntry(j);
    //  cout << " jgc  " << j << " d " << TMath::Abs(nsep-nsep_old) << endl;
    if (TMath::Abs(nsep-nsep_old)>small) {
      // new separation 
      idx_separation_start = j;
      time_separation_start = (Long_t) time;
      separation = nsep;
      sep_array[counter] = nsep;
      counter++;
      nsep_old = nsep;
      while(TMath::Abs(nsep-nsep_old)<small && j<IdxEnd) {
	j++;g_vdm_Tree->GetEntry(j);
      }
      j--;g_vdm_Tree->GetEntry(j);      
      idx_separation_end = j;
      time_separation_end = (Long_t) time;
      // slopes
      if (counter>2) {
	slope_change = (sep_array[counter-1]-sep_array[counter-2])*(sep_array[counter-2]-sep_array[counter-3]);
	//	cout << " slope_change  " << slope_change << "  c " << counter << endl;
      }
      if ((counter<3) || (slope_change>0))  {
	sep_info_tree->Fill();
	//  cout << scan_type<< " New separation = " << nsep << " starts at " << time_separation_start
	// << " ends at " << time_separation_end << " counter " << counter << endl;
      }
    }
    nsep_old = nsep;
  }

  // create file name
  char *file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/NomSep_x_Scan_%d.root",g_vdm_Fill,scan_num);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/NomSep_y_Scan_%d.root",g_vdm_Fill,scan_num);
  TFile *ScanFile = new TFile(file_name,"recreate");
  // create tree with separations
  TTree *sep_tree = new TTree("Separations","Separations");
  char *txt_tmp = new char[kg_string_size];
  sprintf(txt_tmp,"separation[%d]/D",counter-1);
  sep_tree->Branch("separation",sep_array,txt_tmp);

  // fill the same info for each bunch crossing
  // (the info is repeated, because it for nominal separation it is the same for all bunches)
  // -- number of bc
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  ScanFile->cd();
  for(Int_t k=0;k<nIBC;k++) sep_tree->Fill();
  sep_tree->Write();
  // fill separation info
  sep_info_tree->SetDirectory(ScanFile);
  sep_info_tree->Write();
  ScanFile->Close();

  // clean up
  delete [] file_name;
  delete [] txt_tmp;  
  delete [] sep_array;

  
}

//-------------------------------------------------------
// Create root files with the information of
// the nominal separations
//-------------------------------------------------------

void Create_nominal_separation_file(Int_t Fill)
{

  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // find indices for start and end of scans
  Find_start_and_end_of_scans();

  // create nominal separation files
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    //x-scans
    Find_separations(1,i,g_Idx_Start_Scan_x[i],g_Idx_End_Scan_x[i]);
    //y-scans
    Find_separations(2,i,g_Idx_Start_Scan_y[i],g_Idx_End_Scan_y[i]);
  }

}
