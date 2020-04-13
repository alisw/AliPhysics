#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// get index of a given unix time
//-------------------------------------------------------

Int_t get_index(Long_t unix_time)
{
  // set time branch
  Double_t time;
  g_vdm_Tree->ResetBranchAddresses();
  g_vdm_Tree->SetBranchAddress("time",&time);

  // get time of first index
  g_vdm_Tree->GetEntry(0);
  Long_t time0 = ((Long_t) time);

  // get first approximation of index of requested time
  // in the tree the unix time is stored every 2 seconds, but on ocassions one entry is skipped
  Int_t diff = ((Int_t) (unix_time -  time0));
  Int_t idx = 0;
  (diff%2? idx = diff>>1 : idx = (diff-1)>>1);
  g_vdm_Tree->GetEntry(idx);
  Long_t time_idx = ((Long_t) time);
  Int_t diff_first = ((Int_t) (unix_time -  time_idx));
  if (TMath::Abs(diff_first)>1) {
    idx = idx + (diff_first>>1);
    g_vdm_Tree->GetEntry(idx);
    time_idx = ((Long_t) time);
  }

  return idx;
}

//-------------------------------------------------------
// get info on first X scan
// (info from file provided by Christoph
//-------------------------------------------------------

void GetInfoFirstXScan(Int_t *idx_start, Int_t *idx_end,
		       Long_t *time_start, Long_t *time_end,
		       Double_t *separation)
{
  Long_t t_start[25] = {
    1543480289,
    1543480328,
    1543480366,
    1543480404,
    1543480442,
    1543480480,
    1543480518,
    1543480556,
    1543480594,
    1543480633,
    1543480671,
    1543480709,
    1543480747,
    1543480785,
    1543480824,
    1543480862,
    1543480900,
    1543480938,
    1543480976,
    1543481014,
    1543481053,
    1543481091,
    1543481129,
    1543481167,
    1543481205
  };
  Long_t t_end[25] = {
    1543480319,
    1543480358,
    1543480396,
    1543480434,
    1543480472,
    1543480510,
    1543480549,
    1543480586,
    1543480624,
    1543480663,
    1543480701,
    1543480739,
    1543480777,
    1543480815,
    1543480854,
    1543480892,
    1543480930,
    1543480967,
    1543481006,
    1543481044,
    1543481083,
    1543481121,
    1543481159,
    1543481197,
    1543481235
  };
  Double_t sep[25] = {
    -0.09729876,
    -0.08919053,
    -0.08108230,
    -0.07297407,
    -0.06486584,
    -0.05675761,
    -0.04864938,
    -0.04054115,
    -0.03243292,
    -0.02432469,
    -0.01621646,
    -0.00810823,
    0.00000000,
    0.00810823,
    0.01621646,
    0.02432469,
    0.03243292,
    0.04054115,
    0.04864938,
    0.05675761,
    0.06486584,
    0.07297407,
    0.08108230,
    0.08919053,
    0.09729876
  };
  // now fill the info
  for(Int_t i=0;i<25;i++) {
    time_start[i] = t_start[i];
    time_end[i] = t_end[i];
    separation[i] = sep[i];
    idx_start[i] = get_index(t_start[i]);
    idx_end[i] = get_index(t_end[i]);    
  }
}

//-------------------------------------------------------
// get info on first Y scan
// (info from file provided by Christoph
//-------------------------------------------------------

void GetInfoFirstYScan(Int_t *idx_start, Int_t *idx_end,
		       Long_t *time_start, Long_t *time_end,
		       Double_t *separation)
{
  Long_t t_start[25] = {
    1543482080,
    1543482118,
    1543482156,
    1543482195,
    1543482233,
    1543482271,
    1543482310,
    1543482348,
    1543482386,
    1543482424,
    1543482462,
    1543482500,
    1543482538,
    1543482577,
    1543482615,
    1543482653,
    1543482691,
    1543482729,
    1543482767,
    1543482805,
    1543482844,
    1543482882,
    1543482920,
    1543482958,
    1543482996
   };
  Long_t t_end[25] = {
    1543482110,
    1543482148,
    1543482186,
    1543482225,
    1543482263,
    1543482301,
    1543482340,
    1543482378,
    1543482416,
    1543482454,
    1543482492,
    1543482530,
    1543482568,
    1543482607,
    1543482645,
    1543482683,
    1543482721,
    1543482759,
    1543482797,
    1543482835,
    1543482874,
    1543482912,
    1543482950,
    1543482988,
    1543483026
  };
  Double_t sep[25] = {
    -0.09729876,
    -0.08919053,
    -0.08108230,
    -0.07297407,
    -0.06486584,
    -0.05675761,
    -0.04864938,
    -0.04054115,
    -0.03243292,
    -0.02432469,
    -0.01621646,
    -0.00810823,
    0.00000000,
    0.00810823,
    0.01621646,
    0.02432469,
    0.03243292,
    0.04054115,
    0.04864938,
    0.05675761,
    0.06486584,
    0.07297407,
    0.08108230,
    0.08919053,
    0.09729876
  };
  // now fill the info
  for(Int_t i=0;i<25;i++) {
    time_start[i] = t_start[i];
    time_end[i] = t_end[i];
    separation[i] = sep[i];
    idx_start[i] = get_index(t_start[i]);
    idx_end[i] = get_index(t_end[i]);    
  }
}

//-------------------------------------------------------
// get info on second X scan
// (info from file provided by Christoph
//-------------------------------------------------------

void GetInfoSecondXScan(Int_t *idx_start, Int_t *idx_end,
			Long_t *time_start, Long_t *time_end,
			Double_t *separation)
{
  Long_t t_start[25] = {
    1543486621,
    1543486659,
    1543486697,
    1543486735,
    1543486773,
    1543486812,
    1543486850,
    1543486888,
    1543486926,
    1543486964,
    1543487002,
    1543487040,
    1543487078,
    1543487116,
    1543487155,
    1543487193,
    1543487231,
    1543487269,
    1543487307,
    1543487346,
    1543487384,
    1543487422,
    1543487460,
    1543487498,
    1543487536
   };
  Long_t t_end[25] = {
    1543486651,
    1543486689,
    1543486727,
    1543486765,
    1543486803,
    1543486842,
    1543486880,
    1543486918,
    1543486956,
    1543486994,
    1543487032,
    1543487070,
    1543487108,
    1543487146,
    1543487185,
    1543487223,
    1543487261,
    1543487299,
    1543487337,
    1543487376,
    1543487414,
    1543487452,
    1543487490,
    1543487528,
    1543487566
  };
  Double_t sep[25] = {
    -0.09729876,
    -0.08919053,
    -0.08108230,
    -0.07297407,
    -0.06486584,
    -0.05675761,
    -0.04864938,
    -0.04054115,
    -0.03243292,
    -0.02432469,
    -0.01621646,
    -0.00810823,
    0.00000000,
    0.00810823,
    0.01621646,
    0.02432469,
    0.03243292,
    0.04054115,
    0.04864938,
    0.05675761,
    0.06486584,
    0.07297407,
    0.08108230,
    0.08919053,
    0.09729876
  };
  // now fill the info
  for(Int_t i=0;i<25;i++) {
    time_start[i] = t_start[i];
    time_end[i] = t_end[i];
    separation[i] = sep[i];
    idx_start[i] = get_index(t_start[i]);
    idx_end[i] = get_index(t_end[i]);    
  }
}
//-------------------------------------------------------
// get info on second Y scan
// (info from file provided by Christoph
//-------------------------------------------------------

void GetInfoSecondYScan(Int_t *idx_start, Int_t *idx_end,
			Long_t *time_start, Long_t *time_end,
			Double_t *separation)
{
  Long_t t_start[25] = {
    1543487708,
    1543487746,
    1543487784,
    1543487822,
    1543487860,
    1543487898,
    1543487936,
    1543487974,
    1543488013,
    1543488051,
    1543488090,
    1543488128,
    1543488166,
    1543488204,
    1543488242,
    1543488280,
    1543488319,
    1543488357,
    1543488395,
    1543488433,
    1543488472,
    1543488510,
    1543488548,
    1543488586,
    1543488624
  };
  Long_t t_end[25] = {
    1543487737,
    1543487776,
    1543487814,
    1543487852,
    1543487890,
    1543487928,
    1543487966,
    1543488004,
    1543488043,
    1543488082,
    1543488120,
    1543488158,
    1543488196,
    1543488234,
    1543488272,
    1543488310,
    1543488349,
    1543488387,
    1543488425,
    1543488463,
    1543488502,
    1543488540,
    1543488578,
    1543488616,
    1543488654
  };
  Double_t sep[25] = {
    -0.09729876,
    -0.08919053,
    -0.08108230,
    -0.07297407,
    -0.06486584,
    -0.05675761,
    -0.04864938,
    -0.04054115,
    -0.03243292,
    -0.02432469,
    -0.01621646,
    -0.00810823,
    0.00000000,
    0.00810823,
    0.01621646,
    0.02432469,
    0.03243292,
    0.04054115,
    0.04864938,
    0.05675761,
    0.06486584,
    0.07297407,
    0.08108230,
    0.08919053,
    0.09729876
  };
  // now fill the info
  for(Int_t i=0;i<25;i++) {
    time_start[i] = t_start[i];
    time_end[i] = t_end[i];
    separation[i] = sep[i];
    idx_start[i] = get_index(t_start[i]);
    idx_end[i] = get_index(t_end[i]);    
    
  }
}


//-------------------------------------------------------
// Create the nominal separation file for a given scan
//-------------------------------------------------------
void CreateFile(Int_t scan_type, Int_t scan_num,
		Int_t *idx_start, Int_t *idx_end,
		Long_t *time_start, Long_t *time_end,
		Double_t *sep)
// scan_type: 1 => x-scan; 2 => y-scan
{

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

  for(Int_t i=0;i<25;i++) {
    idx_separation_start = idx_start[i];
    idx_separation_end = idx_end[i];    
    time_separation_start = time_start[i];
    time_separation_end = time_end[i];
    sep_info_tree->Fill();
  }
  
  // create file name
  char *file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(file_name,"../Fill-7483/NomSep_x_Scan_%d.root",scan_num);
  if (scan_type == 2) sprintf(file_name,"../Fill-7483/NomSep_y_Scan_%d.root",scan_num);
  TFile *ScanFile = new TFile(file_name,"recreate");
  // create tree with separations
  TTree *sep_tree = new TTree("Separations","Separations");
  sep_tree->Branch("separation",sep,"separation[25]/D");

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
  
}


//-------------------------------------------------------
// create the nominal separation file for fill 7483
// this case is special because (a) data taking was stopped
// and restarted and (b) there is a diagonal scan in the middle
// data used was provided by Christoph
//-------------------------------------------------------

void Create_nominal_separation_file_Fill_7483()
{
  // get name of files and set pointers to trees
  Set_input_file_names(7483);
  Set_pointers_to_input_files_and_trees();

  // reserve space for info
  Int_t idx_start[25];
  Int_t idx_end[25];
  Long_t time_start[25];
  Long_t time_end[25];
  Double_t separation[25];

  // get info first X scan
  GetInfoFirstXScan(idx_start, idx_end,time_start,time_end,separation);
  CreateFile(1,0,idx_start, idx_end,time_start,time_end,separation);
  // get info first Y scan
  GetInfoFirstYScan(idx_start, idx_end,time_start,time_end,separation);
  CreateFile(2,0,idx_start, idx_end,time_start,time_end,separation);
  // get info second X scan
  GetInfoSecondXScan(idx_start, idx_end,time_start,time_end,separation);
  CreateFile(1,1,idx_start, idx_end,time_start,time_end,separation);
  // get info second Y scan
  GetInfoSecondYScan(idx_start, idx_end,time_start,time_end,separation);
  CreateFile(2,1,idx_start, idx_end,time_start,time_end,separation);

}

