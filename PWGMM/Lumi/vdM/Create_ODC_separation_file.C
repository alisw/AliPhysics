

#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//------------------------------------------------------------
// Evaluate the orbit-drift correction for fill 4690
//------------------------------------------------------------

Double_t Get_ODC_x_scan_4690(Double_t mid_time)
{
  Double_t   p0 = 1449192576;
  Double_t   p1 = 0.4323;
  Double_t   p2 = 5.296e-5;
  Double_t x = (mid_time-p0);
  // divide by 1000 to get right units
  return (p1+p2*x)/1000; 
}

Double_t Get_ODC_y_scan_4690(Double_t mid_time)
{
  Double_t   p0 = 1449192576;
  Double_t   p1 = 0.6993;
  Double_t   p2 = 2.273e-5;
  Double_t x = (mid_time-p0);
  // divide by 1000 to get right units
  return (p1+p2*x)/1000; 
}

//------------------------------------------------------------
// Evaluate the orbit-drift correction for fill 4937
//------------------------------------------------------------

Double_t Get_ODC_x_scan_4937(Double_t mid_time)
{
  Double_t   p0 = 1463505614;
  Double_t   p1 = 0.197442;
  Double_t   p2 = 0.00156;
  Double_t   p3 = 6.50027e-7;
  Double_t x = (mid_time-p0);
  // divide by 1000 to get right units
  return (p1+p2*x+p3*x*x)/1000; 
}

Double_t Get_ODC_y_scan_4937(Double_t mid_time)
{
  Double_t   p0 = 1463505614;
  Double_t   p1 = -0.464639;
  Double_t   p3 = 9.45379e-8;
  Double_t   p4 = 2.13991e-11;
  Double_t x = (mid_time-p0);
  // divide by 1000 to get right units
  return (p1+p3*x*x+p4*x*x*x)/1000; 
}


//------------------------------------------------------------
// Evaluate the orbit-drift correction for fill 5533
//------------------------------------------------------------

Double_t Get_ODC_x_scan_5533(Double_t mid_time)
{
  Double_t   p0 = 1479918502;
  Double_t   p1 = 2.0641;
  Double_t   p2 = 0.00161356;
  Double_t   p3 = -2.1392e-8;
  Double_t   p4 = -1.12519e-10;
  Double_t x = (mid_time-p0);
  // divide by 1000 to get right units
  return (p1+p2*x+p3*x*x+p4*x*x*x)/1000; 
}

Double_t Get_ODC_y_scan_5533(Double_t mid_time)
{
  Double_t   p0 = 1479918502;
  Double_t   p1 = -0.30024;
  Double_t   p2 = 0.000314952;
  Double_t   p3 = 1.34475e-7;
  Double_t x = (mid_time-p0);
  // divide by 1000 to get right units
  return (p1+p2*x+p3*x*x)/1000; 
}

//------------------------------------------------------------
// Evaluate the orbit-drift correction for fill 5568
//------------------------------------------------------------

Double_t Get_ODC_x_scan_5568(Double_t mid_time)
{
  Double_t   p0 = 1480695439;
  Double_t   p1 = -2.80692;
  Double_t   p2 = -0.000869452;
  Double_t x = (mid_time-p0);
  // divide by 1000 to get right units
  return (p1+p2*x)/1000; 
}

Double_t Get_ODC_y_scan_5568(Double_t mid_time)
{
  Double_t   p0 = 1480695439;
  Double_t   p1 = 3.13328;
  Double_t   p2 = 0.00094365;
  Double_t x = (mid_time-p0);
  // divide by 1000 to get right units
  return (p1+p2*x)/1000; 
}

//------------------------------------------------------------
// Evaluate the orbit-drift correction for fill 7483
//------------------------------------------------------------

Double_t Get_ODC_x_scan0_7483(Double_t mid_time)
{
 
  Double_t   p0 = -1863070;
  Double_t   p1 = 0.00120706;
  // divide by 1000 to get right units
  return ((p0+p1*mid_time)/1000); 
}

Double_t Get_ODC_y_scan0_7483(Double_t mid_time)
{
  Double_t   p0 = 791330;
  Double_t   p1 = -0.000512693;
  // divide by 1000 to get right units
  return ((p0+p1*mid_time)/1000); 
}

Double_t Get_ODC_x_scan1_7483(Double_t mid_time)
{
  Double_t   p0 = 109870;
  Double_t   p1 = -7.11833e-5;
  // divide by 1000 to get right units
  return ((p0+p1*mid_time)/1000); 
}

Double_t Get_ODC_y_scan1_7483(Double_t mid_time)
{
  Double_t   p0 = -183147;
  Double_t   p1 = 0.000118658;
  // divide by 1000 to get right units
  return ((p0+p1*mid_time)/1000); 
}

//------------------------------------------------------------
// Evaluate the orbit-drift correction for the different fills
//------------------------------------------------------------

Double_t Get_ODC_separations_Fill(Double_t mid_time, Int_t scan_type, Int_t scan)
// scan_type: 1 => x-scan; 2 => y-scan
{
  if (g_vdm_Fill == 4690) {
    if (scan_type == 1) return Get_ODC_x_scan_4690(mid_time);    
    if (scan_type == 2) return Get_ODC_y_scan_4690(mid_time);
    return 0;
  } else if (g_vdm_Fill == 4937) {
    if (scan_type == 1) return Get_ODC_x_scan_4937(mid_time);    
    if (scan_type == 2) return Get_ODC_y_scan_4937(mid_time);
    return 0;
  } else if (g_vdm_Fill == 5533) {
    if (scan_type == 1) return Get_ODC_x_scan_5533(mid_time);    
    if (scan_type == 2) return Get_ODC_y_scan_5533(mid_time);
    return 0;
  } else if (g_vdm_Fill == 5568) {
    if (scan_type == 1) return Get_ODC_x_scan_5568(mid_time);    
    if (scan_type == 2) return Get_ODC_y_scan_5568(mid_time);
    return 0;
  } else if (g_vdm_Fill == 7483) {
    if (scan_type == 1 && scan == 0) return Get_ODC_x_scan0_7483(mid_time);
    if (scan_type == 1 && scan == 1) return Get_ODC_x_scan1_7483(mid_time);
    if (scan_type == 2 && scan == 0) return Get_ODC_y_scan0_7483(mid_time);
    if (scan_type == 2 && scan == 1) return Get_ODC_y_scan1_7483(mid_time);
    return 0;
  } else {
    cout << " Fill not implemented: returning zero " << endl;
    return 0;
  }
}

//-------------------------------------------------------
// Correct the separations for orbit drift using the
// fit of orbit-drift data
//-------------------------------------------------------

void Get_ODC_separations(Int_t scan_type, Int_t scan)
// scan_type: 1 => x-scan; 2 => y-scan
{
  // --> create nominal separation file name
  char *sep_file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(sep_file_name,"../Fill-%d/NomSep_x_Scan_%d.root",g_vdm_Fill,scan);
  if (scan_type == 2) sprintf(sep_file_name,"../Fill-%d/NomSep_y_Scan_%d.root",g_vdm_Fill,scan);

  // open separation file and get the tree
  TFile *sep_file = new TFile(sep_file_name);
  TTree *sep_tree = (TTree *) sep_file->Get("Separations");
  TTree *sepinfo_tree = (TTree *) sep_file->Get("SepInfo");  
  
  // set branches
  // -- get number of separations
  Int_t n_sep = FindNumberSeparations(scan_type, scan);
  Double_t *sep = new Double_t[n_sep];
  sep_tree->ResetBranchAddresses();
  sep_tree->SetBranchAddress("separation",sep);
  Long64_t time_start;
  Long64_t time_end;  
  sepinfo_tree->ResetBranchAddresses();
  sepinfo_tree->SetBranchAddress("time_separation_start",&time_start);
  sepinfo_tree->SetBranchAddress("time_separation_end",&time_end);  

  // compute ODC (the same for all BC)
  Double_t *odc_sep = new Double_t[n_sep];
  sep_tree->GetEntry(0);
  for(Int_t i=0;i<n_sep;i++) {
    sepinfo_tree->GetEntry(i);
    Double_t mid_time = (Double_t) (time_end+time_start);
    mid_time *= 0.5;
    Double_t corr = Get_ODC_separations_Fill(mid_time,scan_type,scan);
    odc_sep[i] = sep[i]+corr;
    // cout << " i " << i << " end " << time_end << " start " << time_start << endl;
    //   cout << " mid_time " << ((Long_t) mid_time) << " corr " << corr << " odc_sep " << odc_sep[i] << " sep " << sep[i] << endl;
  }
  
  // now prepare name of output and reserve space to store it
  if (scan_type == 1) sprintf(sep_file_name,"../Fill-%d/ODCSep_x_Scan_%d.root",
			      g_vdm_Fill,scan);
  if (scan_type == 2) sprintf(sep_file_name,"../Fill-%d/ODCSep_y_Scan_%d.root",
			      g_vdm_Fill,scan);
  TFile *ODCFile = new TFile(sep_file_name,"recreate");
  // -- create tree with separations
  TTree *odc_sep_tree = new TTree("Separations","Separations");
  char *txt_tmp = new char[kg_string_size];
  sprintf(txt_tmp,"separation[%d]/D",n_sep);
  odc_sep_tree->Branch("separation",odc_sep,txt_tmp);

  // fill the same info for each bunch crossing
  // (the info is repeated, because it for ODC separation it is the same for all bunches)
  // -- number of bc
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  ODCFile->cd();
  for(Int_t k=0;k<nIBC;k++) odc_sep_tree->Fill();
  odc_sep_tree->Write();
  ODCFile->Close();

  // clean up
  delete [] sep_file_name;
  delete [] sep;
  delete [] odc_sep;  
  delete [] txt_tmp;  

}


//-------------------------------------------------------
// Create root files with the information of
// the orbit-drift corrected separations
//-------------------------------------------------------

void Create_ODC_separation_file(Int_t Fill)
{

  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // find indices for start and end of scans
  Find_start_and_end_of_scans();

  // create orbit-drift corrected separation files
   for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    //x-scans
    Get_ODC_separations(1,i);
    //y-scans
    Get_ODC_separations(2,i);
   }

}
