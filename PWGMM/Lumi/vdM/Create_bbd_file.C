
#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"

//-------------------------------------------------------
// these constants have to be modified for each Fill ...
//-------------------------------------------------------

Double_t kg_beta_x = 0.0;
Double_t kg_beta_y = 0.0;  
Double_t kg_Qx = 0.0;
Double_t kg_Qy = 0.0;
Double_t kg_Energy = 0.0;

void Init_parameters(Int_t Fill)
{
  if (Fill == 2852) {
    kg_beta_x = kg_beta_y = 19.0; // [m]
    kg_Qx = 64.31;
    kg_Qy = 59.32;
    kg_Energy = 4000.0; // [GeV]
  } else if (Fill == 4269) {
    kg_beta_x = kg_beta_y = 19.0; // [m]
    kg_Qx = 64.31;
    kg_Qy = 59.32;
    kg_Energy = 6500.0; // [GeV]
  } else if (Fill == 4634) {
    kg_beta_x = kg_beta_y = 10.0; // [m]
    kg_Qx = 64.31;
    kg_Qy = 59.32;
    kg_Energy = 2510.0; // [GeV]
  } else if (Fill == 4690) {
    kg_beta_x = kg_beta_y = 0.8; // [m]
    kg_Qx = 64.31;
    kg_Qy = 59.32;
    kg_Energy = 6369.0; // [GeV]
  } else if (Fill == 4937) {
    kg_beta_x = kg_beta_y = 18.91; // [m]
    kg_Qx = 64.31;
    kg_Qy = 59.32;
    kg_Energy = 6500.0; // [GeV]
  } else if (Fill == 5533 || Fill == 5568) {
    kg_beta_x = kg_beta_y = 2.0; // [m]
    kg_Qx = 64.31;
    kg_Qy = 59.32;
    kg_Energy = 6499.0; // [GeV]
  } else if (Fill == 7440 || Fill == 7483) {
    kg_beta_x = kg_beta_y = 0.51; // [m]
    kg_Qx = 64.31;
    kg_Qy = 59.32;
    kg_Energy = 6369.0; // [GeV]
  }
}

//-------------------------------------------------------
// Produces the files used by Ivan to compute the
// beam-beam deflection correction files
//-------------------------------------------------------

void Create_bbd_file(Int_t Fill, Int_t scan, const char *rate_name,
		     const char *rate_type,
		     const char *sep_type,
		     const char *intensity_type)
{
  // initialize
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  Init_parameters(Fill);

  // create name of output file
  char *file_name = new char[kg_string_size];
  sprintf(file_name,"../Fill-%d/BBD_%s_%s_%s_scan_%d.txt",
	  Fill,rate_name,intensity_type,sep_type,(scan+1));
  ofstream bbd_file;
  bbd_file.open(file_name);

  //-------------------------------------------------------
  // -- header part
  //-------------------------------------------------------
  bbd_file << "1) Fill, Scan :"<<Fill<<"," << (scan+1) << endl;
  bbd_file << "2) Betax*, Betay* [m,m] :"<<kg_beta_x<<","<<kg_beta_y << endl;
  bbd_file << "3) tunes Qx, Qy :"<<kg_Qx<<","<<kg_Qy << endl;
  bbd_file << "4) Energy [GeV] :"<< kg_Energy << endl;  
  //-------------------------------------------------------
  // -- separation part
  //-------------------------------------------------------
  // get separations in x
  char *sepx_file_name = new char[kg_string_size];
  sprintf(sepx_file_name,"../Fill-%d/%sSep_x_Scan_%d.root",g_vdm_Fill,sep_type,scan);
  TFile *sepx_file = new TFile(sepx_file_name);
  TTree *sepx_tree = (TTree *) sepx_file->Get("Separations");
  Int_t n_sepx = FindNumberSeparations(1, scan);
  Double_t *sepx = new Double_t[n_sepx];
  sepx_tree->ResetBranchAddresses();
  sepx_tree->SetBranchAddress("separation",sepx);
  sepx_tree->GetEntry(0);

  // get separations in y
  char *sepy_file_name = new char[kg_string_size];
  sprintf(sepy_file_name,"../Fill-%d/%sSep_y_Scan_%d.root",g_vdm_Fill,sep_type,scan);
  TFile *sepy_file = new TFile(sepy_file_name);
  TTree *sepy_tree = (TTree *) sepy_file->Get("Separations");
  Int_t n_sepy = FindNumberSeparations(1, scan);
  Double_t *sepy = new Double_t[n_sepy];
  sepy_tree->ResetBranchAddresses();
  sepy_tree->SetBranchAddress("separation",sepy);
  sepy_tree->GetEntry(0);

  // write out separations
  for(Int_t i=0;i<n_sepx;i++) 
    bbd_file << "5) Separation sepx, sepy [mm,mm]:" << sepx[i] << ",0" << endl;
  bbd_file << "5) Separation sepx, sepy [mm,mm] :0,0" << endl;
  bbd_file << "5) Separation sepx, sepy [mm,mm] :0,0" << endl;
  for(Int_t i=0;i<n_sepy;i++) 
    bbd_file << "5) Separation sepx, sepy [mm,mm]:0," << sepy[i]  << endl;


  //-------------------------------------------------------
  // -- open hx, hy files
  //-------------------------------------------------------
  // build names
  char *hx_file_name = new char[kg_string_size];
  sprintf(hx_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_x_Scan_%d_Fit_NUM.root",
	  g_vdm_Fill,rate_type,rate_name,sep_type,scan);
  char *hy_file_name = new char[kg_string_size];
  sprintf(hy_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_y_Scan_%d_Fit_NUM.root",
	  g_vdm_Fill,rate_type,rate_name,sep_type,scan);
  // open files and get trees
  TFile *hx_file = new TFile(hx_file_name);
  TFile *hy_file = new TFile(hy_file_name);
  TTree *hx_tree = (TTree *) hx_file->Get("AreaRate");
  TTree *hy_tree = (TTree *) hy_file->Get("AreaRate");  
  // set branches
  Double_t chi2_dof_x = -1;
  Double_t chi2_dof_y = -1;  
  Double_t *area_x = new Double_t[2]; // area and its error
  Double_t *rate_zero_x = new Double_t[2]; // rate at zero and its error
  Double_t *area_y = new Double_t[2]; // area and its error
  Double_t *rate_zero_y = new Double_t[2]; // rate at zero and its error
  // --> set branches for hx, hy
  hx_tree->ResetBranchAddresses();
  hy_tree->ResetBranchAddresses();
  hx_tree->SetBranchAddress("area",area_x);
  hy_tree->SetBranchAddress("area",area_y);
  hx_tree->SetBranchAddress("rate_zero",rate_zero_x);
  hy_tree->SetBranchAddress("rate_zero",rate_zero_y);
  hx_tree->SetBranchAddress("chi2_dof",&chi2_dof_x);
  hy_tree->SetBranchAddress("chi2_dof",&chi2_dof_y);  

  //-------------------------------------------------------
  // -- open intensity file
  //-------------------------------------------------------
  char *intensity_file_name = new char[kg_string_size];
  sprintf(intensity_file_name,"../Fill-%d/%s_Scan_%d.root",g_vdm_Fill,intensity_type,scan);
  TFile *intensity_file = new TFile(intensity_file_name);
  TTree *intensity_tree = (TTree *) intensity_file->Get("Bunch_Intensity");
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  Int_t *bunches = new Int_t [nIBC];
  GetBunchIndices(bunches);

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

  //-------------------------------------------------------
  // -- filling bc information
  //-------------------------------------------------------
  for(Int_t k=0;k<nIBC;k++) { // loop over bunches
    hx_tree->GetEntry(k);
    hy_tree->GetEntry(k);
    if(chi2_dof_x != 1 && chi2_dof_y != -1) {
      cout << " chi2 is not one: There is a problem somewhere! " << endl;
      exit(-1);
    }
    if (rate_zero_x[0] == 0 || rate_zero_y[0] == 0) {
      cout << " Rate is zero: There is a problem somewhere! " << endl;
      exit(-1);
    }
    bbd_file << "6) BC Sx,Sy,Np1,Np2 [um,um,prot,prot] :"
	     << bunches[k] << ","
	     << (1000*area_x[0]/(rate_zero_x[0]*TMath::Sqrt(TMath::TwoPi()))) << ","
	     << (1000*area_y[0]/(rate_zero_y[0]*TMath::Sqrt(TMath::TwoPi()))) << ","
	     << (cf_dcct_1*bunch_intensity_1[k]) << ","
	       << (cf_dcct_2*bunch_intensity_2[k]) << endl;
  }

  // last info in file
  bbd_file << "0) : 0" << endl;  
  bbd_file << "X) : 0" << endl;  
 
  //-------------------------------------------------------
  // -- closing/cleaning
  //-------------------------------------------------------
  
  // close file
  bbd_file.close();
  
  // clean-up
  delete [] file_name;
  delete [] sepx_file_name;
  delete [] sepy_file_name;
  delete [] sepx;
  delete [] sepy;  
  delete [] hx_file_name;
  delete [] hy_file_name;
  delete [] area_x;
  delete [] area_y;
  delete [] rate_zero_x;
  delete [] rate_zero_y;
  delete [] intensity_file_name;
  delete [] bunches;
  delete [] bunch_intensity_1;
  delete [] bunch_intensity_2;  

}

